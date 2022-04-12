#install.packages("tidyverse")
#install.packages("stringr")
#install.packages("rjson")
library(tidyverse)
library(stringr)
library(rjson)

PaDEL = function(SMILES_list) {
  descs = tibble() #empty data table of descriptors
  command = "java -jar descriptor-cli-0.1a-SNAPSHOT-all.jar" #file name where they will be calculated
  for (i in 1:dim(SMILES_list)[1]) {
    smiles = SMILES_list$SMILES[i] #takes the respective SMILES
    print(smiles) 
    command_final = paste(command, smiles, sep =" ") #makes text for command prompt
    javaOutput = system(command_final, intern = TRUE) #goes into commant prompt
    output_string = str_split(javaOutput, pattern = " ") #from here processing output
    descs_json = paste(output_string[1][[1]], output_string[2][[1]], output_string[3][[1]], output_string[4][[1]], sep = "")
    descs_this_smiles = fromJSON(descs_json)
    descs_this_smiles = data.frame(descs_this_smiles) #get descs in a table format
    descs_this_smiles = descs_this_smiles %>%
      mutate(SMILES = smiles) #adds compound structure to the descs
    descs = descs %>% #add descs for this SMILES to the collection of all of the SMILES
      bind_rows(descs_this_smiles)
  }
  return(descs)
}

PaDEL_original = function(standards) {
  SMILES_list = standards %>% 
    select(SMILES) %>% 
    unique() %>%
    na.omit() %>%
    mutate(Name = row_number())
  standards = standards %>%
    left_join(SMILES_list)
  
  write_delim(SMILES_list %>% select(SMILES) %>% unique(),
              "SMILES.smi",
              col_names = FALSE)
  
  command = "java -jar PaDEL-Descriptor/PaDEL-Descriptor.jar -dir" #file name where they will be calculated
  command_final = paste(command, "SMILES.smi", "-file", "descs_calc.csv", "-2d", sep =" ") #makes text for command prompt
  javaOutput = system(command_final, intern = TRUE) #goes into commant prompt
  #PaDEL saves the descriptors to the local folder
  descs = read_delim("descs_calc.csv",
                     delim = ",",
                     col_names = TRUE)
  
  descs = descs %>%
    group_by(Name) %>%
    mutate(Name = str_split(Name, pattern = "_")[[1]][length(str_split(Name, pattern = "_")[[1]])]) %>%
    ungroup() %>%
    mutate(Name = as.numeric(Name))
  
  descs = descs %>%
    left_join(standards) %>%
    select(colnames(standards), everything()) %>%
    select(-Name)
  
  write_delim(descs,
              "data/descs_calc.csv",
              delim = ",")
  
  return(descs)
}
