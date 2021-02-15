library("tidyverse")
source("PaDEL_descs_calculator.R")
source("compound_eluent.R")

#here set the folder where you have all thesefiles as working directory

regressor_neg = readRDS("regressor_neg_new.rds")
descs_names = readRDS("negative_descs.rds")

standards = read_delim("PCB_standards.csv",
                       delim = ",",
                       col_names = TRUE)


descriptor_calculated = PaDEL_original(standards)

descriptor_calculated = descriptor_calculated %>%
  select(SMILES, descs_names) %>% 
  mutate(
    organic_modifier = "MeCN",
    organic = 80,
    pH.aq. = 8.0,
    NH4 = 1, #1 if th ebuffer contains NH¤ ions , 0 if not. 
    viscosity =  viscosity(organic,organic_modifier),
    surface_tension = surfacetension(organic,organic_modifier),
    polarity_index = polarityindex(organic,organic_modifier)) 

IE_pred = descriptor_calculated %>% 
  mutate(logIE_pred = 0) %>%
  na.omit()
  
prediction =  predict(regressor_neg, newdata = IE_pred, predict.all = TRUE)
prediction = prediction$aggregate
IE_pred <- IE_pred %>%
  mutate(logIE_pred = prediction) %>%
  select(SMILES,logIE_pred)

write_delim(IE_pred,
            "Ionization_predictions.csv",
            delim = ",") # either ";" or ","

