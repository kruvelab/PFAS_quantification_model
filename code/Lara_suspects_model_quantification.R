
# Change the admin directory
admin <- "C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant"

setwd(admin)
source("code/functions.R")
library(caTools)
library(tidyverse)
library(caret)
library(plotly)
library(cowplot)


#---------------------
# Testing Lara's data
#---------------------

#### When calibrants and suspects are in the same TraceFinder file; same for SMILES ####
setwd(admin)
logIE_pred_model <- readRDS(paste0(admin,"/models/model_PFAS_logIE.Rdata"))


lara_concentrations_pred <- concentration_forAnalytes_model(filename_data = paste0(admin,"/Lara_suspects_model_quantification/20220208_Suspect_screening_TU pools.xlsx"),
                                                            filename_smiles = paste0(admin,"/Lara_suspects_model_quantification/Smiles_for_Target_PFAS_semicolon_lara.csv"),
                                                            filename_eluent = "data/eluent.csv",
                                                            pred_model =  logIE_pred_model,
                                                            compounds_to_be_removed_as_list = c(
                                                              "HFPO-DA", "MeFOSE", "EtFOSE", "10:2 mono PAP", "4:2 mono PAP", "6:2 mono PAP", "8:2 mono PAP")) 



#### When calibrants and suspects are in separate TraceFinder and SMILES files ####
setwd(admin)
cal_filename_data <-  paste0(admin,"/Lara_suspects_model_quantification/Target PFAS for semiquant_120722.xlsx")
cal_filename_smiles <- paste0(admin,"/Lara_suspects_model_quantification/Target_PFAS_Calibrants.csv")
sus_filename_data <- paste0(admin,"/Lara_suspects_model_quantification/20220208_Suspect_screening_TU pools.xlsx")
sus_filename_smiles <- paste0(admin,"/Lara_suspects_model_quantification/suspects_smiles_neutral.csv")
logIE_pred_model <- readRDS(paste0(admin,"/models/model_PFAS_logIE.Rdata"))

lara_concentrations_pred <- concentration_forAnalytes_model_cal_separateFile(cal_filename_data, 
                                                                             cal_filename_smiles, 
                                                                             sus_filename_data,
                                                                             sus_filename_smiles,
                                                                             filename_eluent = "data/eluent.csv",
                                                                             pred_model =  logIE_pred_model,
                                                                             compounds_to_be_removed_as_list = c("PFHpS-br", "PFPeS", "PFHpS", "PFNS", "PFPeDA", "10:2 mono PAP", "4:2 mono PAP", "6:2 mono PAP", "8:2 mono PAP"))


####################################################
#### Save out shorter version of concentrations ####

lara_pred <- lara_concentrations_pred$data

lara_pred <- lara_pred %>%
  select(-c(IC, Molecular_weight, area_IC)) %>%
  rename(Predicted_RF = slope_pred,
         Predicted_conc_uM = conc_pred)

#write_delim(lara_pred, paste0(admin,"/Lara_suspects_model_quantification/Lara_pred_conc_neutralSuspectsSMILES.csv"), delim = ",")



###################################################################
#### Plot calibration graph & predicted vs real concentrations ####

ggplotly(lara_concentrations_pred$plot_predicted_theoretical_conc)
ggplotly(lara_concentrations_pred$plot_predictedIE_slope)



###########################################
#### barplot of suspect concentrations ####

#bar plot
data_short = lara_pred %>%
  filter(Compound %in% c("d/C PFSA n=8", "eecec PFSA n=8", "ether PFSA n=4", "ether PFSA n=8", "H-PFDoDA", "H-PFDS", "NMe-FBSAA")) %>%
  #filter(Compound %in% c("NMe-FBSAA")) %>% 
  mutate(unique_compound = paste0(Compound, " (", SMILES, ")")) %>%
  filter(grepl("Pool", Filename, fixed = TRUE))

barplot <- ggplot(data = data_short, aes( x = factor( unique_compound ), y = ((Predicted_conc_uM)), fill = Compound ) ) +    # print bar chart
  geom_bar( stat = 'identity', position = 'dodge') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("") +
  ylab("Concentration (nM)")+
  my_theme +
  facet_wrap(~ Filename, ncol = 6)
  #coord_flip()

#ggsave(barplot, filename = paste0(admin,"/Lara_suspects_model_quantification/suspects_summary_byFilename.svg"), width=40, height=80, units = "cm")



