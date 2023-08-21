
# Change the admin directory
admin <- "C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant"
setwd(admin)
source(paste0(admin,"/02_code/functions.R"))
library(caTools)
library(tidyverse)
library(caret)
library(plotly)
library(cowplot)


#------------------------------------------------
# Semi-quantification of new suspects 
#------------------------------------------------

#### When calibrants and suspects are in separate TraceFinder and SMILES files ####
cal_filename_data <-  paste0(admin,"/05_suspect_quantification/Batch 1 Semi Quant w frag.xlsx")
cal_filename_smiles <- paste0(admin,"/05_suspect_quantification/Smiles_for_Target_PFAS_semicolon.csv")
sus_filename_data <- paste0(admin,"/05_suspect_quantification/20210810_Suspect_Screening_TF.xlsx")
sus_filename_smiles <- paste0(admin,"/05_suspect_quantification/suspects_smiles_updated_semicolon2.csv")
logIE_pred_model <- readRDS(paste0(admin,"/03_models/230619_logIE_model_withPFAS_allData.Rdata"))


concentrations_pred <- concentration_forAnalytes_model_cal_separateFile(cal_filename_data,
                                                                             cal_filename_smiles,
                                                                             sus_filename_data,
                                                                             sus_filename_smiles,
                                                                             filename_eluent = paste0(admin, "/05_suspect_quantification/eluent.csv"),
                                                                             pred_model =  logIE_pred_model,
                                                                             compounds_to_be_removed_as_list = c("HFPO-DA", "MeFOSE", "EtFOSE", "PFHpS-br", "PFPeS", "PFHpS", "PFNS", "PFPeDA", "10:2 mono PAP", "4:2 mono PAP", "6:2 mono PAP", "8:2 mono PAP"))


#################################
#### Save out concentrations ####

concentrations_pred_conc <- concentrations_pred$data

concentrations_pred_conc <- concentrations_pred_conc %>%
  select(-c(IC, Molecular_weight, area_IC)) %>%
  rename(Predicted_RF = slope_pred)


###################################################################
#### Plot calibration graph & predicted vs real concentrations ####

# Calibration plots to all cal compounds - can check if something needs to be removed that does not have a nice cal graph
concentrations_pred$all_calibration_plots

# Model calibration - response factors and ionization efficiencies for calibration compounds
ggplotly(concentrations_pred$plot_predictedIE_slope)

# Model results based on calibration compounds - predicted vs real concentrations (uM)
ggplotly(concentrations_pred$plot_predicted_theoretical_conc)

