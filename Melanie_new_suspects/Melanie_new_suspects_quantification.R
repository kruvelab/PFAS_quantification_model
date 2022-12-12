
# Change the admin directory
admin <- "C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant"

setwd(admin)
source("code/functions.R")
library(caTools)
library(tidyverse)
library(caret)
library(plotly)
library(cowplot)



#------------------------------------------------
# Semi-quantification of new suspects (Melanie)
#------------------------------------------------

#### When calibrants and suspects are in separate TraceFinder and SMILES files ####
setwd(admin)
cal_filename_data <-  paste0(admin,"/data/Batch 1 Semi Quant w frag.xlsx")
cal_filename_smiles <- paste0(admin,"/data/Smiles_for_Target_PFAS_semicolon.csv")
sus_filename_data <- paste0(admin,"/Melanie_new_suspects/20210810_Melanie_Suspect _Screening_TF.xlsx")
sus_filename_smiles <- paste0(admin,"/Melanie_new_suspects/suspects_smiles_melanie_updated_semicolon.csv")
logIE_pred_model <- readRDS(paste0(admin,"/models/221205_model_PFAS_allData_logIE.Rdata"))

melanie_concentrations_pred <- concentration_forAnalytes_model_cal_separateFile(cal_filename_data,
                                                                             cal_filename_smiles,
                                                                             sus_filename_data,
                                                                             sus_filename_smiles,
                                                                             filename_eluent = "data/eluent.csv",
                                                                             pred_model =  logIE_pred_model,
                                                                             compounds_to_be_removed_as_list = c("HFPO-DA", "MeFOSE", "EtFOSE", "PFHpS-br", "PFPeS", "PFHpS", "PFNS", "PFPeDA", "10:2 mono PAP", "4:2 mono PAP", "6:2 mono PAP", "8:2 mono PAP"))


####################################################
#### Save out shorter version of concentrations ####

melanie_concentrations_pred_conc <- melanie_concentrations_pred$data

melanie_concentrations_pred_conc <- melanie_concentrations_pred_conc %>%
  select(-c(IC, Molecular_weight, area_IC)) %>%
  rename(Predicted_RF = slope_pred)

#write_delim(melanie_concentrations_pred_conc, paste0(admin,"/Melanie_new_suspects/221205_new_suspects_pred_conc.csv"), delim = ";")



###################################################################
#### Plot calibration graph & predicted vs real concentrations ####

# Calibration plots to all cal compounds - can check if something needs to be removed that does not have a nice cal graph
melanie_concentrations_pred$all_calibration_plots

# Model calibration - response factors and ionization efficiencies for calibration compounds
ggplotly(melanie_concentrations_pred$plot_predictedIE_slope)

# Model results based on calibration compounds - predicted vs real concentrations (uM)
ggplotly(melanie_concentrations_pred$plot_predicted_theoretical_conc)




###########################################
#### barplot of suspect concentrations ####

#bar plot
data_short = melanie_concentrations_pred_conc %>%
  #filter(Compound %in% c("d/C PFSA n=8", "eecec PFSA n=8", "ether PFSA n=4", "ether PFSA n=8", "H-PFDoDA", "H-PFDS", "NMe-FBSAA")) %>% #all suspects
  filter(Compound %in% c("ClPFLCAs n=8")) %>%
  mutate(unique_compound = paste0(Compound, " (", SMILES, ")"))


barplot <- ggplot(data = data_short, aes( x = factor( unique_compound ), y = ((conc_pred_pg_uL)), fill = Compound ) ) +    # print bar chart
  geom_bar( stat = 'identity', position = 'dodge') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("") +
  ylab("Concentration (pg/uL)")+
  facet_wrap(~ Filename, ncol = 6)


#ggsave(barplot, filename = paste0(admin,"/Melanie_new_suspects/compound_barplot_concentration.svg"), width=40, height=80, units = "cm")

