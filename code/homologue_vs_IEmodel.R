
#-----------------------------------------------------------
# Homologous series vs ML model predictions - concentrations
#-----------------------------------------------------------

# find homolog series compounds
# save out actual concentrations
# find homolog series quantifications
# train model by leaving one out each time
# Comparison in table


# Find homolog series compounds
data = Orbitrap_dataset_raw %>%
  left_join(SMILES_data) %>%
  drop_na(SMILES)

SMILES_forHomolog <- data %>%
  select(Compound, SMILES) %>%
  unique()

SMILEs_match <- SMILES_forHomolog %>%
  left_join(SMILES_forHomolog %>%
              select(SMILES), by = character())

homologs <- SMILEs_match %>%
  group_by(SMILES.x, SMILES.y) %>%
  mutate(pattern_CF2 = is_homologue(SMILES.x, SMILES.y, "F[C+2]F"),
         pattern_CF2CF2 = is_homologue(SMILES.x, SMILES.y, "F[C+](F)[C+](F)F")) %>% #  "F[C+2]F"
  ungroup()

homologs_CF2 <- homologs %>%
  filter(!is.na(pattern_CF2))

homologs_CF2CF2 <- homologs %>%
  filter(!is.na(pattern_CF2CF2))


# save out actual concentrations

Orbitrap_dataset_raw = read_excel_allsheets(filename = "data_for_modelling/Batch 1 Semi Quant w frag.xlsx")

Orbitrap_dataset_raw = Orbitrap_dataset_raw %>% 
  group_by(Compound) %>%
  mutate(Theoretical_amt = case_when(
    Filename == "2020071205-cal21" ~ mean(Theoretical_amt[Filename=="2020071205-cal22"]),
    TRUE ~ Theoretical_amt))%>%
  ungroup()

data_real_conc = Orbitrap_dataset_raw %>%
  mutate(Theoretical_amt = replace(Theoretical_amt , grepl("NaN", Theoretical_amt, fixed = TRUE), NA)) %>%
  group_by(Compound) %>%
  mutate(Theoretical_amt = case_when(
    Filename == "2020071205-cal21" ~ mean(Theoretical_amt[Filename=="2020071205-cal22"]),
    TRUE ~ Theoretical_amt))%>%
  ungroup() %>%
  left_join(SMILES_data) %>%
  drop_na(SMILES) %>%
  mutate(RT = as.numeric(RT),
         area_IC = Area*IC,
         Theoretical_conc_uM = Theoretical_amt/Molecular_weight) 
#  %>%
# group_by(SMILES, Compound) %>%
# mutate(slope = linear_regression(area_IC, Theoretical_conc_uM)$slope,
#        intercept = linear_regression(area_IC, Theoretical_conc_uM)$intercept) %>%
# ungroup()



# find homolog series quantifications

data_homolog_conc_CF2 <- concentration_forAnalytes_homolog(filename_data = "data_for_modelling/Batch 1 Semi Quant w frag.xlsx",
                                                           filename_smiles = "data_for_modelling/Smiles_for_Target_PFAS_semicolon.csv",
                                                           homolog_pattern_SMILES = "F[C+2]F",
                                                           findHomolog_onlyForAnalytes = FALSE) 


data_homolog_conc_CF2_intercept <- concentration_forAnalytes_homolog_withIntercept(filename_data = "data_for_modelling/Batch 1 Semi Quant w frag.xlsx",
                                                                                   filename_smiles = "data_for_modelling/Smiles_for_Target_PFAS_semicolon.csv",
                                                                                   homolog_pattern_SMILES = "F[C+2]F",
                                                                                   findHomolog_onlyForAnalytes = FALSE) 



#CF2CF2

data_homolog_conc_CF2CF2 <- concentration_forAnalytes_homolog(filename_data = "data_for_modelling/Batch 1 Semi Quant w frag.xlsx",
                                                              filename_smiles = "data_for_modelling/Smiles_for_Target_PFAS_semicolon.csv",
                                                              homolog_pattern_SMILES = "F[C+](F)[C+](F)F",
                                                              findHomolog_onlyForAnalytes = FALSE) 


data_homolog_conc_CF2CF2_intercept <- concentration_forAnalytes_homolog_withIntercept(filename_data = "data_for_modelling/Batch 1 Semi Quant w frag.xlsx",
                                                                                   filename_smiles = "data_for_modelling/Smiles_for_Target_PFAS_semicolon.csv",
                                                                                   homolog_pattern_SMILES = "F[C+](F)[C+](F)F",
                                                                                   findHomolog_onlyForAnalytes = FALSE) 



# ----------------------------------------
# Use LOO models and predict concentrations (using other 32 PFAS as calibrants)
# ----------------------------------------

# for CF2 compounds
SMILES_list_homolog_CF2 <- homologs_CF2 %>%
  select(Compound, SMILES.x) %>%
  rename(SMILES = SMILES.x) %>% 
  unique()


SMILES_names_with_homologues_CF2 = SMILES_list_homolog_CF2
table_with_PFAS_LOO_pred_logIEs = read_delim("results/modelling_results/PFAS_pred_logIEs_with_leave_one_out_approach.csv")
directory_with_LOO_models = "C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant/models/leave_one_out_approach"
data_detected_PFAS = data_real_conc

model_pred_CF2_homologues = concentration_forPFAS_pretrained_models(SMILES_names_with_homologues = SMILES_names_with_homologues_CF2,
                                                                    table_with_PFAS_LOO_pred_logIEs = table_with_PFAS_LOO_pred_logIEs,
                                                                    directory_with_LOO_models = directory_with_LOO_models,
                                                                    data_detected_PFAS = data_detected_PFAS)


# for CF2CF2 compounds
SMILES_list_homolog_CF2CF2 <- homologs_CF2CF2 %>%
  select(Compound, SMILES.x) %>%
  rename(SMILES = SMILES.x) %>% 
  unique()

SMILES_names_with_homologues_CF2CF2 = SMILES_list_homolog_CF2CF2

model_pred_CF2CF2_homologues = concentration_forPFAS_pretrained_models(SMILES_names_with_homologues = SMILES_names_with_homologues_CF2CF2,
                                                                       table_with_PFAS_LOO_pred_logIEs = table_with_PFAS_LOO_pred_logIEs,
                                                                       directory_with_LOO_models = directory_with_LOO_models,
                                                                       data_detected_PFAS = data_detected_PFAS)



# Comparison in table

#CF2
# summary_table_CF2  <- data_real_conc %>%
#   left_join(predicted_concentrations %>%
#               select(slope_pred, conc_pred, Compound, Filename, Area, SMILES), 
#             by = c("Compound", "Filename", "Area", "SMILES")) %>%
#   left_join(data_homolog_conc_CF2 %>%
#               rename(Compound = Compound.y,
#                      SMILES = SMILES.y,
#                      Compound_homolog = Compound.x,
#                      SMILES_homolog = SMILES.x,
#                      slope_homolog = slope,
#                      conc_homolog = conc))

#CF2
summary_table_CF2  <- model_pred_CF2_homologues$predicted_conc %>%
  left_join(data_homolog_conc_CF2 %>%
              rename(Compound = Compound.y,
                     SMILES = SMILES.y,
                     Compound_homolog = Compound.x,
                     SMILES_homolog = SMILES.x,
                     slope_homolog = slope,
                     conc_homolog_uM = conc)) %>%
  left_join(data_homolog_conc_CF2_intercept %>%
              rename(Compound = Compound.y,
                     SMILES = SMILES.y,
                     Compound_homolog = Compound.x,
                     SMILES_homolog = SMILES.x,
                     slope_homolog = slope,
                     conc_homolog_withIntercept_uM = conc,
                     intercept_homolog = intercept))




# write_delim(summary_table_CF2, "results/homologue_vs_IEmodel_results/summary_table_CF2_concentrations.csv")

#CF2CF2
summary_table_CF2CF2  <- model_pred_CF2CF2_homologues$predicted_conc %>%
  left_join(data_homolog_conc_CF2CF2 %>%
              rename(Compound = Compound.y,
                     SMILES = SMILES.y,
                     Compound_homolog = Compound.x,
                     SMILES_homolog = SMILES.x,
                     slope_homolog = slope,
                     conc_homolog_uM = conc)) %>%
  left_join(data_homolog_conc_CF2CF2_intercept %>%
              rename(Compound = Compound.y,
                     SMILES = SMILES.y,
                     Compound_homolog = Compound.x,
                     SMILES_homolog = SMILES.x,
                     slope_homolog = slope,
                     conc_homolog_withIntercept_uM = conc,
                     intercept_homolog = intercept))



#write_delim(summary_table_CF2CF2, "results/homologue_vs_IEmodel_results/summary_table_CF2CF2_concentrations.csv")



# plots

IE_c_plot = ggplot(data = summary_table_CF2CF2)+
  geom_point(mapping = aes(x = Theoretical_conc_uM,
                           y = conc_pred_uM,
                           color = Compound)) +
  scale_y_log10(limits = c(10^-5, 10^0)) +
  scale_x_log10(limits = c(10^-5, 10^0)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_abline(slope = 1, intercept = 1) +
  geom_abline(slope = 1, intercept = -1) +
  theme(aspect.ratio = 1)

homolog_c_plot = ggplot(data = summary_table_CF2CF2)+
  geom_point(mapping = aes(x = Theoretical_conc_uM,
                           y = conc_homolog_uM,
                           color = Compound,
                           text = Compound_homolog)) +
  scale_y_log10(limits = c(10^-5, 10^0)) +
  scale_x_log10(limits = c(10^-5, 10^0)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_abline(slope = 1, intercept = 1) +
  geom_abline(slope = 1, intercept = -1) +
  geom_abline(slope = 1, intercept = 1) + 
  theme(aspect.ratio = 1)#,
#legend.position = "none")
ggplotly(homolog_c_plot)


plot_comp <- plot_grid(IE_c_plot, homolog_c_plot)

plot_comp

# Checking cal graphs and intercepts differences

cal_graph = ggplot(data = summary_table_CF2_filtered %>%
                     filter(Compound == "PFNA"))+
  geom_point(mapping = aes(x = Theoretical_conc_uM,
                           y = Area)) +
  geom_abline(slope = summary(lm(summary_table_CF2_filtered$Area ~ summary_table_CF2_filtered$Theoretical_conc_uM))$coefficients[2],
              intercept = summary(lm(summary_table_CF2_filtered$Area ~ summary_table_CF2_filtered$Theoretical_conc_uM))$coefficients[1])+
  geom_abline(slope = summary_table_CF2_filtered$slope_homolog, intercept = summary_table_CF2_filtered$intercept) +
  theme(aspect.ratio = 1)
cal_graph



data_forCal <- data_real_conc %>%
  filter(Compound == "PFNA")

data_homolog_quant <- data_real_conc %>%
  filter(Compound == "PFDA")

cal_graph = ggplot(data = data_forCal)+
  geom_point(mapping = aes(x = Theoretical_conc_uM,
                           y = area_IC)) +
  # geom_abline(slope = summary(lm(summary_table_CF2_filtered$Area ~ summary_table_CF2_filtered$Theoretical_conc_uM))$coefficients[2],
  #             intercept = summary(lm(summary_table_CF2_filtered$Area ~ summary_table_CF2_filtered$Theoretical_conc_uM))$coefficients[1])+
  geom_abline(slope = data_forCal$slope, intercept = data_forCal$intercept) +
  geom_abline(slope = data_homolog_quant$slope, intercept = data_homolog_quant$intercept, colour='#E41A1C') +
  theme(aspect.ratio = 1)
cal_graph
print(data_forCal$slope[1])
print(data_forCal$intercept[1])



# Error calculations
summary_table_CF2_filtered = summary_table_CF2_filtered %>%
  mutate(error_IE = case_when(
    Theoretical_conc_uM > conc_pred ~ Theoretical_conc_uM/conc_pred,
    TRUE ~ conc_pred/Theoretical_conc_uM),
    error_homolog = case_when(
      Theoretical_conc_uM > conc_homolog ~ Theoretical_conc_uM/conc_homolog,
      TRUE ~ conc_homolog/Theoretical_conc_uM),)

summary_table_CF2_filtered %>%
  na.omit() %>%
  group_by(pattern) %>%
  summarize(error_IE = mean(error_IE),
            error_homolog = mean(error_homolog)) %>%
  ungroup()






#---------------------
# Testing Lara's data
#---------------------

#alldata_lara <- read_excel_allsheets("C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/Topics_exp_codes/Lara_PFAS/20220208_Suspect_screening_TU pools.xlsx")

#SMILES_data <- read_SMILES(filename = "C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/Topics_exp_codes/Lara_PFAS/Smiles_for_Target_PFAS_semicolon_lara.csv")

lara_concentrations_pred <- concentration_forAnalytes_model(filename_data = "C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/Topics_exp_codes/Lara_PFAS/20220208_Suspect_screening_TU pools.xlsx",
                                                            filename_smiles = "C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/Topics_exp_codes/Lara_PFAS/Smiles_for_Target_PFAS_semicolon_lara.csv",
                                                            filename_eluent = "data/eluent.csv",
                                                            pred_model =  logIE_pred_model,
                                                            compounds_to_be_removed_as_list = c(
                                                              "HFPO-DA", "MeFOSE", "EtFOSE", "10:2 mono PAP", "4:2 mono PAP", "6:2 mono PAP", "8:2 mono PAP")) #by names, all of those exist in both Lara's and Thomas' data

setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant_HS")
cal_filename_data <-  "C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/projects_measurements/Lara_PFAS/Target PFAS for semiquant_120722.xlsx"
cal_filename_smiles <- "C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/projects_measurements/Lara_PFAS/Target_PFAS_Calibrants.csv"
sus_filename_data <- "C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/projects_measurements/Lara_PFAS/20220208_Suspect_screening_TU pools.xlsx"
sus_filename_smiles <- "C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/projects_measurements/Lara_PFAS/suspects_smiles_neutral.csv"
logIE_pred_model <- readRDS("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant_HS/model_PFAS_logIE.Rdata")

lara_concentrations_pred <- concentration_forAnalytes_model_cal_separateFile(cal_filename_data, 
                                                                             cal_filename_smiles, 
                                                                             sus_filename_data,
                                                                             sus_filename_smiles,
                                                                             filename_eluent = "data/eluent.csv",
                                                                             pred_model =  logIE_pred_model,
                                                                             compounds_to_be_removed_as_list = c("PFPeS", "PFHpS", "PFNS", "PFPeDA", "10:2 mono PAP", "4:2 mono PAP", "6:2 mono PAP", "8:2 mono PAP"))

#data_ions <- lara_concentrations_pred$data
#data_neutral <- lara_concentrations_pred$data
#saveRDS(lara_concentrations_pred, file="C:/Users/karpa/OneDrive - Kruvelab/Helen_phd/projects_measurements/Lara_PFAS/lara_PFAS_data2.RData")

lara_pred <- lara_concentrations_pred$data

lara_pred <- lara_pred %>%
  select(-c(IC, Molecular_weight, area_IC)) %>%
  rename(Predicted_RF = slope_pred,
         Predicted_conc_uM = conc_pred)
#write_delim(lara_pred, "C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/projects_measurements/Lara_PFAS/Lara_pred_conc_neutralSuspectsSMILES.csv", delim = ",")

ggplotly(lara_concentrations_pred$plot_predicted_theoretical_conc)
ggplotly(lara_concentrations_pred$plot_predictedIE_slope)


data_both <- data_neutral %>%
  rename(conc_pred_neutral = conc_pred) %>%
  #filter(Compound == "eecec PFSA n=8") %>% 
  left_join(data_ions %>%
              select(Compound, Filename, Area, RT, conc_pred) %>% 
              rename(conc_pred_ion = conc_pred)) %>% 
  mutate(conc_difference = abs(conc_pred_neutral-conc_pred_ion))
  filter(conc_pred_neutral != conc_pred_ion)


ggplot(data = data_both) +
  geom_point(mapping = aes(x = log10(conc_pred_neutral), 
                           y = log10(conc_pred_ion),
                           color = Compound))+
  xlim(c(-8, -4)) +
  ylim(c(-8,-4)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 1, slope = 1) +
  geom_abline(intercept = -1, slope = 1) +
  
  my_theme




# barplot of suspect concentrations?

#bar plot
data_short = lara_pred %>%
  filter(Compound %in% c("d/C PFSA n=8", "eecec PFSA n=8", "ether PFSA n=4", "ether PFSA n=8", "H-PFDoDA", "H-PFDS", "NMe-FBSAA")) %>%
  #filter(Compound %in% c("NMe-FBSAA")) %>% 
  mutate(unique_compound = paste0(Compound, " (", SMILES, ")")) %>%
  filter(grepl("Pool", Filename, fixed = TRUE))

barplot <- ggplot(data = data_short, aes( x = factor( unique_compound ), y = ((Predicted_conc_uM)), fill = Compound ) ) +    # print bar chart
  geom_bar( stat = 'identity', position = 'dodge') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #scale_fill_manual(values=c("#959D95", "#515251", "#7CB368")) +
  xlab("") +
  ylab("Concentration (nM)")+
  my_theme +
  facet_wrap(~ Filename, ncol = 6)
  #coord_flip()
#scale_log_x

#ggsave(barplot, filename = "C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/projects_measurements/Lara_PFAS/suspects_summary_byFilename.svg", width=40, height=80, units = "cm")



















