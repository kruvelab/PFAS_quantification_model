


#-----------------------------------------------------------
# Homologous series vs ML model predictions - concentrations
#-----------------------------------------------------------

# find homolog series compounds
# save out actual concentrations
# find homolog series quantifications
# train model by leaving one out each time
# Comparison in table
# profit


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

Orbitrap_dataset_raw = read_excel_allsheets(filename = "data/Batch 1 Semi Quant w frag.xlsx")

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

data_homolog_conc_CF2 <- concentration_forAnalytes_homolog(filename_data = "data/Batch 1 Semi Quant w frag.xlsx",
                                                           filename_smiles = "data/Smiles_for_Target_PFAS_semicolon.csv",
                                                           homolog_pattern_SMILES = "F[C+2]F",
                                                           findHomolog_onlyForAnalytes = FALSE) 


data_homolog_conc_CF2_intercept <- concentration_forAnalytes_homolog_withIntercept(filename_data = "data/Batch 1 Semi Quant w frag.xlsx",
                                                                                   filename_smiles = "data/Smiles_for_Target_PFAS_semicolon.csv",
                                                                                   homolog_pattern_SMILES = "F[C+2]F",
                                                                                   findHomolog_onlyForAnalytes = FALSE) 



#CF2CF2

data_homolog_conc_CF2CF2 <- concentration_forAnalytes_homolog(filename_data = "data/Batch 1 Semi Quant w frag.xlsx",
                                                              filename_smiles = "data/Smiles_for_Target_PFAS_semicolon.csv",
                                                              homolog_pattern_SMILES = "F[C+](F)[C+](F)F",
                                                              findHomolog_onlyForAnalytes = FALSE) 





# train model by leaving one out each time

# for CF2 compounds
SMILES_list_homolog_CF2 <- homologs_CF2 %>%
  select(Compound, SMILES.x) %>%
  unique()

predicted_concentrations <- tibble()

for (i in 1:length(SMILES_list_homolog_CF2$SMILES.x)) {
  data_forTraining <- data_clean %>%
    filter(!grepl(SMILES_list_homolog_CF2[i,2], SMILES, fixed = TRUE))
  
  setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant_HS/code/homologue_series_models")
  logIE_pred_model_new = training_logIE_pred_model(data = data_forTraining,
                                                   save_model_name =  paste("without_homologue", i, sep = "_"))
  
  setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant_HS")
  conc_forAll_compounds <- concentration_forAnalytes_model(filename_data = "data/Batch 1 Semi Quant w frag.xlsx", 
                                                           filename_smiles = "data/Smiles_for_Target_PFAS_semicolon.csv", 
                                                           filename_eluent = "data/eluent.csv", 
                                                           logIE_pred_model_new,
                                                           compounds_to_be_removed_as_list = c("HFPO-DA", "MeFOSE", "EtFOSE", "10:2 mono PAP", "4:2 mono PAP", "6:2 mono PAP", "8:2 mono PAP"))
  
  
  conc_pred_forSuspect <- conc_forAll_compounds$data %>%
    filter(grepl(SMILES_list_homolog_CF2[i,2], SMILES, fixed = TRUE))
  
  predicted_concentrations <- predicted_concentrations %>%
    bind_rows(conc_pred_forSuspect)
}

# for CF2CF2 compounds
SMILES_list_homolog_CF2CF2 <- homologs_CF2CF2 %>%
  select(Compound, SMILES.x) %>%
  filter(!Compound %in% c("PFTriDA", "PFNA")) %>%
  unique()

predicted_concentrations_CF2CF2 <- tibble()

for (i in 4:length(SMILES_list_homolog_CF2CF2$SMILES.x)) {
  data_forTraining <- data_clean %>%
    filter(!grepl(SMILES_list_homolog_CF2CF2[i,2], SMILES, fixed = TRUE))
  
  logIE_pred_model_new = training_logIE_pred_model(data = data_forTraining)
  
  conc_forAll_compounds <- concentration_forAnalytes_model(filename_data = "data/Batch 1 Semi Quant w frag.xlsx", 
                                                           filename_smiles = "data/Smiles_for_Target_PFAS_semicolon.csv", 
                                                           filename_eluent = "data/eluent.csv", 
                                                           logIE_pred_model_new,
                                                           compounds_to_be_removed_as_list = c("HFPO-DA", "MeFOSE", "EtFOSE", "10:2 mono PAP", "4:2 mono PAP", "6:2 mono PAP", "8:2 mono PAP"))
  
  
  conc_pred_forSuspect <- conc_forAll_compounds$data %>%
    filter(grepl(SMILES_list_homolog_CF2CF2[i,2], SMILES, fixed = TRUE))
  
  predicted_concentrations_CF2CF2 <- predicted_concentrations_CF2CF2 %>%
    bind_rows(conc_pred_forSuspect)
}





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
summary_table_CF2  <- data_real_conc %>%
  left_join(predicted_concentrations %>%
              select(slope_pred, conc_pred, Compound, Filename, Area, SMILES), 
            by = c("Compound", "Filename", "Area", "SMILES")) %>%
  left_join(data_homolog_conc_CF2 %>%
              rename(Compound = Compound.y,
                     SMILES = SMILES.y,
                     Compound_homolog = Compound.x,
                     SMILES_homolog = SMILES.x,
                     slope_homolog = slope,
                     conc_homolog = conc)) %>%
  left_join(data_homolog_conc_CF2_intercept %>%
              rename(Compound = Compound.y,
                     SMILES = SMILES.y,
                     Compound_homolog = Compound.x,
                     SMILES_homolog = SMILES.x,
                     slope_homolog = slope,
                     conc_homolog_withIntercept = conc,
                     intercept_homolog = intercept))



summary_table_CF2_filtered <- summary_table_CF2 %>%
  filter(!is.na(Compound_homolog))

write.csv2(summary_table_CF2_filtered, "C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/Topics_exp_codes/Lara_PFAS/Melanie/summary_table_CF2_withintercept_filtered.csv")

#CF2CF2
summary_table_CF2CF2  <- data_real_conc %>%
  left_join(predicted_concentrations_CF2CF2 %>%
              select(slope_pred, conc_pred, Compound, Filename, Area, SMILES), 
            by = c("Compound", "Filename", "Area", "SMILES")) %>%
  left_join(data_homolog_conc_CF2CF2 %>%
              rename(Compound = Compound.y,
                     SMILES = SMILES.y,
                     Compound_homolog = Compound.x,
                     SMILES_homolog = SMILES.x,
                     slope_homolog = slope,
                     conc_homolog = conc))


summary_table_CF2CF2_filtered <- summary_table_CF2CF2 %>%
  filter(!is.na(Compound_homolog))

write.csv2(summary_table_CF2CF2_filtered, "C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/Topics_exp_codes/Lara_PFAS/Melanie/summary_table_CF2CF2_filtered.csv")

# plots

IE_c_plot = ggplot(data = summary_table_CF2_filtered)+
  geom_point(mapping = aes(x = Theoretical_conc_uM,
                           y = conc_pred,
                           color = Compound)) +
  scale_y_log10(limits = c(10^-5, 10^0)) +
  scale_x_log10(limits = c(10^-5, 10^0)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_abline(slope = 1, intercept = 1) +
  geom_abline(slope = 1, intercept = -1) +
  theme(aspect.ratio = 1)

homolog_c_plot = ggplot(data = summary_table_CF2_filtered)+
  geom_point(mapping = aes(x = Theoretical_conc_uM,
                           y = conc_homolog,
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
                                                            compounds_to_be_removed_as_list = c("HFPO-DA", "MeFOSE", "EtFOSE", "10:2 mono PAP", "4:2 mono PAP", "6:2 mono PAP", "8:2 mono PAP")) #by names, all of those exist in both Lara's and Thomas' data

lara_pred <- lara_concentrations_pred$data

ggplotly(lara_concentrations_pred$plot_predicted_theoretical_conc)









