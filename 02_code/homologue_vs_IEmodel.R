
#-----------------------------------------------------------
# Homologous series vs ML model predictions - concentrations
#-----------------------------------------------------------

# find homologue series compounds
# save out actual concentrations
# find homolog series quantifications
# train model by leaving one out each time
# Comparison in table

admin = "C:/Users/..."
setwd(admin)

#-------------------------------
# Find homolog series compounds
#-------------------------------
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
# write_delim(homologs_CF2 %>%
#               left_join(SMILES_forHomolog %>%
#                           rename(compound_homolog = Compound,
#                                  SMILES.y  = SMILES)), "results/homologue_vs_IEmodel_results/homologues_CF2.csv", delim = ",")

homologs_CF2CF2 <- homologs %>%
  filter(!is.na(pattern_CF2CF2))
# write_delim(homologs_CF2CF2%>%
#               left_join(SMILES_forHomolog %>%
#                           rename(compound_homolog = Compound,
#                                  SMILES.y  = SMILES)), "results/homologue_vs_IEmodel_results/homologues_CF2CF2.csv", delim = ",")


#----------------------------------
# save out real concentrations
#----------------------------------

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



#-----------------------------------------------------
# find concentrations using homologue response factor 
#-----------------------------------------------------

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
directory_with_LOO_models = "Documents/GitHub/PFOA_semi_quant/models/leave_one_out_approach"
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




# write_delim(summary_table_CF2, "results/homologue_vs_IEmodel_results/230703_summary_table_CF2_concentrations.csv")

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



#write_delim(summary_table_CF2CF2, "results/homologue_vs_IEmodel_results/230703_summary_table_CF2CF2_concentrations.csv")


summary_table_CF2 <- read_delim("results/homologue_vs_IEmodel_results/230703_summary_table_CF2_concentrations.csv")


# Error calculations
summary_table_CF2_filtered = summary_table_CF2 %>%
  mutate(error_IE = case_when(
    Theoretical_conc_uM > conc_pred_uM ~ Theoretical_conc_uM/conc_pred_uM,
    TRUE ~ conc_pred_uM/Theoretical_conc_uM),
    error_homolog = case_when(
      Theoretical_conc_uM > conc_homolog_uM ~ Theoretical_conc_uM/conc_homolog_uM,
      TRUE ~ conc_homolog_uM/Theoretical_conc_uM),)


summary_table_CF2_filtered %>%
  na.omit() %>%
  #group_by(pattern) %>%
  summarize(error_IE_mean = mean(error_IE),
            #error_IE_geommean = exp(mean(log(error_IE))),
            #error_IE_median = median(error_IE),
            error_homolog_mean = mean(error_homolog),
            #error_homolog_geommean = exp(mean(log(error_homolog))),
            #error_homolog_median = median(error_homolog)
            ) %>%
  ungroup()



