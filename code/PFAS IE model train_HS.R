setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant")
source("code/functions.R")
library(caTools)
library(tidyverse)
library(caret)
library(plotly)
library(cowplot)

#-------------------------------------------------------------
# Building model for IE pred in neg mode (joining Liigand's data)
#-------------------------------------------------------------

## ---- Reading in LC-MS data of calibration solutions ----
Orbitrap_dataset_raw = read_excel_allsheets(filename = "data_for_modelling/Batch 1 Semi Quant w frag.xlsx")

Orbitrap_dataset_raw = Orbitrap_dataset_raw %>%
  group_by(Compound) %>%
  mutate(Theoretical_amt = case_when(
    Filename == "2020071205-cal21" ~ mean(Theoretical_amt[Filename=="2020071205-cal22"]),
    TRUE ~ Theoretical_amt))%>%
  ungroup() %>%
  filter(Theoretical_amt != "NaN")

## ---- Reading in SMILES for calibration compounds, removing NAs and adducts, mono PAPs, HFPO-DA ----
SMILES_data = read_SMILES(filename = "data_for_modelling/Smiles_for_Target_PFAS_semicolon.csv",
                          compounds_to_be_removed_as_list = c("HFPO-DA", "MeFOSE", "EtFOSE", "10:2 mono PAP", "4:2 mono PAP", "6:2 mono PAP", "8:2 mono PAP"))

## ---- Joining all collected data to one tibble, removing missing values, calculating slopes ----
data = Orbitrap_dataset_raw %>%
  left_join(SMILES_data) %>%
  drop_na(SMILES) %>%
  mutate(RT = as.numeric(RT),
         area_IC = Area*IC,
         Theoretical_conc_uM = Theoretical_amt/Molecular_weight) %>%
  group_by(SMILES, Compound) %>%
  summarize(slope = linear_regression(area_IC, Theoretical_conc_uM, remove_lowest_conc = TRUE)$slope,
            RT = mean(RT)) %>%
  ungroup()

data = add_mobile_phase_composition(data = data,
                                    eluent_file_name = "data_for_modelling/eluent.csv")

## ---- Converting slopes to logIE with PFOS as anchor ----
training = anchoring(data_to_be_anchored = data,
                     data_containing_anchor = "data_for_modelling/IE_training_data/190714_negative_model_logIE_data.csv")

## If only modelling with original data from Liigand et al
# training = training %>%
#   filter(data_type == "non-PFAS")

## ---- Calculating PaDEL descriptors to all compounds based on SMILES ----
data_all_binded = PaDEL_original(training)

## ---- Cleaning data ----
data_clean = cleaning_data(data_all_binded)

## ---- Training the model with train/test split of 0.8 ----
logIE_pred_model_train_test = training_logIE_pred_model(data = data_clean,
                                                        split = 0.8)

# metrics of trained model
logIE_pred_model_train_test$metrics

#---- correlation plot ----
IE_slope_cor = ggplot() +
  geom_point(data = logIE_pred_model_train_test$data$training_set,
             mapping = aes(logIE, logIE_predicted)) +
  geom_point(data = logIE_pred_model_train_test$data$test_set,
             mapping = aes(logIE, logIE_predicted),color = "red") +
  geom_abline(intercept = -1, slope = 1) +
  geom_abline(intercept =  1, slope = 1) +
  geom_abline(intercept =  0, slope = 1) +
  xlab(substitute(paste("log", italic("IE"))["predicted"]))  +
  ylab(substitute(paste("log", italic("IE"))["measured"])) +
  theme_classic() +
  theme(aspect.ratio = 1,
        legend.key = element_blank()) +
  facet_wrap(~data_type)
 
IE_slope_cor


# Save model, data, metrics and varImp
#saveRDS(logIE_pred_model_train_test, file="models/230619_logIE_model_withPFAS_train_test.RData")


#-------------------------------------------
# ---- Training the model with all data ----
#-------------------------------------------
logIE_pred_model = training_logIE_pred_model(data = data_clean,
                                             bestTune_optimized_model = logIE_pred_model_train_test$model$bestTune,
                                             split = NULL)

logIE_pred_model$metrics
#saveRDS(logIE_pred_model, file="models/230619_logIE_model_withPFAS_allData.RData")


# ----------------------------------------------
# Modelling with leave-one-out approach for PFAS
# ----------------------------------------------

# train model by leaving one out each time

SMILES_list_PFAS <- data %>%
  select(SMILES, Compound) %>%
  unique()

predicted_PFAS_IEs <- tibble()
dim_data = tibble()

for (i in 1:length(SMILES_list_PFAS$SMILES)) {
  #remove one PFAS from the training set
  #NB, PFOS needs to be reomoved from old and new datasets!
  if(SMILES_list_PFAS[i,2][[1]] == "PFOS") {
    data_forTraining <- data_clean %>%
      filter(data_clean$SMILES != SMILES_list_PFAS[i,1][[1]]) %>%  #!grepl(SMILES_list_PFAS[i,1], SMILES, fixed = TRUE)) %>% 
      filter(!grepl("perfluorooctanesulfonic acid", name, fixed = TRUE))
  } else {
    data_forTraining <- data_clean %>%
    filter(SMILES != SMILES_list_PFAS[i,1][[1]])
  }
  print(SMILES_list_PFAS[i,2][[1]])
  dim(data_forTraining)[1]
  test_object = tibble(compound = c(SMILES_list_PFAS[i,2][[1]]),
                       dimensions = c(dim(data_forTraining)[1]))
  dim_data=dim_data %>% 
    bind_rows(test_object)
  
  #model
  logIE_pred_model_new = training_logIE_pred_model(data = data_forTraining,
                                                   split = 1,
                                                   save_model_name =  paste("models/leave_one_out_approach/model_", i, sep = ""))
  
  #predict
  IE_pred_for_PFAS = data_clean %>%
    filter(SMILES == SMILES_list_PFAS[i,1][[1]])#grepl(SMILES_list_PFAS[i,1], SMILES, fixed = TRUE)) 
  IE_pred_for_PFAS = IE_pred_for_PFAS %>% 
    mutate(logIE_predicted = predict(logIE_pred_model_new$model, newdata = IE_pred_for_PFAS)) %>% 
    mutate(model_nr = paste("model_", i, sep = "")) %>% 
    select(logIE_predicted, everything())
  
  predicted_PFAS_IEs <- predicted_PFAS_IEs %>%
    bind_rows(IE_pred_for_PFAS)
}

# #need to add slope as it is needed later
# predicted_PFAS_IEs = predicted_PFAS_IEs %>% 
#   left_join(data %>% select(SMILES, slope, area_IC, Theoretical_conc_uM))

#write_delim(predicted_PFAS_IEs, "results/modelling_results/PFAS_pred_logIEs_with_leave_one_out_approach.csv", delim =",")



# ------------------------------------------------------------
# ---- Calculating additional model evaluation parameters ----
# ------------------------------------------------------------

#logIE_pred_model_train_test <- readRDS(file="models/230329_logIE_model_withPFAS_train_test.RData")

#only PFAS
data_for_error_calc_train = logIE_pred_model_train_test$data$training_set %>%
  filter(data_type == "PFAS")

data_for_error_calc_test = logIE_pred_model_train_test$data$test_set %>%
  filter(data_type == "PFAS")

rmse_train <- rmse(data_for_error_calc_train$logIE, data_for_error_calc_train$logIE_predicted)
rmse_test <- rmse(data_for_error_calc_test$logIE, data_for_error_calc_test$logIE_predicted)

  
# mean error
logIE_pred_model_train_test_error <- logIE_pred_model_train_test$data$test_set %>% #PFAS_LOO_data %>%  
  #filter(data_type == "PFAS") %>% 
  mutate(pred_error = case_when(10^logIE > 10^logIE_predicted ~ 10^logIE/10^logIE_predicted,
                                TRUE ~ 10^logIE_predicted/10^logIE)) %>%
  # group_by(name) %>%
  # mutate(mean_pred_error = mean(pred_error)) %>%
  # ungroup() %>%
  select(pred_error, everything())


# mean pred error
mean(logIE_pred_model_train_test_error$pred_error)

# geometric mean
exp(mean(log(logIE_pred_model_train_test_error$pred_error)))

# median pred error
median(logIE_pred_model_train_test_error$pred_error)



# mean error of PFAS
logIE_pred_model_test_error_pfas <- logIE_pred_model_train_test$data$test_set %>%
  filter(data_type == "PFAS") %>% 
  mutate(pred_error = case_when(10^logIE > 10^logIE_predicted ~ 10^logIE/10^logIE_predicted,
                                TRUE ~ 10^logIE_predicted/10^logIE))

mean(logIE_pred_model_test_error_pfas$pred_error)












