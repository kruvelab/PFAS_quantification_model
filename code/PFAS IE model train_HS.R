
#setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFAS_semi_quant_HS")
setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant_HS")
#setwd("/GitHub/PFAS_semi_quant_HS")
source("code/functions.R")
library(caTools)
library(tidyverse)
library(caret)
library(plotly)
library(cowplot)

#-------------------------------------------------------------
# Building model for IE pred in neg mode (joining Jaanus's data)
#-------------------------------------------------------------

## ---- Reading in LC-MS data of calibration solutions ----
Orbitrap_dataset_raw = read_excel_allsheets(filename = "data/Batch 1 Semi Quant w frag.xlsx")

Orbitrap_dataset_raw = Orbitrap_dataset_raw %>% 
  group_by(Compound) %>%
  mutate(Theoretical_amt = case_when(
    Filename == "2020071205-cal21" ~ mean(Theoretical_amt[Filename=="2020071205-cal22"]),
    TRUE ~ Theoretical_amt))%>%
  ungroup() %>%
  filter(Theoretical_amt != "NaN")

## ---- Reading in SMILES for calibration compounds, removing NAs and adducts, mono PAPs, HFPO-DA ----
SMILES_data = read_SMILES(filename = "data/Smiles_for_Target_PFAS_semicolon.csv",
                          compounds_to_be_removed_as_list = c("HFPO-DA", "MeFOSE", "EtFOSE", "10:2 mono PAP", "4:2 mono PAP", "6:2 mono PAP", "8:2 mono PAP"))

## ---- Joining all collected data to one tibble, removing missing values, calculating slopes ----
data = Orbitrap_dataset_raw %>%
  left_join(SMILES_data) %>%
  drop_na(SMILES) %>%
  mutate(RT = as.numeric(RT),
         area_IC = Area*IC,
         Theoretical_conc_uM = Theoretical_amt/Molecular_weight) %>%       
  group_by(SMILES, Compound) %>%
  summarize(slope = linear_regression(area_IC, Theoretical_conc_uM)$slope,
            RT = mean(RT)) %>%
  ungroup()

data = add_mobile_phase_composition(data = data,
                                    eluent_file_name = "data/eluent.csv")

## ---- Converting slopes to logIE with PFOS as anchor ----
training = anchoring(data_to_be_anchored = data,
                     data_containing_anchor = "data/IE_training_data/190714_negative_model_logIE_data.csv")

## ---- Calculating PaDEL descriptors to all compounds based on SMILES ----
data_all_binded = PaDEL_original(training)

## ---- Cleaning data ----
data_clean = cleaning_data(data_all_binded)

## ---- Training the model with train/test split of 0.8 ----
logIE_pred_model_train_test = training_logIE_pred_model(data = data_clean,
                                                        split = 0.8)

IE_slope_cor = ggplot() +
  geom_point(data = logIE_pred_model_train_test$data$training_set,
             mapping = aes(logIE, logIE_predicted),
             color = "light grey",
             alpha = 0.5,
             size = 3) +
  geom_point(data = logIE_pred_model_train_test$data$test_set,
             mapping = aes(logIE, logIE_predicted),
             color = "blue",
             alpha = 0.5,
             size = 3) +
  labs(title = "Training set", 
       x = "Measured logIE",
       y = "Predicted logIE")+
  theme(axis.text = element_text(size=12),
        plot.title = element_text(size = 12),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.y = element_line(size = 1, color = "black"),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))+
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~data_type)

IE_slope_cor


## ---- Training the model with all data ----
logIE_pred_model = training_logIE_pred_model(data = data_clean)












