
setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant")
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

logIE_pred_model_train_test$metrics

IE_slope_cor = ggplot() +
  geom_point(data = logIE_pred_model_train_test$data$training_set,
             mapping = aes(logIE, logIE_predicted),
             color = "#515251",
             alpha = 0.8,
             size = 3) +
  geom_point(data = logIE_pred_model_train_test$data$test_set,
             mapping = aes(logIE, logIE_predicted),
             color = "#7CB368",
             alpha = 0.8,
             size = 3) +
  geom_abline(intercept = -1, slope = 1) +
  geom_abline(intercept = 1, slope = 1) +
  labs(#title = "Training set", 
       x = "Measured logIE",
       y = "Predicted logIE")+
  theme(plot.title = element_text(size = 14),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.y = element_line(size = 1, color = "black"),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        aspect.ratio = 1,
        axis.text = element_text(family = font,
                                size = fontsize,
                               color = basecolor),
        legend.key = element_blank(),
        strip.background = element_blank(),
        text = element_text(family = font,
                            size = fontsize,
                            color = basecolor))+
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~data_type)

IE_slope_cor
ggsave(IE_slope_cor, "logIE_test_train.png", device = NULL)

ggplotly(IE_slope_cor)


## ---- Training the model with all data ----
logIE_pred_model = training_logIE_pred_model(data = data_clean)

logIE_pred_model$metrics

#saveRDS(logIE_pred_model, file="model_PFAS_logIE.RData")


# mean error 

logIE_pred_model_train_test_error <- logIE_pred_model_train_test$data$test_set %>%
  mutate(pred_error = case_when(10^logIE > 10^logIE_predicted ~ 10^logIE/10^logIE_predicted,
                                TRUE ~ 10^logIE_predicted/10^logIE)) %>%
  group_by(name) %>% 
  mutate(mean_pred_error = mean(pred_error)) %>% 
  ungroup() %>% 
  select(pred_error, mean_pred_error, everything())

# mean pred error
mean(logIE_pred_model_train_test_error$pred_error)

mean(unique(logIE_pred_model_train_test_error$mean_pred_error))







