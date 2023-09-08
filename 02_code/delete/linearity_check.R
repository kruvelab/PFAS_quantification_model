setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/PFOA_semi_quant")
source("code/functions.R")
library(caTools)
library(tidyverse)
library(caret)
library(plotly)
library(cowplot)

# linear reg with log values
linear_regression <- function(y, x) {
  # if multiple replicas with same concentrations - average
  df <- tibble(x = x,
               y = y) %>%
    group_by(x) %>%
    mutate(y = mean(y)) %>%
    ungroup() %>%
    unique() %>%
    na.omit()
  
  y = log10(df$y)
  x = log10(df$x)
  
  plot_slopes_list = list()
  
  if(length(x) == 0){
    slope = NA
    intercept = NA
    regression_parameters <- list("slope" = slope, "intercept" = intercept)
    return(regression_parameters)
  } else if (length(y) > 5) {
    for (i in length(y):5){
      lm_summary = summary(lm(y ~ x))
      slope = lm_summary$coefficients[2]
      intercept = lm_summary$coefficients[1]
      residuals = (y - (slope*x +intercept))/y*100

      plot_here = ggplot() +
        geom_point(mapping = aes(x = log10(df$x), 
                                 y = log10(df$y))) +
        geom_abline(slope = slope, intercept = intercept) +
        theme_bw() +
        labs(x = "log-conc_uM",
             y = "log-area")
      plot_slopes_list[[i]] = plot_here
      
      y_non_log = 10^y
      x_non_log = 10^x
      reg_non_log = summary(lm(y_non_log ~ x_non_log))
      regression_parameters <- list("slope" = reg_non_log$coefficients[2], 
                                    "intercept" = reg_non_log$coefficients[1])
      if (max(abs(residuals)) < 10) {
        #return(regression_parameters)
        break
      }
      y = y[1:(length(y)-1)]
      x = x[1:(length(x)-1)]
    }
    return(list("param" = regression_parameters,
                "plots" = plot_slopes_list,
                "resid" = residuals))
  } else {
    y_non_log = 10^y
    x_non_log = 10^x
    reg_non_log = summary(lm(y_non_log ~ x_non_log)) 
    slope = reg_non_log$coefficients[2]
    intercept = reg_non_log$coefficients[1]
    regression_parameters <- list("slope" = slope, "intercept" = intercept)
    return(regression_parameters)
  }
}

#linear regression with leave out from high conc
linear_regression <- function(y, x) {
  # if multiple replicas with same concentrations - average
  df <- tibble(x = x,
               y = y) %>%
    group_by(x) %>%
    mutate(y = mean(y)) %>%
    ungroup() %>%
    unique() %>%
    na.omit()
  
  y = df$y
  x = df$x
  
  plot_slopes_list = list()
  
  if(length(x) == 0){
    slope = NA
    intercept = NA
    regression_parameters <- list("slope" = slope, "intercept" = intercept)
    return(regression_parameters)
  } else if (length(y) > 5) {
    for (i in length(y):5){
      lm_summary = summary(lm(y ~ x))
      slope = lm_summary$coefficients[2]
      intercept = lm_summary$coefficients[1]
      residuals = (y - (slope*x +intercept))/y*100
      
      plot_here = ggplot() +
        geom_point(mapping = aes(x = df$x, 
                                 y = df$y)) +
        geom_abline(slope = slope, intercept = intercept) +
        theme_bw() +
        labs(x = "conc_uM",
             y = "area")
      plot_slopes_list[[i]] = plot_here

      # regression_parameters <- list("slope" = slope, 
      #                               "intercept" = intercept)
      if (max(abs(residuals)) < 10) {
        break
      }
      y = y[1:(length(y)-1)]
      x = x[1:(length(x)-1)]
    }
    return(list("slope" = slope,
                "intercept" = intercept,
                "plots" = plot_slopes_list,
                "resid" = residuals))
  } else {
    lm_summary = summary(lm(y ~ x))
    slope = lm_summary$coefficients[2]
    intercept = lm_summary$coefficients[1]
    regression_parameters <- list("slope" = slope, "intercept" = intercept)
    return(regression_parameters)
  }
}


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
         Theoretical_conc_uM = Theoretical_amt/Molecular_weight) # %>%


# in loop, generate cal plots to all
plot_list = list()

for (i in 1:length(unique(data$Compound))) {
  
  cal_data = data %>% 
    filter(Compound == unique(data$Compound)[i]) %>% 
    mutate(cal_slope = linear_regression(area_IC, Theoretical_conc_uM)$slope)

  plot = ggplot(data = cal_data) +
    geom_point(mapping = aes(Theoretical_conc_uM, area_IC),
               color = "#274c77",
               alpha = 0.7,
               size = 2) +
    geom_abline(slope = cal_data$cal_slope, intercept = 0) +
    theme_bw()
  
  plot_list[[i]] = plot
}


y = cal_data$area_IC
x = cal_data$Theoretical_conc_uM

results = linear_regression(y,x)

results$resid
results$plots

save_plot <-  plot_grid(results$plots[[8]], results$plots[[7]], results$plots[[6]], results$plots[[5]], ncol = 1)

ggsave(save_plot, filename="C:/Users/HelenSepman/OneDrive - Kruvelab/Helen_phd/projects_measurements/Melanie/linearity/plot_remove_last_log10.png", height = 7, width = 3.21)






