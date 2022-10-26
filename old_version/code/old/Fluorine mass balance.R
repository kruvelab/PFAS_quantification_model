library(tidyverse)
library(plotly)
library(forcats)
setwd("~/GitHub/PFAS_semi_quant")
source("code/PaDEL_descs_calculator.R")
source("code/reading_excel.R")
source("code/compound_eluent.R")

###fluorine mass balance###

target_data = read_delim("data/Calculated_amounts_target.csv",
                         delim = ",",
                         col_names = TRUE)
#fetching absolute responses from tracefinder data
target_data = target_data %>%
  na.omit()%>%
  
  
  
  CIC_data = read_delim("data/Smiles_for_Target_PFAS_semicolon.csv",
                        delim = ";",
                        col_names = TRUE)


#convert to predicted ng F/g
blended <- blended %>%
  mutate(mass_F = f.atoms*19,
         percent.F=mass_F/MW.x,
         "Predicted ng_F/uL"= (percent.F*pred_conc_pg_uL)/1000,
         pred_ng_F_g = (`Predicted ng_F/uL`*`Ex vol.`)/`sample weight`)

blended <- blended %>%
  group_by(`Sample ID`)%>%
  mutate(sum_suspects_ng_F_g = sum(pred_ng_F_g))%>%
  ungroup()