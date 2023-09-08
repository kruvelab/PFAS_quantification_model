library(tidyverse)
library(plotly)
library(forcats)
source("code/compound_eluent.R")
setwd("~/GitHub/PFAS_semi_quant")


target_data = read_delim("data/Calculated_amounts_target.csv",
                         delim = ",",
                         col_names = TRUE)

SMILES_data = read_delim("data/Smiles_for_Target_PFAS_semicolon.csv",
                         delim = ";",
                         col_names = TRUE)

SMILES_data = SMILES_data %>%
  rename(Compound = ID) %>%
  select(Compound, SMILES, Class) %>%
  na.omit()

suspects_semi_quant = read_delim("data/suspect semi-quant amounts.csv",
                         delim = ",",
                         col_names = TRUE)

CIC_data = read_delim("data/Mass_balance_cetaceans.csv",
                      delim = ";",
                      col_names = TRUE)



#convert to predicted ng F/g

x <- SMILES_data$SMILES
f.atoms <- lengths(regmatches(x, gregexpr("F",x)))

fatoms. <- data.frame(f.atoms)

SMILES_data <- SMILES_data%>%
  cbind(fatoms.)%>%
  group_by(SMILES) %>%
  mutate(MW = molecularmass(SMILES))%>%
  ungroup()%>%
  mutate(mass_F = f.atoms*19,
         percent.F=mass_F/MW)



#cleaning data + calculating RSD
target_data = target_data %>%
  na.omit()


