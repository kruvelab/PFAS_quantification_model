library(tidyverse)
library(plotly)
setwd("~/GitHub/PFOA_semi_quant")
source("code/PaDEL_descs_calculator.R")
source("code/reading_excel.R")
source("code/compound_eluent.R")

#setwd("C:/Users/annel/Nextcloud/mudeli script ja failid/PFOA_semi_quant/PFOA_semi_quant")

#regressor----

regressor = readRDS("regressors/regressor_neg_new.rds")
descs_names = readRDS("regressors/negative_descs.rds")

#lcms data ----

filename = "data/Batch 1 Semi Quant w frag.xlsx"
Orbitrap_dataset_raw = read_excel_allsheets(filename)
#Orbitrap_dataset_raw %>% select(Compound) %>% unique()
Orbitrap_dataset_raw = Orbitrap_dataset_raw %>%
  na.omit(Area) %>%
  filter(Area != "N/F") %>%
  mutate(Area = as.numeric(Area))

#Orbitrap_dataset_raw %>% select(Compound) %>% unique()

Orbitrap_dataset_raw = Orbitrap_dataset_raw %>%
  group_by(Compound) %>%
  mutate(`Theoretical Amt`=case_when(
    Filename=="2020071205-cal21"~mean(`Theoretical Amt`[Filename=="2020071205-cal22"]),
    TRUE~`Theoretical Amt`
  ))%>%ungroup()

#Orbitrap_dataset_raw %>% select(Compound) %>% unique()

#smiles and descriptors----

SMILES_data = read_delim("data/Smiles_for_Target_PFAS_semicolon.csv",
                         delim = ";",
                         col_names = TRUE)

#how many unique SMILES are there
#SMILES_data %>% select(SMILES) %>% unique()

SMILES_data = SMILES_data %>%
  rename(Compound = ID) %>%
  select(Compound, SMILES, Class) %>%
  na.omit()

#descs_calc_PFOA = PaDEL_original(SMILES_data)

descs_calc_PFOA = read_delim("data/descs_calc.csv",
                              delim = ",",
                              col_names = TRUE)

descs_calc_PFOA = descs_calc_PFOA %>%
  select(Compound, SMILES, 
         #all_of(descs_names),
        everything())

#check number of unique analytes
#descs_calc_PFOA %>% select(SMILES) %>% unique()

descs_calc_PFOA = descs_calc_PFOA %>%
  group_by(SMILES) %>%
  mutate(IC = isotopedistribution(SMILES),
         MW = molecularmass(SMILES)) %>%
  ungroup()

#check number of unique analytes
#descs_calc_PFOA %>% select(SMILES) %>% unique()

#eluent---
eluent = read_delim("data/eluent.csv",
                    delim = ",",
                    col_names = TRUE)

organic_modifier = "MeCN"
pH.aq. = 7.0

#descs_calc_PFOA %>% select(SMILES) %>% unique()

data = Orbitrap_dataset_raw %>%
  left_join(descs_calc_PFOA)

#check number of unique analytes    
#Orbitrap_dataset_raw %>% select(Compound) %>% unique()
#data %>% select(SMILES) %>% unique()
  

data = data %>%
  mutate(
    RT = as.numeric(RT),
    area_IC = Area*IC,
    organic_modifier = organic_modifier,
    pH.aq. = pH.aq.,
    NH4 = 1, #1 if th ebuffer contains NH¤ ions , 0 if not. 
    organic = organicpercentage(eluent,RT),
    viscosity = viscosity(organic,organic_modifier),
    surface_tension = surfacetension(organic,organic_modifier),
    polarity_index = polarityindex(organic,organic_modifier)) 

#data %>% select(SMILES) %>% unique()
training = data %>%
  filter(!is.na(SMILES)) %>%
  filter(`Theoretical Amt` != "N/F") %>%
  filter(`Theoretical Amt` != "N/A") %>%
  mutate(`Theoretical Amt` = as.numeric(`Theoretical Amt`)) %>%
  mutate(`Theoretical Amt` = `Theoretical Amt`/MW) #correct with MW
#print(training %>% select(Compound) %>% unique(), n=40)

ggplot(data = training) +
  geom_point(mapping = aes(x = `Theoretical Amt`,
                           y = area_IC)) +
  facet_wrap(~Compound, scales = "free") +
  scale_x_log10() +
  scale_y_log10()

training = training %>%
  group_by(SMILES) %>%
  mutate(slope = linear_regression(area_IC, `Theoretical Amt`)$slope) %>%
  ungroup()
#print(training %>% select(Compound) %>% unique(), n=40)

IE_pred = training %>% 
  mutate(logIE_pred = 0) %>%
  na.omit()

#print(IE_pred %>% select(Compound) %>% unique(), n=40)

prediction =  predict(regressor, newdata = IE_pred, predict.all = TRUE)
prediction = prediction$aggregate
IE_pred <- IE_pred %>%
  mutate(logIE_pred = prediction) %>%
  select(SMILES,logIE_pred, everything())

#print(IE_pred %>% select(Compound) %>% unique(), n=40)

IE_pred = IE_pred %>%
  left_join(SMILES_data) %>%
  select(Compound, SMILES, Class, logIE_pred, slope, everything())

#print(IE_pred %>% select(Compound) %>% unique(), n=40)

IE_slope_cor = ggplot(data = IE_pred) +
  geom_point(mapping = aes(x = logIE_pred,
                           y = slope, 
                           text = Compound, 
                           color = Class)) +
  scale_y_log10() +
  facet_wrap(~Class)

IE_slope_cor

graph_1sttryPFAScal=ggplotly(IE_slope_cor)
graph_1sttryPFAScal

htmlwidgets::saveWidget(plotly::as_widget(graph_1sttryPFAScal), "1stryPFAScal.html")

slope_RT_cor = ggplot(data = IE_pred) +
  geom_point(mapping = aes(x = as.numeric(RT),
                           y = slope, 
                           text = Compound)) +
  scale_y_log10()

ggplotly(slope_RT_cor)

graph_slope_logP = ggplot(data = IE_pred) +
  geom_point(mapping = aes(x = ALogP,
                           y = slope, 
                           text = Compound)) +
  facet_wrap(~Class)

ggplotly(graph_slope_logP)
