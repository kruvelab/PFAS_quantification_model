library(tidyverse)
library(plotly)
setwd("~/GitHub/PFOA_semi_quant")
source("code/PaDEL_descs_calculator.R")
source("code/reading_excel.R")
source("code/compound_eluent.R")

#setwd("C:/Users/annel/Nextcloud/mudeli script ja failid/PFOA_semi_quant/PFOA_semi_quant")

#regressor----

regressor = readRDS("regressors/PFAS_FOREST.rds")
descs_names = readRDS("regressors/descs_PFASadd.rds")
#lcms data ----

filename = "data/Batch 1 Semi Quant w frag.xlsx"
Orbitrap_dataset_raw = read_excel_allsheets(filename)
Orbitrap_dataset_raw = Orbitrap_dataset_raw %>%
  na.omit(Area) %>%
  filter(Area != "N/F") %>%
  mutate(Area = as.numeric(Area))

Orbitrap_dataset_raw = Orbitrap_dataset_raw %>%
  group_by(Compound) %>%
  mutate(`Theoretical Amt`=case_when(
    Filename=="2020071205-cal21"~mean(`Theoretical Amt`[Filename=="2020071205-cal22"]),
    TRUE~`Theoretical Amt`
  ))%>%ungroup()

#smiles----

SMILES_data = read_delim("data/Smiles_for_Target_PFAS_semicolon.csv",
                         delim = ";",
                         col_names = TRUE)

SMILES_data = SMILES_data %>%
  rename(Compound = ID) %>%
  select(Compound, SMILES, Class) %>%
  na.omit()

descs_calc_PFAS = PaDEL_original(SMILES_data)

descs_calc_PFAS = read_delim("data/descs_calc.csv",
                             delim = ",",
                             col_names = TRUE)

descs_calc_PFAS = descs_calc_PFAS %>%
  select(Compound, SMILES, 
         #all_of(descs_names),
         everything())

descs_calc_PFAS = descs_calc_PFAS %>%
  group_by(SMILES) %>%
  mutate(IC = isotopedistribution(SMILES),
         MW = molecularmass(SMILES)) %>%
  ungroup()

#eluent---
eluent = read_delim("data/eluent.csv",
                    delim = ",",
                    col_names = TRUE)

organic_modifier = "MeCN"
pH.aq. = 7.0

data = Orbitrap_dataset_raw %>%
  left_join(descs_calc_PFAS)

data = data %>%
  group_by(Compound) %>%
  mutate(RT = mean(RT)) %>%
  ungroup()

data = data %>%
  mutate(
    RT = as.numeric(RT),
    area_IC = Area*IC,
    organic_modifier = organic_modifier,
    pH.aq. = pH.aq.,
    NH4 = 1, #1 if the buffer contains NH¤ ions , 0 if not. 
    organic = organicpercentage(eluent,RT),
    viscosity = viscosity(organic,organic_modifier),
    surface_tension = surfacetension(organic,organic_modifier),
    polarity_index = polarityindex(organic,organic_modifier)) 

training = data %>%
  filter(!is.na(SMILES)) %>%
  filter(`Theoretical Amt` != "N/F") %>%
  filter(`Theoretical Amt` != "N/A") %>%
  mutate(`Theoretical Amt` = as.numeric(`Theoretical Amt`)) %>%
  mutate(`Theoretical Amt` = `Theoretical Amt`/MW) #correct with MW

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

IE_pred = training %>% 
  mutate(logIE_pred = 0) %>%
  rename(organic_modifier_percentage = organic,
         "name" = Compound)%>%
  mutate(additive = "ammoniumacetate",
       additive_concentration_mM = 2,
       instrument = "Orbitrap",
       source = "ESI",
       solvent = "MeCN",#placeholder
       SPLIT = "TRUE")%>%#placeholder
  select(name,SMILES, 1452:1468,everything(),
         -Filename,-RT,-Class,-`Theoretical Amt`,-SPLIT)%>%
  na.omit()

IE_pred = IE_pred%>%
unique()

IE_pred = IE_pred[-c(14,23),]

# prediction =  predict(regressor, newdata = IE_pred, predict.all = TRUE)
# prediction = prediction$aggregate
IE_pred <- IE_pred %>%
  mutate(logIE_pred = predict(regressor, newdata = IE_pred, predict.all = TRUE)) %>%
  select(SMILES,logIE_pred, everything())

IE_pred = IE_pred %>%
  left_join(SMILES_data) %>%
  select(name, SMILES, Class, logIE_pred, slope, everything(),-Compound)

correlation_factor = lm(log10(slope)~logIE_pred,data = IE_pred)#its a model?

IE_pred = IE_pred %>%
mutate(response_factor = correlation_factor$coefficients[2]*logIE_pred +
       correlation_factor$coefficients[1])%>%
  mutate(pred_conc. = area_IC/(10^response_factor))%>%
  left_join(training)%>%
  mutate(error_factor = case_when(pred_conc. < `Theoretical Amt` ~ `Theoretical Amt`/pred_conc.,
                                  pred_conc. >`Theoretical Amt` ~ pred_conc./`Theoretical Amt`))%>%
  mutate(less_than_ten = case_when(error_factor < 10 ~ TRUE,
                                   error_factor > 10 ~ FALSE))

Conc._error_boxplots = ggplot(data = IE_pred) +
  geom_boxplot(mapping = aes(x = Class,
                             y = error_factor))
  
Conc._error_boxplots +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
