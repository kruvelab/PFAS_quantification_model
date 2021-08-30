library(tidyverse)
library(plotly)
setwd("~/GitHub/PFOA_semi_quant")
source("code/PaDEL_descs_calculator.R")
source("code/reading_excel.R")
source("code/compound_eluent.R")

#regressor----

regressor = readRDS("regressors/PFAS_FOREST.rds")

#importing lcms data with background signal----

filename = "data/Batch 1 Semi Quant w frag PL4 rep.xlsx"
Orbitrap_dataset_raw = read_excel_allsheets(filename)

Spiked_samples = Orbitrap_dataset_raw %>%
  filter(Filename == "QCN-CL" | Filename == "QCN-BL" | Filename == "QCN-AL" | Filename == "PL-4") %>%
  select(-`Theoretical Amt`) %>%
  na.omit()

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

#remove adducts + HPFODA
SMILES_data = SMILES_data[-c(17,23,26),]

# descs_calc_PFAS = PaDEL_original(SMILES_data)
# 
# descs_calc_PFAS = read_delim("data/descs_calc.csv",
#                              delim = ",",
#                              col_names = TRUE)
# 
# descs_calc_PFAS = descs_calc_PFAS %>%
#   select(Compound, SMILES, 
#          #all_of(descs_names),
#          everything())
# 
# descs_calc_PFAS = descs_calc_PFAS %>%
#   group_by(SMILES) %>%
#   mutate(IC = isotopedistribution(SMILES),
#          MW = molecularmass(SMILES)) %>%
#   ungroup()

#eluent---
eluent = read_delim("data/eluent.csv",
                    delim = ",",
                    col_names = TRUE)

organic_modifier = "MeCN"
pH.aq. = 7.0

descs <- read_delim("data/descs_recalc.csv")

descs = descs %>%
group_by(SMILES) %>%
mutate(IC = isotopedistribution(SMILES),
MW = molecularmass(SMILES)) %>%
ungroup()

data = Orbitrap_dataset_raw %>%
  left_join(descs)

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
  select(name,SMILES, 1451:1460,everything(),
         -Filename,-RT,)%>%#-`Theoretical Amt`)%>%
  na.omit()

IE_pred = IE_pred%>%
unique()

IE_pred <- IE_pred %>%
  mutate(logIE_pred = predict(regressor, newdata = IE_pred, predict.all = TRUE)) %>%
  select(SMILES,logIE_pred, everything())

IE_pred = IE_pred %>%
  left_join(SMILES_data) %>%
  select(name, SMILES, Class, logIE_pred, slope, everything(),-Compound)

correlation_factor = lm(log10(slope)~logIE_pred,data = IE_pred)

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

IE_pred = IE_pred %>%
  rename("Theoretical_pg/uL" = "Theoretical Amt",
         "pred_pg/uL" = pred_conc.)

concentration_cor_Thomas = ggplot(data = IE_pred %>%
                               filter(instrument == "Orbitrap")) +
  geom_point(mapping = aes(x = "Theoretical_pg/uL",
                           y = "pred_pg/uL", 
                           #text = name)) + 
                           color = name)) +
  scale_y_log10() +
  scale_x_log10() +
  theme(legend.position="none")+
  geom_abline(slope = 1, intercept = 0)+
  facet_wrap(~Class)

concentration_cor_Thomas

#predict concentrations for spiked QC samples

Spiked_samples = Spiked_samples %>%
  rename("Sample ID" = Filename,
         name = Compound)


Spiked_samples = Spiked_samples[-c(36:38,42,68:70,112:114,119:122,59:61),] #all spiked PFAS
#Spiked_samples = Spiked_samples[-c(1:3,7:12,17:25,29:42,46:52,59:73,74:138,59:61),] #only fluorotelomers

IE_prededit <- IE_pred %>%
  select(name, SMILES, logIE_pred,IC,MW)

IE_prededit <- IE_prededit %>%
  unique()

blank_areas = Spiked_samples %>%
  select(name) %>%
  unique() %>%
  left_join(Spiked_samples %>% 
              filter(`Sample ID` == "PL-4") %>%
              select(name, Area)) %>%
  mutate(Area = case_when(
    is.na(Area) ~ 0,
    TRUE ~ Area)) %>%
  rename(Area_background = Area) 

Spiked_samples = Spiked_samples %>%
  left_join(IE_prededit) %>%
  left_join(blank_areas) %>%
  mutate(Area = Area - Area_background,
         area_IC = Area*IC)

Spiked_samples = Spiked_samples %>%
  mutate(response_factor3 = correlation_factor$coefficients[2]*logIE_pred +
           correlation_factor$coefficients[1])%>%
  mutate(pred_conc. = area_IC/(10^response_factor3))

Spiked_samples = Spiked_samples %>%
  mutate(pred_conc_pg_uL = pred_conc.*MW)%>%
  select(pred_conc_pg_uL, everything())%>%  
  mutate(predicted_mass_ng = case_when(
    `Sample ID` == "QCN-AL" ~ (pred_conc_pg_uL*1097.2010178117)/1000,
    `Sample ID` == "QCN-BL" ~ (pred_conc_pg_uL*961.704834605598)/1000,
    `Sample ID` == "QCN-CL" ~ (pred_conc_pg_uL*952.162849872773)/1000,
    TRUE ~ (pred_conc_pg_uL*1000)/1000))%>%
  select(predicted_mass_ng, everything())

ggplot(data = Spiked_samples %>%
         filter(`Sample ID` != "PL-4"), 
       mapping = aes(x = name,
                     y = predicted_mass_ng,
                     fill = `Sample ID`)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_abline(slope = 0, intercept = 5) +
  labs(x = "") +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5,
                                   hjust = 1))



Spiked_samples = Spiked_samples %>%
  filter(Area != 0)%>%
  mutate(error_factor = case_when(predicted_mass_ng < 5 ~ 5/predicted_mass_ng,
                                  predicted_mass_ng >5 ~ predicted_mass_ng/5))%>%
  mutate(less_than_ten = case_when(error_factor < 10 ~ TRUE,
                                   error_factor > 10 ~ FALSE))


mean(Spiked_samples$error_factor, na.rm = TRUE)


QCN_error_boxplots = ggplot(data = IE_pred) +
  geom_boxplot(mapping = aes(x = Class,
                             y = error_factor))

QCN_error_boxplots +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))



#predicting concentrations for detected suspects
suspSMILES_data = read_delim("data/sus_data.csv",
                         delim = ",",
                         col_names = TRUE,)

descs_PFAS_suspects = PaDEL_original(suspSMILES_data %>% select(SMILES) %>% unique())

write_delim(descs_PFAS_suspects,
            "data/descs_PFAS_suspects.csv",
            delim = ",")

suspSMILES_data = suspSMILES_data %>%
  group_by(SMILES) %>%
  mutate(IC = isotopedistribution(SMILES),
         MW = molecularmass(SMILES)) %>%
  ungroup()
  
suspSMILES_data = suspSMILES_data %>%
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

suspSMILES_data =  suspSMILES_data %>%
  left_join(descs_PFAS_suspects, by = "SMILES")%>%
  rename(organic_modifier_percentage = organic)%>%
  mutate(additive = "ammoniumacetate",
         additive_concentration_mM = 2,
         instrument = "Orbitrap",
         source = "ESI",
         solvent = "MeCN",#placeholder
         SPLIT = "TRUE")%>%#placeholder
  select(name,SMILES,everything())%>%
  na.omit()

suspSMILES_data =  suspSMILES_data %>%
  mutate(logIE_pred2 = predict(regressor, newdata = suspSMILES_data, predict.all = TRUE)) %>%
  select(SMILES,logIE_pred2, everything())

suspSMILES_data =  suspSMILES_data %>%
   mutate(response_factor2 = correlation_factor$coefficients[2]*logIE_pred2 +
            correlation_factor$coefficients[1])%>%
   mutate(pred_conc. = area_IC/(10^response_factor2))

suspSMILES_data =  suspSMILES_data %>%
  mutate(pred_conc_pg_uL = pred_conc.*MW.x)%>%
  select(pred_conc_pg_uL, everything())

predicted_conc.= suspSMILES_data %>%
  select(pred_conc_pg_uL, name, `Sample ID`,SMILES,MW.x,`Ex vol.`, `sample weight`)

#convert to F equivalent

#count number of F atoms
x <- suspSMILES_data$SMILES
f.atoms <- lengths(regmatches(x, gregexpr("F",x)))

fatoms. <- data.frame(f.atoms)

blended <- cbind(fatoms., predicted_conc.)

#convert to predicted ng F/g
blended <- blended %>%
  mutate(mass_F = f.atoms*19,
         percent.F=mass_F/MW.x,
         "Predicted ng_F/uL"= (percent.F*pred_conc_pg_uL)/1000,
         pred_ng_g = (`Predicted ng_F/uL`*`Ex vol.`)/`sample weight`)

