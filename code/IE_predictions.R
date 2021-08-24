library(tidyverse)
library(plotly)
setwd("~/GitHub/PFOA_semi_quant")
source("code/PaDEL_descs_calculator.R")
source("code/reading_excel.R")
source("code/compound_eluent.R")

#loading regressor----

regressor = readRDS("regressors/PFAS_FOREST.rds")

#lcms data from model training----

training = read_delim("data/cal_exp_data.csv",
           delim = ",",
           col_names = TRUE)

IE_pred <- training

#predicting log_IE for cal. standards

IE_pred <- IE_pred %>%
  mutate(logIE_pred = predict(regressor, newdata = IE_pred, predict.all = TRUE)) %>%
  select(SMILES,logIE_pred, everything())

slopes = read_delim("data/slopes.csv",
                    delim = ",",
                    col_names = TRUE)

correlation_factor = lm(log10(slope)~logIE_pred,data = IE_pred)

#Comparing predicted concentrations and theoretical concentrations

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

filename = "data/Batch 1 Semi Quant w frag PL4 rep.xlsx"
Orbitrap_dataset_raw = read_excel_allsheets(filename)

Spiked_samples = Orbitrap_dataset_raw %>%
  filter(Filename == "QCN-CL" | Filename == "QCN-BL" | Filename == "QCN-AL" | Filename == "PL-4") %>%
  select(-`Theoretical Amt`) %>%
  na.omit()

Spiked_samples = Spiked_samples %>%
  rename("Sample ID" = Filename,
         name = Compound)

#Spiked_samples = Spiked_samples[-c(1:3,7:12,17:25,29:42,46:52,59:73,127:138,112:114,119:122),] #structurally relevant
Spiked_samples = Spiked_samples[-c(36:38,42,68:70,112:114,119:122),] #all spiked PFAS
#Spiked_samples = Spiked_samples[-c(1:3,7:12,17:25,29:42,46:52,59:73,74:138),] #only fluorotelomers

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
suspSMILES_data = read_delim("data/sus data.csv",
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
  select(pred_conc_pg_uL, name, `Sample ID`)

write.xlsx(predicted_conc.,"data/pred_conc_sus.xlsx")


#convert to F equivalent (in excel)

#dataset[Filename == "Spiked sample"]$Area - dataset[Filename == "Blank"]$Area 




