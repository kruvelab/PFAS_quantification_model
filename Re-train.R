library(caret)
library(xgboost)
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

#smiles and descriptors----

 SMILES_data = read_delim("data/Smiles_for_Target_PFAS_semicolon.csv",
                          delim = ";",
                          col_names = TRUE)

 SMILES_data = SMILES_data %>%
   rename(Compound = ID) %>%
   select(Compound, SMILES, Class) %>%
   na.omit()

 descs_calc_PFOA = PaDEL_original(SMILES_data)

 write_delim(descs_calc_PFOA,
             "data/descs_calc.csv",
             delim = ",")

##############################################################3
#descs_calc_PFOA = read_delim("data/descs_calc.csv",
#delim = ",",
#col_names = TRUE)
#############################################################
descs_calc_PFOA = descs_calc_PFOA %>%
  select(Compound, SMILES,
         #all_of(descs_names),
         everything())

descs_calc_PFOA = descs_calc_PFOA %>%
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
  left_join(descs_calc_PFOA)

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

####convert slopes to logIE

filename = "data/IE_training_data/190714_negative_model_logIE_data.csv"
old_training_data = read_delim(filename,
                               delim = ";",
                               col_names = TRUE)

old_training_data = old_training_data%>%
  select(-pH_aq, -logIE_pred, -error_abs)


old_training_data_filtered = old_training_data %>%
  filter(organic_modifier=="MeCN",
  name=="perfluorooctanesulfonic acid",
  additive=="ammonium acetate",
  pH.aq.==7.8)

IEPFOSvalue = old_training_data_filtered$logIE

Anchor_slope = training %>% 
  select(Compound,slope) %>%
  unique()%>%
  filter(Compound == "PFOS")

slopePFOS <- Anchor_slope$slope #make slopePFOS a single value

training = training %>%
  mutate(RIE =slope/slopePFOS,
  RIE =log(RIE),
  IE =RIE + IEPFOSvalue)%>%
  select(Compound,3:6,1452:1463, IE,slope,everything(),-RIE)

####prepping + joining old dataset to new

datarbind = training%>%
  select(1:17,-"Theoretical Amt", -Area, -area_IC, -IC, -Filename)%>%
  unique()%>%
  group_by(Compound) %>%
  mutate(RT=mean(RT))%>%
  ungroup()

datarbind = datarbind%>%
  group_by(Compound)%>%
  mutate(polarity_index=mean(polarity_index),
  organic=mean(organic),
  viscosity=mean(viscosity),
  surface_tension=mean(surface_tension))%>%
  ungroup()%>%
  unique()

datarbindedit = datarbind%>%
  rename(organic_modifier_percentage = organic,
  "name" = Compound)%>%
  select(-RT, -slope)
 
datarbindedit = datarbindedit%>%
  mutate(additive = "ammoniumacetate",
         additive_concentration_mM = 2,
         logIE = IE,
         instrument = "Orbitrap",
         source = "ESI",
         solvent = "MeCN",#placeholder
         SPLIT = "TRUE")%>%#placeholder
  select(-IE)
  
colorder = colnames(datarbindedit)
old_training_data = old_training_data[,colorder]

old_training_data$SPLIT = as.character(old_training_data$SPLIT)

datarbindedit = datarbindedit%>%
  rbind(old_training_data)%>%
  unique()%>%
  select(-SPLIT)

forsplit <- datarbindedit%>%
  select(name) %>%
  unique()

######re-calculating and joining descriptors

descs_recalc = datarbindedit %>%
  rename(Compound = name) %>%
  select(Compound, SMILES) %>%
  unique()%>%
  na.omit()

descs_recalc_go = PaDEL_original(descs_recalc)

datarbindedit = datarbindedit %>%
  left_join(descs_recalc_go)

write_delim(descs_recalc_go,
            "data/descs_recalc.csv",
            delim = ",")

datarbindedit = datarbindedit %>%
  select(Compound, SMILES, 
         #all_of(descs_names),
         everything())

datarbindedit = datarbindedit %>%
  group_by(SMILES) %>%
  mutate(IC = isotopedistribution(SMILES),
         MW = molecularmass(SMILES)) %>%
  ungroup()

#cleaning data

datarbindeditAA = datarbindedit %>%
  dplyr::select(-c(Compound, SMILES, name, organic_modifier,additive, instrument, source, solvent))%>%
  select_if(~sum(is.na(.))< 10,)%>%
  drop_na()
 
datarbindeditAAA = datarbindeditAA %>%
  select(-c(nearZeroVar(datarbindeditAA))) 
    
correlationMatrix <- cor(datarbindeditAAA, use = "complete.obs")

highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)

datarbindeditclean <- datarbindeditAAA %>%
  dplyr::select(-highlyCorrelated)%>%
  dplyr::mutate(name = datarbindedit$name,
                SMILES = datarbindedit$SMILES,
                organic_modifier = datarbindedit$organic_modifier,
                additive = datarbindedit$additive,
                instrument = datarbindedit$instrument,
                source = datarbindedit$source,
                solvent = datarbindedit$solvent)%>%
  select(name,SMILES,organic_modifier,organic_modifier_percentage,
         additive, additive_concentration_mM, instrument,source,solvent, everything())

#splitting

set.seed(123) 
split_first <- sample.split(forsplit$name, SplitRatio = 0.8)
train <- forsplit %>%
  filter(split_first == TRUE)%>%
  mutate(split_first = "TRUE") %>%
  left_join(datarbindeditclean) %>%
  unique()%>%
  na.omit()
test <- forsplit %>%
  filter(split_first == FALSE)%>%
  mutate(split_first = "FALSE") %>%
  left_join(datarbindeditclean) %>%
  unique()%>%
  na.omit()
datarbindeditclean <- rbind(train,test) %>%
  na.omit()

#Training the model
set.seed(123)
folds = groupKFold(train$name, k = 5) #k - how many times 
fitControl <- trainControl(method = "boot", index = folds)


#this fit control is for when there is only one value per compound

# fitControl <- trainControl(
#   method = "repeatedcv",
#   number = 2,  
#   repeats = 1)

RFR <- 
  train(`logIE`~ ., data = datarbindeditclean%>%
          select(-instrument, -source, -split_first,-name),
        method = "xgbTree",
        trControl = fitControl)

datarbind_with_predicted <- datarbindeditclean %>%
  mutate( logIE_pred = predict(RFR, newdata = datarbindeditclean))


#IE_pred = datarbindedit %>% 
  #mutate(logIE_pred = 0) %>%
  #na.omit()

#print(IE_pred %>% select(Compound) %>% unique(), n=40)

#prediction =  predict(RFR, newdata = IE_pred, predict.all = TRUE)
#prediction = prediction$aggregate
# IE_pred <- IE_pred %>%
#   mutate(logIE_pred = prediction) %>%
#   select(SMILES,logIE_pred, everything())

#print(IE_pred %>% select(Compound) %>% unique(), n=40)

# IE_pred = IE_pred %>%
#   left_join(SMILES_data) %>%
#   select(Compound, SMILES, Class, logIE_pred, slope, everything())

#print(IE_pred %>% select(Compound) %>% unique(), n=40)

IE_slope_cor = ggplot(data = datarbind_with_predicted) +
  geom_point(mapping = aes(x = logIE,#see if pred and actual are correlated
                           y = logIE_pred, 
                           #text = name)) + 
                           color = name)) +
  #scale_y_log10() +
  theme(legend.position="none")+
  geom_abline(slope = 1, intercept = 0)+
  facet_wrap(~split_first)

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


