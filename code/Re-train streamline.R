library(caret)
library(xgboost)
library(tidyverse)
library(plotly)
setwd("~/GitHub/PFOA_semi_quant")
source("code/PaDEL_descs_calculator.R")
source("code/reading_excel.R")
source("code/compound_eluent.R")

#importing PFAS calibration training data (lcms) ----

filename = "Data/Batch 1 Semi Quant w frag.xlsx"
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
  ))%>%
  ungroup()

#importing Smiles codes for PFAS in calibration set----

 SMILES_data = read_delim("data/Smiles_for_Target_PFAS_semicolon.csv",
                          delim = ";",
                          col_names = TRUE)

 SMILES_data = SMILES_data %>%
   rename(Compound = ID) %>%
   select(Compound, SMILES, Class) %>%
   na.omit()

 SMILES_data = SMILES_data %>%
  group_by(SMILES) %>%
  mutate(IC = isotopedistribution(SMILES),
         MW = molecularmass(SMILES)) %>%
  ungroup()
 
 #removing adducts from dataset
 SMILES_data = SMILES_data[-c(17,26),]

#eluent---
eluent = read_delim("data/eluent.csv",
                    delim = ",",
                    col_names = TRUE)

organic_modifier = "MeCN"
pH.aq. = 7.0

data = Orbitrap_dataset_raw %>%
  left_join(SMILES_data)

#adding eluent data to dataset
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

#getting slopes

training = training %>%
  group_by(SMILES) %>%
  mutate(slope = linear_regression(area_IC, `Theoretical Amt`)$slope) %>%
  ungroup()

#####converting slopes to logIE########

#importing previous trainng set
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

#Anchor IE value
IEPFOSvalue = old_training_data_filtered$logIE

Anchor_slope = training %>% 
  select(Compound,slope) %>%
  unique()%>%
  filter(Compound == "PFOS")

slopePFOS <- Anchor_slope$slope #make slopePFOS a single value

#Converting all slopes to values relative to slope of PFOS,
#converting to log form, adding anchor value
training = training %>%
  mutate(RIE =slope/slopePFOS,
  RIE =log(RIE),
  IE =RIE + IEPFOSvalue)%>%
  select(Compound,3:6, IE,slope,everything(),-RIE)

####prepping + joining old dataset to new

datarbind = training%>%
  select(everything(),-"Theoretical Amt", -Area, -area_IC, -IC, -Filename, -MW, -Class, )%>%
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

######Calculating and joining descriptors for all compounds in training set

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
saveRDS(descs_recalc_go,
        "descs_PFASadd.rds")

datarbindedit = datarbindedit %>%
  select(Compound, SMILES,everything())

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
  select(-c(nearZeroVar(datarbindeditAA, freqCut = 80/20))) 
    
correlationMatrix <- cor(datarbindeditAAA, use = "complete.obs")

highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)

datarbindeditclean <- datarbindeditAAA %>%
  dplyr::select(-highlyCorrelated)%>%
  bind_cols(datarbindedit %>%
              select(SMILES, name, organic_modifier,additive, instrument, source, solvent))%>%
  select(name,SMILES,organic_modifier,organic_modifier_percentage,
         additive, additive_concentration_mM, instrument,source,solvent, everything())

#splitting test and training set

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

#Training model
set.seed(123)
folds = groupKFold(train$name, k = 5)
fitControl <- trainControl(method = "boot", index = folds)

RFR <- 
  train(`logIE`~ ., data = train %>%
          select(-instrument, -source, -split_first,-name,-SMILES),
        method = "xgbTree",
        trControl = fitControl)

saveRDS(RFR,
        "regressors/PFAS_FOREST.rds")

#predicting IE values 
datarbind_with_predicted <- datarbindeditclean %>%
  mutate( logIE_predicted = predict(RFR, newdata = datarbindeditclean))

datarbind_with_predicted <- datarbind_with_predicted %>%
  arrange(instrument)

datarbind_with_predicted$split_first = factor(datarbind_with_predicted$split_first,
                                  levels = c(TRUE, FALSE),
                                  labels = c("Training set", "Test set"))
#plotting predicted vs. calculated IE values for PFAS in calibration set
IE_slope_cor = ggplot() +
  geom_point(data = datarbind_with_predicted %>%
               filter(instrument != "Orbitrap"),
             mapping =aes(logIE, logIE_predicted),
             color = "light grey",
             alpha = 0.5,
             size = 3) +
  geom_point(data = datarbind_with_predicted %>%
               filter(instrument == "Orbitrap"),
             mapping =aes(logIE, logIE_predicted),
             color = "blue",
             alpha = 0.5,
             size = 3) +
  #scale_y_log10() +
  theme(legend.position="none")+
  theme_classic(base_size = 20)+
  geom_abline(slope = 1, intercept = 0)+
  facet_wrap(~split_first)

IE_slope_cor

data_for_split = datarbind_with_predicted %>%
  select(name, split_first, instrument, logIE, logIE_predicted)%>%
  rename(Compound =name)

data43254 = data_for_split %>% filter(instrument == "Orbitrap")

data43254 = data43254%>%
  left_join(SMILES_data)

write.xlsx(data43254,"data/test_train_split.xlsx")

IE_slope_cor_PFAS = ggplot(data = data43254 %>%
                        filter(instrument == "Orbitrap")) +
  geom_point(mapping = aes(x = logIE,y = logIE_predicted,
                           color = Compound,
                           size = 0.25,
                           shape = Class)) +

  #scale_y_log10() +
  theme(legend.position="right")+
  geom_abline(slope = 1, intercept = 0)+
  facet_wrap(~split_first)

IE_slope_cor_PFAS

rmse((datarbind_with_predicted %>% filter(split_first == "Training set"))$logIE,
     (datarbind_with_predicted %>% filter(split_first == "Training set"))$logIE_predicted)
#0.46
rmse((datarbind_with_predicted %>% filter(split_first == "Test set"))$logIE,
     (datarbind_with_predicted %>% filter(split_first == "Test set"))$logIE_predicted)
#0.74


rmse((datarbind_with_predicted %>% filter(instrument == "Orbitrap" & split_first == "Training set"))$logIE,
     (datarbind_with_predicted %>% filter(instrument == "Orbitrap" & split_first == "Training set"))$logIE_predicted)
#0.69
rmse((datarbind_with_predicted %>% filter(instrument == "Orbitrap" & split_first == "Test set"))$logIE,
     (datarbind_with_predicted %>% filter(instrument == "Orbitrap" & split_first == "Test set"))$logIE_predicted)
#0.50

graph_retrainPFAS=ggplotly(IE_slope_cor_PFAS)
graph_retrainPFAS

htmlwidgets::saveWidget(plotly::as_widget(graph_retrainPFAS), "PFAS_IE_predvscalc.html")


