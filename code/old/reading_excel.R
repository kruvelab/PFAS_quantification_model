library(janitor)
library(caret)
library(enviPat)
library(gbm)
library(gsubfn)
library(rJava)
library(rcdk)
library(rcdklibs)
library(tidyverse)
library(OrgMassSpecR)
require(mgcv)
library(stringr)

molecularmass <- function(smiles){
  #convert SMILES to chemical formula
  molecule <- parse.smiles(smiles)[[1]]
  formula <- get.mol2formula(molecule,charge=0)
  formula <- formula@string
  #calcuate molecular weight
  MW <- MolecularWeight(formula = ListFormula(formula))
  return(MW)
}

isotopedistribution <- function(smiles){
  #convert SMILES to chemical formula
  molecule <- parse.smiles(smiles)[[1]]
  formula <- get.mol2formula(molecule,charge=0)
  formula <- formula@string
  
  # Chemical formula to isotope distribution
  data(isotopes)
  pattern<-isopattern(isotopes,
                      formula,
                      threshold=0.1,
                      plotit=FALSE,
                      charge=FALSE,
                      algo=1)
  isotopes <- as.data.frame(pattern[[1]])
  isotope_dist <- as.numeric(sum(isotopes$abundance))
  return(isotope_dist)
}

organicpercentage <- function(eluent_parameters,ret_time){
  ApproxFun <- approxfun(x = eluent_parameters$time, y = eluent_parameters$B)
  organic <- ApproxFun(ret_time)
  return(organic)
}

polarityindex <- function(organic,organic_modifier){
  polarity_index <- case_when(
    organic_modifier == "MeCN" ~ (organic/100)*5.1+((100-organic)/100)*10.2,
    organic_modifier == "MeOH" ~ (organic/100)*5.1+((100-organic)/100)*10.2)
  return(polarity_index)
}


surfacetension <- function(organic,organic_modifier){
  surface_tension <- case_when(
    organic_modifier == "MeCN" ~ 71.76-2.906*71.76*(organic/100)+(7.138*27.86+2.906*71.76-71.76)*(organic/100)^2+(27.86-7.138*27.86)*(organic/100)^3,
    organic_modifier == "MeOH" ~ 71.76-2.245*71.76*(organic/100)+(5.625*22.12+2.245*71.76-71.76)*(organic/100)^2+(22.12-5.625*22.12)*(organic/100)^3)
  return(surface_tension)
}

viscosity <- function(organic,organic_modifier){
  viscosity <- case_when(
    organic_modifier == "MeCN" ~ (-0.000103849885417527)*organic^2+0.00435719229180079*organic+0.884232851261593,
    organic_modifier == "MeOH" ~ (-0.00035908)*organic^2+0.031972067*organic+0.90273943)
  return(viscosity)
}

linear_regression <- function(y, x) {
  if(length(y) > 5) {
    for (i in length(y):5){
      y = y[2:length(y)]
      x = x[2:length(x)]
      slope = summary(lm(y ~ x))$coefficients[2]
      intercept = summary(lm(y ~ x))$coefficients[1]
      residuals = (y - (slope*x +intercept))/y*100
      regression_parameters <- list("slope" = slope, "intercept" = intercept)
      if (max(abs(residuals)) < 10) {
        return(regression_parameters)
        break
      }
    }
    return(regression_parameters)
  }
  else {
    slope = summary(lm(y ~ x))$coefficients[2]
    intercept = summary(lm(y ~ x))$coefficients[1]
    regression_parameters <- list("slope" = slope, "intercept" = intercept)
    return(regression_parameters)
    
  }
}

read_excel_allsheets <- function(filename) {
  data = tibble()
  sheets <- readxl::excel_sheets(filename)
  for (sheet in sheets) {
    data_this_sheet = readxl::read_excel(filename, sheet = sheet, skip = 0)
    data_this_sheet = data_this_sheet %>%
      mutate(across(everything(), as.character))
    #print(data_this_sheet)
    data_this_sheet = data_this_sheet %>%
      group_by(Compound, Filename) %>%
    dplyr::summarise(Area = max(as.double(Area) %>% na.omit()),
             RT = mean(as.double(`Actual RT`) %>% na.omit()),
            `Theoretical Amt` = mean(as.double(str_replace(`Theoretical Amt`, pattern = ",", replacement = ".")) %>% na.omit() )
            ) %>%
      ungroup()
    data = data %>%
      bind_rows(data_this_sheet)
  }
  data <- data %>%
    rename(Theoretical_amt = `Theoretical Amt`) %>%
    filter(Area != "-Inf") %>%                          
    mutate(Area = as.numeric(Area))
  
  return(data)
}

read_SMILES <- function(filename, 
                        compounds_to_be_removed_as_list = c()
                        ) {
  
  SMILES_data = read_delim(filename,
                            delim = ";",
                            col_names = TRUE)
  
  SMILES_data = SMILES_data %>%
    rename(Compound = ID) %>%
    select(Compound, SMILES) %>%
    na.omit() %>%
    group_by(SMILES) %>%
    mutate(IC = isotopedistribution(SMILES),
           MW = molecularmass(SMILES)) %>%
    ungroup() 
  if (length(compounds_to_be_removed_as_list) > 0) {
    SMILES_data = SMILES_data%>%
      filter(!Compound %in% compounds_to_be_removed_as_list)
  }
  return(SMILES_data)
}

add_mobile_phase_composition = function(data,
                                    eluent_file_name,
                                    organic_modifier = "MeCN",
                                    pH.aq. = 7.0,
                                    NH4 = 1,
                                    additive = "ammoniumacetate",
                                    additive_concentration_mM = 2,
                                    instrument = "Orbitrap",
                                    source = "ESI") {
  eluent = read_delim(eluent_file_name,
                      delim = ",",
                      col_names = TRUE)
  
  ## Joining all collected data to one tibble, removing missing values, calculating slopes ----
  data = data %>%
    mutate(organic_modifier_percentage = organicpercentage(eluent,RT),
           viscosity = viscosity(organic_modifier_percentage,organic_modifier),
           surface_tension = surfacetension(organic_modifier_percentage,organic_modifier),
           polarity_index = polarityindex(organic_modifier_percentage,organic_modifier), 
           organic_modifier = organic_modifier,
           pH.aq. = pH.aq.,
           NH4 = 1,
           additive = "ammoniumacetate",
           additive_concentration_mM = 2,
           instrument = "Orbitrap",
           source = "ESI")
  return(data)
}

anchoring = function(data_to_be_anchored,
                     data_containing_anchor,
                     binding = TRUE,
                     anchor = "perfluorooctanesulfonic acid",
                     anchor_in_new_dataset = "PFOS",
                     organic_modifier = "MeCN",
                     additive = "ammonium acetate",
                     pH.aq. = 7.8) {
  #training = data
  old_training_data = read_delim(data_containing_anchor,
                                 delim = ";",
                                 col_names = TRUE)
  
  # Anchor compound value from previous data set
  old_training_data_IEPFOSvalue = old_training_data %>%
    filter(organic_modifier == organic_modifier,
           name == anchor,
           additive == additive,
           pH.aq. == pH.aq.)
  
  # Anchor compound value from this dataset
  Anchor_slope = data_to_be_anchored %>% 
    select(Compound, slope) %>%
    unique()%>%
    filter(Compound == anchor_in_new_dataset)
  
  data_to_be_anchored = data_to_be_anchored %>%
    mutate(logRIE = log(slope/Anchor_slope$slope),
           logIE  = logRIE + old_training_data_IEPFOSvalue$logIE) %>% 
    select(Compound, SMILES, logIE, pH.aq., polarity_index, viscosity, surface_tension, NH4)
  
  ## Prepping + joining old dataset to new ----
  if (binding) {
    data_to_be_anchored = data_to_be_anchored %>%
      rename(name = Compound) %>%
      mutate(data_type = "new") %>%
      bind_rows(old_training_data %>%
                  select(name, SMILES, logIE, pH.aq., polarity_index, viscosity, surface_tension, NH4) %>%
                  mutate(data_type = "old")
      )
  }
  return(data_to_be_anchored)
}

cleaning_data = function(data,
                         nearZeroVar_freqCut = 80/20,
                         highlyCorrelated_cutoff = 0.75) {
  data = data %>%
    drop_na()
  name = data %>%
    select(name)
  data_type = data %>%
    select(data_type)
  # Removing columns with missing values
  data = data %>%
    dplyr::select(-c(SMILES, name, data_type)) %>%
    select_if(~ sum(is.na(.))< 10,) 
  
  # Checking that any of the categorical values would not have more than 80% of existing/non-existing values
  data = data %>%  
    select(-c(nearZeroVar(data, 
                          freqCut = nearZeroVar_freqCut)))
  
  # Removing columns that are correlated to each other
  correlationMatrix <- cor(data, use = "complete.obs")
  highlyCorrelated <- findCorrelation(correlationMatrix, 
                                      cutoff = highlyCorrelated_cutoff)
  
  data <- data %>%
    dplyr::select(-highlyCorrelated) %>%
    bind_cols(name) %>%
    bind_cols(data_type)
  
  return(data)
}

training_logIE_pred_model = function(data,
                                     folds = 5,
                                     fitControlMethod = "boot",
                                     method = "xgbTree",
                                     split = NULL,
                                     save_model_name = NULL) {
  set.seed(123)
  if (!is.null(split)) {
    training_set = tibble()
    test_set = tibble()
    for (data_type_this in levels(factor(data$data_type))) {
      data_this = data %>%
        filter(data_type == data_type_this)
      name = data_this %>% select(name) %>% unique()
      split_train_test = sample.split(name$name, SplitRatio = split)
      name = name %>% mutate(split_train_test = split_train_test)
      data_this = data_this %>%
        left_join(name)
      training_set = training_set %>%
        bind_rows(data_this %>% 
                    filter(split_train_test))
      test_set = test_set %>%
        bind_rows(data_this %>% 
                    filter(!split_train_test))
    }
  } else {
    training_set = data
  }
  
  folds = groupKFold(training_set$name, k = folds) 
  fitControl <- trainControl(method = fitControlMethod, index = folds)
  
  model <- train(logIE ~ ., 
               data = training_set %>% select(-name, -data_type),
               method = method,
               trControl = fitControl)
  training_set <- training_set %>%
    mutate(logIE_predicted = predict(model, newdata = training_set))
  
  RMSE_training_set = rmse(training_set$logIE, training_set$logIE_predicted)
  
  ## Determining most influential descriptors ----
  variable_importance <- varImp(model)
  variable_importance <- as.data.frame(variable_importance$importance)
  
  data = list("training_set" = training_set)
  metrics = list("RMSE_training_set" = RMSE_training_set)
  
  if (!is.null(save_model_name)) {
    saveRDS(model,save_model_name)
  }
  
  if (!is.null(split)) {
    test_set <- test_set %>%
      mutate(logIE_predicted = predict(model, newdata = test_set))
    RMSE_test_set = rmse(test_set$logIE, test_set$logIE_predicted)
    
    data = list("training_set" = training_set,
                "test_set" = test_set)
    metrics = list("RMSE_training_set" = RMSE_training_set,
                   "RMSE_test_set" = RMSE_test_set)
  }
  
  model = list("model" = model,
               "data" = data,
               "metrics" = metrics,
               "variable_importance" = variable_importance)
  return(model)
  
}

