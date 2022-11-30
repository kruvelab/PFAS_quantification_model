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
library(xgboost)
library(Metrics)
library(readxl)
library(plotly)

molecularmass <- function(smiles){
  #convert SMILES to chemical formula
  molecule <- parse.smiles(smiles)[[1]]
  formula <- get.mol2formula(molecule,charge=0)
  formula <- formula@string
  #calcuate molecular weight
  MW <- MolecularWeight(formula = ListFormula(formula))
  return(MW)
}

smiles_to_formula <- function(smiles){
  #convert SMILES to chemical formula
  molecule <- parse.smiles(smiles)[[1]]
  formula <- get.mol2formula(molecule,charge=0)
  formula <- formula@string
  return(formula)
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
  df <- tibble(x = x, #log(x), 
               y = y, #log(y)
  ) %>%
    group_by(x) %>% 
    mutate(y = mean(y)) %>% 
    ungroup() %>% 
    unique() %>% 
    na.omit()
  y = df$y
  x = df$x
  if(length(x) == 0){
    slope = NA
    intercept = NA
    regression_parameters <- list("slope" = slope, "intercept" = intercept)
    return(regression_parameters)
  } else if (length(y) > 5) {
    for (i in length(y):5){
      y = y[2:(length(y))]
      x = x[2:(length(x))]
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
  } else {
    slope = summary(lm(y ~ x))$coefficients[2]
    intercept = summary(lm(y ~ x))$coefficients[1]
    regression_parameters <- list("slope" = slope, "intercept" = intercept)
    return(regression_parameters)
  }
}


# linear_regression <- function(y, x) {
#   df <- tibble(x = x, #log(x), 
#                y = y, #log(y)
#                ) %>%
#     group_by(x) %>% 
#     mutate(y = mean(y)) %>% 
#     ungroup() %>% 
#     unique() %>% 
#     na.omit()
#    y = df$y
#    x = df$x
#   
#   
#   if(length(x) == 0){
#     slope = NA
#     intercept = NA
#     regression_parameters <- list("slope" = slope, "intercept" = intercept)
#     return(regression_parameters)
#   } else if (length(y) > 5) {
#     for (i in length(y):5){
#       y = y[1:(length(y)-1)]
#       x = x[1:(length(x)-1)]
#       slope = summary(lm(y ~ x))$coefficients[2]
#       intercept = summary(lm(y ~ x))$coefficients[1]
#       residuals = (y - (slope*x +intercept))/y*100
#       regression_parameters <- list("slope" = slope, "intercept" = intercept)
#       if (max(abs(residuals)) < 50) {
#         return(regression_parameters)
#         break
#       }
#     }
#     return(regression_parameters)
#   } else {
#     slope = summary(lm(y ~ x))$coefficients[2]
#     intercept = summary(lm(y ~ x))$coefficients[1]
#     regression_parameters <- list("slope" = slope, "intercept" = intercept)
#     return(regression_parameters)
#   }
# }

# model <- lm(y ~ x)
# res <- resid(model)
# plot(fitted(model), res)
# plot(fitted(model), residuals)


read_excel_allsheets <- function(filename) {
  data = tibble()
  sheets <- readxl::excel_sheets(filename)
  for (sheet in sheets) {
    data_this_sheet = readxl::read_excel(filename, sheet = sheet, skip = 0)
    data_this_sheet = data_this_sheet %>%
      mutate(across(everything(), as.character))
    #print(data_this_sheet)
    #suppressMessages(
      data_this_sheet = data_this_sheet %>%
        group_by(Compound, Filename) %>%
        dplyr::summarise(Area = max(as.double(Area) %>% na.omit()),
                         RT = mean(as.double(`Actual RT`) %>% na.omit()),
                         `Theoretical Amt` = mean(as.double(str_replace(`Theoretical Amt`, pattern = ",", replacement = ".")) %>% na.omit())) %>%
        ungroup()
     # )
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
  
  SMILES_data = read.csv(filename, sep=";")
  
  SMILES_data = SMILES_data %>%
    rename(Compound = ID) %>%
    select(Compound, SMILES) %>%
    na.omit() %>%
    group_by(SMILES) %>%
    mutate(IC = isotopedistribution(SMILES),
           Molecular_weight = molecularmass(SMILES)) %>%
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
           NH4 = NH4,
           additive = additive,
           additive_concentration_mM = additive_concentration_mM,
           instrument = instrument,
           source = source)
  return(data)
}


PaDEL_original = function(standards) {
  SMILES_list = standards %>% 
    select(SMILES) %>% 
    unique() %>%
    na.omit() %>%
    mutate(Name = row_number())
  standards = standards %>%
    left_join(SMILES_list)
  
  write_delim(SMILES_list %>% select(SMILES) %>% unique(),
              "SMILES.smi",
              col_names = FALSE)
  
  command = "java -jar PaDEL-Descriptor/PaDEL-Descriptor.jar -dir" #file name where they will be calculated
  command_final = paste(command, "SMILES.smi", "-file", "descs_calc.csv", "-2d", sep =" ") #makes text for command prompt
  javaOutput = system(command_final, intern = TRUE) #goes into commant prompt
  #PaDEL saves the descriptors to the local folder
  descs = read_delim("descs_calc.csv",
                     delim = ",",
                     col_names = TRUE)
  
  descs = descs %>%
    group_by(Name) %>%
    mutate(Name = str_split(Name, pattern = "_")[[1]][length(str_split(Name, pattern = "_")[[1]])]) %>%
    ungroup() %>%
    mutate(Name = as.integer(Name))
  
  descs = descs %>%
    left_join(standards) %>%
    select(colnames(standards), everything()) %>%
    select(-Name)
  
  write_delim(descs,
              "data/descs_calc.csv",
              delim = ",")
  
  return(descs)
}

anchoring = function(data_to_be_anchored,
                     data_containing_anchor,
                     binding = TRUE,
                     anchor = "perfluorooctanesulfonic acid",
                     anchor_in_new_dataset = "PFOS",
                     organic_modifier = "MeCN",
                     additive_name = "ammonium acetate",
                     pH_aqueous = 7.8) {
  #training = data
  old_training_data = read_delim(data_containing_anchor,
                                 delim = ";",
                                 col_names = TRUE)
  
  # Anchor compound value from previous data set
  old_training_data_IEPFOSvalue = old_training_data %>%
    filter(organic_modifier == organic_modifier,
           name == anchor,
           additive == additive_name,
           pH.aq. == pH_aqueous)
  
  # Anchor compound value from this dataset
  Anchor_slope = data_to_be_anchored %>% 
    select(Compound, slope) %>%
    unique()%>%
    filter(Compound == anchor_in_new_dataset)
  
  data_to_be_anchored = data_to_be_anchored %>%
    mutate(#logRIE = log(slope/Anchor_slope$slope),
           logIE  = log(slope/Anchor_slope$slope) #logRIE 
                    + old_training_data_IEPFOSvalue$logIE) %>% 
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
  smiles = data %>%
    select(SMILES)
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
    bind_cols(data_type) %>%
    bind_cols(smiles)
  
  return(data)
}

training_logIE_pred_model = function(data,
                                     folds = 5,
                                     fitControlMethod = "boot",
                                     method = "xgbTree",
                                     split = NULL,
                                     save_model_name = NULL) {
  #data = data_clean
  #split = 0.8
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
                    filter(split_train_test)) %>%
        select(-split_train_test)
      test_set = test_set %>%
        bind_rows(data_this %>% 
                    filter(!split_train_test)) %>%
        select(-split_train_test)
    }
  } else {
    training_set = data
  }
  
  folds = groupKFold(training_set$name, k = folds) 
  fitControl <- trainControl(method = fitControlMethod, index = folds)
  
  model <- train(logIE ~ ., 
               data = training_set %>% select(-name, -data_type, -SMILES),
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

#NB works until count for each atom is under 100!!!
formula_to_atomCount <- function(formula) {
  
  split_string <- (str_split(formula, ""))[[1]]
  atoms <- tibble()
  counts <- tibble()
  
  for (i in 1:length(split_string)) {
    atom <- ""
    number_count <- ""
    if (is.na(as.double(split_string[i])) == TRUE && toupper(split_string[i]) == split_string[i]) { #Letter
      atom <- split_string[i]
      if (i<length(split_string)) {
        if (is.na(as.double(split_string[i+1])) == TRUE) { #Next one also letter
          if(toupper(split_string[i+1]) != split_string[i+1]){
            atom <- paste(atom, split_string[i+1], sep="")
            i <- i+1
          }
          else {
            number_count <- "1"
          }
        }
        
        if (is.na(as.double(split_string[i+1])) == FALSE) {
          number_count <- split_string[i+1]
          if(i+1 < length(split_string) && (is.na(as.double(split_string[i+2])) == FALSE)) {
            number_count <- paste(number_count, split_string[i+2], sep = "")
          }
        }
        
      }
      else {
        number_count <- "1"
      }
      atoms <- c(atoms, atom)
      counts <- c(counts, as.double(number_count)) 
    }
  }
  
  for (i in 1:length(atoms)){
    atoms[[i]] <- paste("atom", atoms[[i]], sep = "_")
  }
  
  summary <- as.data.frame(counts, col.names = atoms)
  return(summary)
}

atomcount_columns_toDF <- function(dataframe) {
  Formulas <- dataframe %>%
    select(molecularFormula) %>%
    unique()
  
  atomCoun_all <- tibble()
  for (i in 1:length(Formulas$molecularFormula)) {
    atomCount <- formula_to_atomCount(Formulas$molecularFormula[i])
    atomCoun_all <- atomCoun_all %>%
      bind_rows(atomCount)
  }
  
  Formulas <- Formulas %>%  
    bind_cols(atomCoun_all)
  
  dataframe <- dataframe %>%
    left_join(Formulas)
  
  return(dataframe)
}


is_homologue = function(analyte_SMILES,
                        standard_SMILES,
                        homologue_pattern_SMART) {
  analyte = get.mol2formula(parse.smiles(analyte_SMILES)[[1]])@isotopes
  standard = get.mol2formula(parse.smiles(standard_SMILES)[[1]])@isotopes
  pattern = get.mol2formula(parse.smiles(homologue_pattern_SMART)[[1]])@isotopes
  homologue_of_analyte = as_tibble(analyte) %>%
    left_join(as_tibble(pattern), by = c("isoto" = "isoto")) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    mutate(smaller = as.numeric(number.x) - as.numeric(number.y),
           bigger = as.numeric(number.x) + as.numeric(number.y)) %>%
    select(isoto, smaller, bigger)
  comparison = homologue_of_analyte %>%
    left_join(as_tibble(standard) %>%
                mutate(number = as.numeric(number)),
              by = c("isoto" = "isoto")) %>%
    mutate_all(~replace(., is.na(.), 0))
  if(all(comparison$smaller == comparison$number)) {
      return("smaller")
  } else if  (all(comparison$bigger == comparison$number)) {
      return("bigger")
  } else {
      return(NA)
  }
}


concentration_forAnalytes_model <- function(filename_data, 
                                            filename_smiles, 
                                            filename_eluent, 
                                            pred_model,
                                            compounds_to_be_removed_as_list = c(),
                                            organic_modifier = "MeCN",
                                            pH.aq. = 7.0,
                                            NH4 = 1,
                                            additive = "ammoniumacetate",
                                            additive_concentration_mM = 2,
                                            instrument = "Orbitrap",
                                            source = "ESI") {
  
  analysis_data <- read_excel_allsheets(filename_data)
  SMILES_data <- read_SMILES(filename_smiles, compounds_to_be_removed_as_list)
  
  analysis_data = analysis_data %>%
    mutate(Theoretical_amt = replace(Theoretical_amt , grepl("NaN", Theoretical_amt, fixed = TRUE), NA)) %>%
    left_join(SMILES_data) %>%
    drop_na(SMILES) %>%
    mutate(RT = as.numeric(RT),
           area_IC = Area*IC,
           Theoretical_conc_uM = Theoretical_amt/Molecular_weight) 
  
  analysis_data_descr <- analysis_data %>%    
    group_by(SMILES, Compound) %>%
    mutate(slope = linear_regression(area_IC, Theoretical_conc_uM)$slope,
           RT = mean(RT)) %>%
    ungroup()
  
  analysis_data_descr = add_mobile_phase_composition(data = analysis_data_descr,
                                                     eluent_file_name = filename_eluent,
                                                     organic_modifier = organic_modifier,
                                                     pH.aq. = pH.aq.,
                                                     NH4 = NH4,
                                                     additive = additive,
                                                     additive_concentration_mM = additive_concentration_mM,
                                                     instrument = instrument,
                                                     source = source)
  
  
  
  analysis_data_descr <- PaDEL_original(analysis_data_descr)
  analysis_data_descr <- analysis_data_descr %>%
    mutate(logIE_predicted = predict(pred_model$model, newdata = analysis_data_descr))
  
  # lm function to find RF of suspect compound
  linMod <- lm(log10(slope) ~ logIE_predicted, data = analysis_data_descr %>%
                                                            drop_na(slope))     # Take this out in final!!
  
  
  # Plot measured vs predicted
  plot_predictedIE_slope = ggplot() +
    geom_point(data = analysis_data_descr,
               mapping = aes(logIE_predicted, log10(slope),
                             text = Compound),
               color = "black",
               alpha = 0.5,
               size = 3) +
    geom_abline(slope = summary(linMod)$coefficients[2], intercept = summary(linMod)$coefficients[1]) 
  ggplotly(plot_predictedIE_slope)
  
  
  # Find RF values from predicted IEs to all non-calibrant analytes
  analytes_concentrations <- analysis_data_descr %>%
    # filter(is.na(slope)) %>%
    mutate(slope_pred = 10^(summary(linMod)$coefficients[2]*logIE_predicted + summary(linMod)$coefficients[1])) %>%
    mutate(conc_pred_uM = area_IC/slope_pred) %>%
    mutate(conc_pred_pg_uL = conc_pred_uM*Molecular_weight) %>% 
    select(colnames(analysis_data), slope_pred, conc_pred_uM, conc_pred_pg_uL)
  
  plot_predicted_theoretical_conc <- ggplot() +
    geom_point(data = analytes_concentrations %>%
                  drop_na(Theoretical_conc_uM) 
               #   group_by(Compound, Theoretical_conc_uM) %>%
               #   mutate(conc_pred_uM = mean(conc_pred_uM)) %>%
               #   ungroup()
               ,
               mapping = aes(Theoretical_conc_uM, conc_pred_uM,
                             color = Compound)) +
    geom_abline(slope = 1, intercept = 0) +
    scale_y_log10(limits = c(10^-5, 10^0)) +
    scale_x_log10(limits = c(10^-5, 10^0)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_abline(slope = 1, intercept = 1) +
    geom_abline(slope = 1, intercept = -1) +
    theme(aspect.ratio = 1)
  ggplotly(plot_predicted_theoretical_conc)
  
  
  #Return list with data, plot
  data_predicted = list("plot_predictedIE_slope" = plot_predictedIE_slope,
                        "plot_predicted_theoretical_conc" = plot_predicted_theoretical_conc,
                        "data" = analytes_concentrations)
  
  return(data_predicted)
}

## Function to convert concentration to M conc?



concentration_forAnalytes_homolog <- function(filename_data, 
                                              filename_smiles,
                                              homolog_pattern_SMILES,
                                              findHomolog_onlyForAnalytes = TRUE) {
  
  analysis_data <- read_excel_allsheets(filename_data)
  SMILES_data <- read_SMILES(filename_smiles)
  
  analysis_data = analysis_data %>%
    left_join(SMILES_data) %>%
    mutate(Theoretical_amt = replace(Theoretical_amt , grepl("NaN", Theoretical_amt, fixed = TRUE), NA)) %>%
    drop_na(SMILES) %>%
    mutate(RT = as.numeric(RT),
           area_IC = Area*IC,
           Theoretical_conc_uM = Theoretical_amt/Molecular_weight) 
  
  analysis_data_descr <- analysis_data %>%    
    group_by(SMILES, Compound) %>%
    mutate(slope = linear_regression(area_IC, Theoretical_conc_uM)$slope,
           RT = mean(RT)) %>%
    ungroup()
  
  if(findHomolog_onlyForAnalytes == TRUE) {
    # Create two datasets of analytes and calibrants - quantification based on homologous calibrants used in this analysis?
    data_training_original <- analysis_data_descr %>%
      filter(!is.na(slope)) %>%
      select(Compound, SMILES, slope) %>%
      unique()
    
    data_analytes <- analysis_data_descr %>%
      filter(is.na(slope)) %>%
      select(Compound, Filename, SMILES, area_IC) %>%
      unique()
    
    data_joined <- data_training_original %>%
      left_join(data_analytes, by = character())
    
  } else if (findHomolog_onlyForAnalytes == FALSE) {
    data_joined <- analysis_data_descr %>%
      select(Compound, SMILES, slope) %>%
      unique() %>%
      left_join(analysis_data_descr %>%
                  select(Compound, Filename, SMILES, area_IC) %>%
                  unique(), 
                by = character())
  }
  
  
  # data_joined <- data_joined %>%
  #   group_by(SMILES.x, SMILES.y) %>%
  #   mutate(pattern = is_homologue(SMILES.x, SMILES.y, homolog_pattern_SMILES)) %>% # Problems with isomeric compounds???  "F[C+2]F"
  #   ungroup() %>%
  #   mutate(conc = case_when((pattern == "smaller" | pattern == "bigger") ~ area_IC/slope))
  
  data_joined_smiles <- data_joined %>%
    select(SMILES.x, SMILES.y) %>%
    unique() %>%
    group_by(SMILES.x, SMILES.y) %>%
    mutate(pattern = is_homologue(SMILES.y, SMILES.x, homolog_pattern_SMILES)) %>%
    ungroup()
  
  data_joined <- data_joined %>%
    left_join(data_joined_smiles) %>%
    mutate(conc = case_when((pattern == "smaller" | pattern == "bigger") ~ (area_IC)/slope)) %>%
    filter(!is.na(pattern))
  
  return(data_joined)
}


concentration_forAnalytes_homolog_withIntercept <- function(filename_data, 
                                                            filename_smiles,
                                                            homolog_pattern_SMILES,
                                                            findHomolog_onlyForAnalytes = TRUE) {
  
  analysis_data <- read_excel_allsheets(filename_data)
  SMILES_data <- read_SMILES(filename_smiles)
  
  analysis_data = analysis_data %>%
    left_join(SMILES_data) %>%
    mutate(Theoretical_amt = replace(Theoretical_amt , grepl("NaN", Theoretical_amt, fixed = TRUE), NA)) %>%
    drop_na(SMILES) %>%
    mutate(RT = as.numeric(RT),
           area_IC = Area*IC,
           Theoretical_conc_uM = Theoretical_amt/Molecular_weight) 
  
  analysis_data_descr <- analysis_data %>%    
    group_by(SMILES, Compound) %>%
    mutate(slope = linear_regression(area_IC, Theoretical_conc_uM)$slope,
           intercept = linear_regression(area_IC, Theoretical_conc_uM)$intercept,
           RT = mean(RT)) %>%
    ungroup()
  
  if(findHomolog_onlyForAnalytes == TRUE) {
    # Create two datasets of analytes and calibrants - quantification based on homologous calibrants used in this analysis?
    data_training_original <- analysis_data_descr %>%
      filter(!is.na(slope)) %>%
      select(Compound, SMILES, slope, intercept) %>%
      unique()
    
    data_analytes <- analysis_data_descr %>%
      filter(is.na(slope)) %>%
      select(Compound, Filename, SMILES, area_IC) %>%
      unique()
    
    data_joined <- data_training_original %>%
      left_join(data_analytes, by = character())
    
  } else if (findHomolog_onlyForAnalytes == FALSE) {
    data_joined <- analysis_data_descr %>%
      select(Compound, SMILES, slope, intercept) %>%
      unique() %>%
      left_join(analysis_data_descr %>%
                  select(Compound, Filename, SMILES, area_IC) %>%
                  unique(), 
                by = character())
  }
  
  
  # data_joined <- data_joined %>%
  #   group_by(SMILES.x, SMILES.y) %>%
  #   mutate(pattern = is_homologue(SMILES.x, SMILES.y, homolog_pattern_SMILES)) %>% # Problems with isomeric compounds???  "F[C+2]F"
  #   ungroup() %>%
  #   mutate(conc = case_when((pattern == "smaller" | pattern == "bigger") ~ area_IC/slope))
  
  data_joined_smiles <- data_joined %>%
    select(SMILES.x, SMILES.y) %>%
    unique() %>%
    group_by(SMILES.x, SMILES.y) %>%
    mutate(pattern = is_homologue(SMILES.y, SMILES.x, homolog_pattern_SMILES)) %>%
    ungroup()
  
  data_joined <- data_joined %>%
    left_join(data_joined_smiles) %>%
    mutate(conc = case_when((pattern == "smaller" | pattern == "bigger") ~ (area_IC-intercept)/slope)) %>%
    filter(!is.na(pattern))
  
  return(data_joined)
}


concentration_forAnalytes_model_cal_separateFile <- function(cal_filename_data, 
                                            cal_filename_smiles, 
                                            sus_filename_data,
                                            sus_filename_smiles,
                                            filename_eluent, 
                                            pred_model,
                                            compounds_to_be_removed_as_list = c(),
                                            organic_modifier = "MeCN",
                                            pH.aq. = 7.0,
                                            NH4 = 1,
                                            additive = "ammoniumacetate",
                                            additive_concentration_mM = 2,
                                            instrument = "Orbitrap",
                                            source = "ESI") {
  
  analysis_data_cal <- read_excel_allsheets(cal_filename_data)
  SMILES_data_cal <- read_SMILES(cal_filename_smiles, compounds_to_be_removed_as_list)
  #suspects data from analysis, filtered by smiles
  analysis_data_sus <- read_excel_allsheets(sus_filename_data)
  SMILES_data_sus <- read_SMILES(sus_filename_smiles)
  analysis_data_sus <- analysis_data_sus %>%
    mutate(Theoretical_amt = replace(Theoretical_amt , grepl("NaN", Theoretical_amt, fixed = TRUE), NA)) %>%
    left_join(SMILES_data_sus) %>%
    drop_na(SMILES)
  
  
  analysis_data = analysis_data_cal %>%
    mutate(Theoretical_amt = replace(Theoretical_amt , grepl("NaN", Theoretical_amt, fixed = TRUE), NA)) %>%
    left_join(SMILES_data_cal) %>%
    drop_na(SMILES) %>%
    bind_rows(analysis_data_sus) %>%
    mutate(RT = as.numeric(RT),
           area_IC = Area*IC,
           Theoretical_conc_uM = Theoretical_amt/Molecular_weight) 
  
  analysis_data_descr <- analysis_data %>% 
    #filter(Compound == "11Cl-PF3OUdS") %>% 
    group_by(Compound) %>%
    mutate(slope = linear_regression(area_IC, Theoretical_conc_uM)$slope,
           intercept = linear_regression(area_IC, Theoretical_conc_uM)$intercept,
           RT = mean(RT)) %>%
    ungroup()
  
  # Calibrants plot

  plot_calgraph <- ggplot() +
    geom_point(data = analysis_data_descr,
               mapping = aes(x = Theoretical_conc_uM,
                             y = area_IC)) +
  facet_wrap(~ Compound, scales = "free") +
    theme_classic()
  # plot_calgraph
  
  analysis_data_descr = add_mobile_phase_composition(data = analysis_data_descr,
                                                     eluent_file_name = filename_eluent,
                                                     organic_modifier = organic_modifier,
                                                     pH.aq. = pH.aq.,
                                                     NH4 = NH4,
                                                     additive = additive,
                                                     additive_concentration_mM = additive_concentration_mM,
                                                     instrument = instrument,
                                                     source = source)
  
  
  
  analysis_data_descr <- PaDEL_original(analysis_data_descr)
  
  
  analysis_data_descr <- analysis_data_descr %>%
    mutate(logIE_predicted = predict(pred_model$model, newdata = analysis_data_descr))
  
  # lm function to find RF of suspect compound
  linMod <- lm(logIE_predicted ~ log10(slope), data = analysis_data_descr %>%
                 drop_na(slope))
  
  # Plot measured vs predicted
  plot_predictedIE_slope = ggplot() +
    geom_point(data = analysis_data_descr,
               mapping = aes(log10(slope),
                             logIE_predicted,
                             text = Compound),
               color = "black",
               alpha = 0.5,
               size = 3) +
    theme_classic() +
    geom_abline(slope = summary(linMod)$coefficients[2], intercept = summary(linMod)$coefficients[1]) +
    geom_abline(slope = summary(linMod)$coefficients[2], intercept = 1+summary(linMod)$coefficients[1]) +
    geom_abline(slope = summary(linMod)$coefficients[2], intercept = -1+summary(linMod)$coefficients[1]) +
    xlab("log10(measured RF)") +
    ylab("log10(predicted IE)")
  ggplotly(plot_predictedIE_slope)
  

  # Find RF values from predicted IEs to all non-calibrant analytes
  analytes_concentrations <- analysis_data_descr %>%
    # filter(is.na(slope)) %>%
    mutate(slope_pred = 10^((logIE_predicted - summary(linMod)$coefficients[1])/summary(linMod)$coefficients[2])) %>%
    mutate(conc_pred_uM = area_IC/slope_pred) %>%
    mutate(conc_pred_pg_uL = conc_pred_uM*Molecular_weight) %>% 
    select(colnames(analysis_data), slope_pred, conc_pred_uM, conc_pred_pg_uL)
  
  plot_predicted_theoretical_conc <- ggplot() +
    geom_point(data = analytes_concentrations %>%
                 drop_na(Theoretical_conc_uM) %>%
                 group_by(Compound, Theoretical_conc_uM) %>%
                 mutate(conc_pred_uM = mean(conc_pred_uM)) %>%
                 ungroup(),
               mapping = aes(Theoretical_conc_uM, conc_pred_uM,
                             color = Compound)) +
    geom_abline(slope = 1, intercept = 0) +
    scale_y_log10(limits = c(10^-5, 10^0)) +
    scale_x_log10(limits = c(10^-5, 10^0)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_abline(slope = 1, intercept = 1) +
    geom_abline(slope = 1, intercept = -1) +
    theme(aspect.ratio = 1) +
    theme_classic()
  ggplotly(plot_predicted_theoretical_conc)
  
  
  #Return list with data, plot
  data_predicted = list("all_calibration_plots" = plot_calgraph,
                        "plot_predictedIE_slope" = plot_predictedIE_slope,
                        "plot_predicted_theoretical_conc" = plot_predicted_theoretical_conc,
                        "data" = analytes_concentrations)
  
  return(data_predicted)
}





