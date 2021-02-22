library(car)
library(caret)
library(caTools)
library(dplyr)
library(enviPat)
library(gbm)
library(gsubfn)
library(janitor)
library(Metrics)
library(randomForest)
library(rcdk)
library(rcdklibs)
library(rJava)
library(rjson)
library(RRF)
library(plotly)
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

