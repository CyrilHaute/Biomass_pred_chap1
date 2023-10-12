# script to produce plots of predicted vs. modelled abundance across all model frameworks 
  
# load libraries ----
libs <- c('tidyverse', 'cowplot', 'RColorBrewer', 'gridExtra', 'data.table', 'pbapply', 'patchwork', 'parallel')
lapply(libs, library, character.only = T, lib.loc = '/home/marbec/R/x86_64-pc-linux-gnu-library/4.1')

# check all packages are loaded
if(sum(libs %in% (.packages())) != length(libs)){
  stop('packages not loaded correctly')}

# source functions ----
  
source('scripts-final/00_functions/models_evaluation_functions.R')
source('scripts-final/00_functions/model_performance_functions.R')
source('scripts-final/00_functions/model_prediction_functions.R')

predictions_all <- readRDS("results/predictions_bind/rls_predictions_biomass.rds")
predictions_all <- do.call(rbind, predictions_all)

# get subsets
model_predictions_all <- predictions_all

# rm(model_predictions)
        
# get seperate data for validations and verifications
validation_data <- model_predictions_all %>% 
  group_by(fitted_model, species_name) %>%
  do(validation_observed = .$validation_observed[[1]][.$validation_observed[[1]] > 0], 
     validation_predict = .$validation_predict[[1]][.$validation_observed[[1]] > 0]) %>% 
  ungroup()
        
# remove species with identical observations (cannot fit a model here)
if(sum(sapply(validation_data$validation_observed, function(x) length(unique(x))==1)) != 0){
   validation_data <- validation_data[-which(sapply(validation_data$validation_observed, function(x) length(unique(x))==1)),]}

# create plots 
observed_predicted_plot(input_data = validation_data, 
                        nbins = 10, 
                        levels = c('GLM', 'GAM', 'SPAMM', 'RF', 'GBM', 'SPRF'))
