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

metrics = c('Amae_rel_mean', 'Intercept', 'Slope', 'Pearson', 'Spearman', 'Pdispersion')
  
# select best fitted model for each model type based on a concensus metrics ----
  
# read and clean assessment data as in the script 01-model-performance-figures.R
  
all_assessments_SCV <- readRDS("results/model_assessment_all_R3/SCV/validation.rds")
all_assessments_SCV <- do.call(rbind, all_assessments_SCV)

# select only the columns to be used later 
all_assessments_SCV <- all_assessments_SCV %>% 
  
  dplyr::select(fitted_model, species_name, metrics, Evaluation_number, Evaluation_message)

# estimate for each species the best model based on performance metrics  
best_models <- all_assessments_SCV %>% 
  # estimate the relative metric performance within a cross validation and dataset
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(., 
                                                                  metrics = c('Amae_rel_mean', 'Intercept', 'Slope', 
                                                                              'Pearson', 'Spearman', 'Pdispersion')))) %>% 
  .$metric_aggregation %>% 
  do.call(rbind, .) %>% 
  # find the best fitting model for each species within each fitted_model
  group_by(species_name) %>% 
  do(best_model = .$fitted_model[which.max(.$discrimination)]) %>% 
  unnest(cols = c('best_model'))

saveRDS(best_models, file = 'results/overall_best_models_R3.rds')
  
predictions_all <- readRDS("results/predictions_all_R3/bind/SCV/rls_predictions_biomass.rds")
predictions_all <- mclapply(1:length(predictions_all), function(i) {
  test <- predictions_all[[i]]
  test <- test[,1:15]
})
predictions_all <- do.call(rbind, predictions_all)

# get subsets
model_predictions_all <- predictions_all

# rm(model_predictions)
        
# get seperate data for validations and verifications
validation_data <- model_predictions_all %>% 
  group_by(fitted_model, species_name) %>%
  do(validation_observed_mean = .$validation_observed_mean[[1]][.$validation_observed_mean[[1]] > 0], 
     validation_predict_mean = .$validation_predict_mean[[1]][.$validation_observed_mean[[1]] > 0]) %>% 
  ungroup()
        
# remove species with identical observations (cannot fit a model here)
if(sum(sapply(validation_data$validation_observed_mean, function(x) length(unique(x))==1)) != 0){
   validation_data <- validation_data[-which(sapply(validation_data$validation_observed_mean, function(x) length(unique(x))==1)),]}

# create plots 
observed_predicted_plot(input_data = validation_data, 
                        nbins = 10, 
                        levels = c('GLM', 'GAM', 'SPAMM', 'RF', 'GBM', 'SPRF'))
