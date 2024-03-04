# source functions ----
  
source("R/04_model_prediction_functions.R")

output_files <- list.files("outputs/biomass_prediction", full.names = T)

load_outputs <- lapply(1:length(output_files), function(i) {
  
  load(output_files[i])
  assign(unique(extracted_predictions$fitted_model), extracted_predictions)
  
})
predictions_all <- do.call(rbind, load_outputs)

# get subsets
model_predictions_all <- predictions_all

# rm(model_predictions)
        
# get seperate data for validations and verifications
validation_data <- model_predictions_all |> 
  dplyr::group_by(fitted_model, species_name) |> 
  dplyr::do(validation_observed = .$validation_observed[[1]][.$validation_observed[[1]] > 0], 
            validation_predict = .$validation_predict[[1]][.$validation_observed[[1]] > 0]) |> 
  dplyr::ungroup()
        
# remove species with identical observations (cannot fit a model here)
if(sum(sapply(validation_data$validation_observed, function(x) length(unique(x))==1)) != 0){
   validation_data <- validation_data[-which(sapply(validation_data$validation_observed, function(x) length(unique(x))==1)),]}

# create plots 
observed_predicted_plot(input_data = validation_data, 
                        nbins = 10, 
                        levels = c('GLM', 'GAM', 'SPAMM', 'RF', 'GBM', 'SPRF'))
