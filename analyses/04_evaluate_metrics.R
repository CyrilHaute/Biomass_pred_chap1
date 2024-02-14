# extract evaluation 

# all functions for evaluating outputs
source("R/04_models_evaluation_functions.R")

output_files <- list.files("outputs/biomass_prediction", full.names = T)

load_outputs <- lapply(1:length(output_files), function(i) {
  
  load(output_files[i])
  assign(unique(extracted_predictions$fitted_model), extracted_predictions)
  
})

# perform validation (out the bag)

model_assessment <- lapply(1:length(load_outputs), function(i){
  
  model_i <- load_outputs[[i]]
  
  model_i |>
    dplyr::rowwise() |>  
    dplyr::do(metrics = biomass_assessment_metrics(predictions   = .$validation_predict, 
                                                   observations  = .$validation_observed,
                                                   scale = NULL)) |> 
    tidyr::unnest(metrics) |> 
    dplyr::bind_cols(model_i[,c("species_name", "fitted_model")])
  
})

model_assessment <- do.call(rbind, model_assessment)

model_assessment <- model_assessment |> 
  dplyr::group_split(species_name)

model_assessment <- lapply(1:length(model_assessment), function(i){
  
  model_assessment <- model_assessment[[i]]
  
  })

dir.create("outputs/model_assessment_validation", recursive = T)

save(model_assessment, file = "outputs/model_assessment_validation/validation.Rdata")
