# extract evaluation 

############## Don't need to execute this, go directly to the bind results ##############

# Merge all files together by species

# Extract names all elements in folder
# path <- "results/predictions_all_species"
# 
# all_files <- list.files(path = path)
# 
# Predictions_biomass <- mclapply(1:length(all_files), function(i) {
#   all_files_full <- list()
#   all_files_full[i] <- paste0(path, '/', all_files[i])
#   readRDS(all_files_full[[i]])
# },mc.cores = 1)
# 
# dir.create("results/predictions_bind")
# 
# saveRDS(Predictions_biomass, file = "results/predictions_bind/rls_predictions_biomass.rds")
# rm(list = ls())
# gc()

# all functions for evaluating outputs
source("R/04_models_evaluation_functions.R")

######################################## For spatial cross validation
output_files <- list.files("outputs/biomass_prediction", full.names = T)

##### For bind_files

load_outputs <- lapply(1:length(output_files), function(i) {
  
  load(output_files[i])
  assign(unique(extracted_predictions$fitted_model), extracted_predictions)
  
})

# outputs_bind <- purrr::reduce(load_outputs, dplyr::full_join)

# perform validation (out the bag)

model_assessment <- pbmcapply::pbmclapply(1:length(load_outputs), function(i){
  
  model_i <- load_outputs[[i]]
  
  model_i |>
    dplyr::rowwise() |>  
    dplyr::do(metrics = biomass_assessment_metrics(predictions   = .$validation_predict, 
                                                   observations  = .$validation_observed,
                                                   scale = NULL)) |> 
    tidyr::unnest(metrics) |> 
    dplyr::bind_cols(model_i[,c("species_name", "fitted_model")])
  
}, mc.cores = parallel::detectCores() - 1)

model_assessment <- do.call(rbind, model_assessment)

model_assessment <- model_assessment |> 
  dplyr::group_split(species_name)

model_assessment <- parallel::mclapply(1:length(model_assessment), function(i){
  
  model_assessment <- model_assessment[[i]]
  
  }, mc.cores = 1)

dir.create("outputs/model_assessment_validation", recursive = T)

save(model_assessment, file = "outputs/model_assessment_validation/validation.Rdata")
