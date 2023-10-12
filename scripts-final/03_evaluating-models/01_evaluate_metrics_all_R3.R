# extract evaluation 

library(tidyverse)
library(parallel)
library(pbmcapply)

# Merge all files together by species

# Extract names all elements in folder
path <- "results/predictions_all_species"

all_files <- list.files(path = path)

Predictions_biomass <- mclapply(1:length(all_files), function(i) {
  all_files_full <- list()
  all_files_full[i] <- paste0(path, '/', all_files[i])
  readRDS(all_files_full[[i]])
},mc.cores = 1)

dir.create("results/predictions_bind")

saveRDS(Predictions_biomass, file = "results/predictions_bind/rls_predictions_biomass.rds")
rm(list = ls())
gc()

# all functions for evaluating outputs
source('scripts-final/00_functions/models_evaluation_functions.R')

######################################## For spatial cross validation
bind_files <- list.files('results/predictions_bind', full.names = T)

##### For bind_files

clean_files <- readRDS(bind_files)

# perform validation (out the bag)

model_assessment <- pbmclapply(1:length(clean_files), function(i){
  clean_files <- clean_files[[i]]
  clean_files %>% 
    rowwise() %>% 
    do(metrics = biomass_assessment_metrics(predictions   = .$validation_predict, 
                                              observations  = .$validation_observed,
                                              scale = NULL)) %>% 
    unnest(metrics) %>% 
    bind_cols(clean_files[,1:5], .)
}, mc.cores = detectCores() -1)

model_assessment <- do.call(bind_rows, model_assessment)
model_assessment <- model_assessment %>% group_split(species_name)
model_assessment <- mclapply(1:length(model_assessment), function(i){
  model_assessment <- model_assessment[[i]]}, mc.cores=1)

dir.create('results/model_assessment_validation', recursive = T)

saveRDS(model_assessment, file = "results/model_assessment_validation/validation.rds")
