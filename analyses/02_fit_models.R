# This script run the six biomass prediction models (glm, gam, rf, sprf, spamm and brt)

# General functions to run models
source("R/01_cross_validation_function.R")
source("R/01_noise_function.R")

# Model functions
source("R/02_glm_function_SCV.R")
source("R/02_gam_function_SCV.R")
source("R/02_rf_function_SCV.R")
source("R/02_spatialrf_function_SCV.R")
source("R/02_spamm_function_SCV.R")
source("R/02_brt_function_SCV.R")

# load fish biomass data and covariates
load("data/new_derived_data/rls_biomass.RData")
load("data/new_derived_data/rls_covariates.RData")
load("data/new_raw_data/00_rls_surveys.Rdata")

rls_surveys$survey_id <- as.character(rls_surveys$survey_id)

species_name <- colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "site_code", "latitude", "longitude")]

# run glm 
print("glm biomass prediction")

base_dir <- "outputs/glm_prediction2/"

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "longitude", "latitude", "site_code", species_name[i])]
  
  glm_function(biomass = rls_biomass_i,
               covariates = rls_covariates,
               species_name = species_name[i],
               base_dir = base_dir)
  
}, mc.cores = parallel::detectCores() - 1)

# run random forest
print("rf biomass prediction")

base_dir <- "outputs/rf_prediction2/"

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "longitude", "latitude", "site_code", species_name[i])]
  
  rf_function(biomass = rls_biomass_i,
              covariates = rls_covariates,
              species_name = species_name[i],
              base_dir = base_dir)
  
}, mc.cores = parallel::detectCores() - 1)

# run boosted regression trees
print("brt biomass prediction")

base_dir <- "outputs/brt_prediction2/"

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "longitude", "latitude", "site_code", species_name[i])]
  
  brt_function(biomass = rls_biomass_i,
               covariates = rls_covariates,
               species_name = species_name[i],
               base_dir = base_dir)
  
}, mc.cores = 1)

# run gam
print("gam biomass prediction")

base_dir <- "outputs/gam_prediction2/"

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "longitude", "latitude", "site_code", species_name[i])]
  
  gam_function(biomass = rls_biomass,
               covariates = rls_covariates,
               species_name = species_name[i],
               base_dir = base_dir)
  
}, mc.cores = 1)

# run spamm (GLMM)
print("spamm biomass prediction")

base_dir <- "outputs/spamm_prediction2/"

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "longitude", "latitude", "site_code", species_name[i])]
  
  spamm_function(biomass = rls_biomass_i,
                 covariates = rls_covariates,
                 species_name = species_name[i],
                 base_dir = base_dir)
  
}, mc.cores = 1)

# run spatial random forest
print("sprf biomass prediction")

base_dir <- "outputs/sprf_prediction2/"

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "longitude", "latitude", "site_code", species_name[i])]
  
  spatialrf_function(biomass = rls_biomass_i,
                     covariates = rls_covariates,
                     species_name = species_name[i],
                     base_dir = base_dir)
  
}, mc.cores = 1)

