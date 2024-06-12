# This script run the six biomass prediction models (glm, gam, rf, sprf, spamm and brt)

source("R/01_cross_validation_function.R")
source("R/01_noise_function.R")
source("R/02_spatialml.R")
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

base_dir <- "outputs/biomass_prediction/"

# run glm 
print("glm biomass prediction")
glm_function(biomass = rls_biomass,
             covariates = rls_covariates,
             species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")],
             base_dir = base_dir)

# run random forest
print("rf biomass prediction")
rf_function(biomass = rls_biomass,
            covariates = rls_covariates,
            species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")],
            base_dir = base_dir)

# run boosted regression trees
print("brt biomass prediction")
brt_function(biomass = rls_biomass,
             covariates = rls_covariates,
             species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")],
             base_dir = base_dir)

# run gam
print("gam biomass prediction")
gam_function(biomass = rls_biomass,
             covariates = rls_covariates,
             species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")],
             base_dir = base_dir)

# run spamm (GLMM)
print("spamm biomass prediction")
spamm_function(biomass = rls_biomass,
               covariates = rls_covariates,
               species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")],
               base_dir = base_dir)

# run spatial random forest
print("sprf biomass prediction")
spatialrf_function(biomass = rls_biomass,
                   covariates = rls_covariates,
                   species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")],
                   base_dir = base_dir)
