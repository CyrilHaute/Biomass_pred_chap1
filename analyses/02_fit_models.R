# This script run the six biomass prediction models (glm, gam, rf, sprf, spamm and brt)

source("R/02_glm_function_SCV.R")
source('scripts-final/00_functions/model-functions/gam_function_SCV.R')
source('scripts-final/00_functions/model-functions/rf_function_SCV.R')
source('scripts-final/00_functions/model-functions/spatialrf_function_SCV.R')
source('scripts-final/00_functions/model-functions/spamm_function_SCV.R')
source('scripts-final/00_functions/model-functions/brt_function_SCV.R')


# load fish biomass data and covariates
load("new_data/new_derived_data/biomass_scv.RData")
load("new_data/new_derived_data/rls_covariates.RData")

# set up base-directory for file saves
base_dir <- 'outputs/biomass_prediction/'

# run glm 
print('glm biomass prediction')
glm_function(biomass = biomass_scv,
             covariates = rls_covariates,
             species_name = colnames(biomass_scv[[1]]$fitting)[!colnames(biomass_scv[[1]]$fitting) %in% c("survey_id", "latitude", "longitude")],
             base_dir = base_dir)

# run gam
print('gam biomass prediction')
gam_function(biomass = rls_biomass_SCV,
             covariates = covariates,
             species_name = colnames(rls_biomass_SCV[[1]]$fitting[,-1]),
             base_dir = base_dir)

# run random forest
print('rf biomass prediction')
rf_function(biomass = rls_biomass_SCV,
            covariates = covariates,
            species_name = colnames(rls_biomass_SCV[[1]]$fitting[,-1]),
            base_dir = base_dir)

# run spatial random forest
print('sprf biomass prediction')
spatialrf_function(biomass = rls_biomass_SCV,
                   covariates = spatial_covariates,
                   species_name = colnames(rls_biomass_SCV[[1]]$fitting[,-1]),
                   base_dir = base_dir)

# run spamm (GLMM)
print('spamm biomass prediction')
spamm_function(biomass = rls_biomass_SCV,
               covariates = spatial_covariates,
               species_name = colnames(rls_biomass_SCV[[1]]$fitting[,-1]),
               base_dir = base_dir)

# run boosted regression trees
print('brt biomass prediction')
brt_function(biomass = rls_biomass_SCV,
             covariates = covariates,
             species_name = colnames(rls_biomass_SCV[[1]]$fitting[,-1]),
             n.cores=1,
             base_dir = base_dir)
