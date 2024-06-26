# This script run the six biomass contribution models (glm, gam, rf, sprf, spamm and brt)

# General functions to run models contribution
source("R/01_noise_function.R")

# Contribution functions
source("R/03_glm_cont_var.R")
source("R/03_rf_cont_var.R")
source("R/03_gam_cont_var.R")
source("R/03_sprf_cont_var.R")
source("R/03_spamm_cont_var.R")
source("R/03_brt_cont_var.R")

# load fish biomass data and covariates
load("data/new_derived_data/rls_biomass.RData")
load("data/new_derived_data/rls_covariates.RData")
load("data/new_raw_data/00_rls_surveys.Rdata")

base_dir_contribution <- "outputs/biomass_contribution/"

# run glm for covariates contribution
print("glm biomass contribution")
glm_function_cont(biomass = rls_biomass,
                  covariates = rls_covariates,
                  species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")],
                  base_dir_cont = base_dir_contribution)

# run random Forest for covariates contribution
print("rf biomass contribution")
rf_function_cont(biomass = rls_biomass,
                 covariates = rls_covariates,
                 species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")],
                 base_dir_cont = base_dir_contribution)

# run gam for covariates contribution
print("gam biomass contribution")
gam_function_cont(biomass = rls_biomass,
                  covariates = rls_covariates,
                  species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")],
                  base_dir_cont = base_dir_contribution)

# run spamm for covariates contribution
print("spamm biomass contribution")
spamm_function_cont(biomass = rls_biomass,
                    covariates = rls_covariates,
                    species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")],
                    base_dir_cont = base_dir_contribution)

# run gbm for covariates contribution
print("gbm biomass contribution")
brt_function_cont(biomass = rls_biomass,
                  covariates = rls_covariates,
                  species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")],
                  base_dir_cont = base_dir_contribution)

# run spatial Random Forest for covariates contribution
print("sprf biomass contribution")

species_name <- colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")]

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "latitude", "longitude", "site_code", species_name[i])]
  
  spatialrf_function_cont(biomass = rls_biomass_i,
                          covariates = rls_covariates,
                          base_dir_cont = base_dir_contribution)
  
}, mc.cores = 1)
