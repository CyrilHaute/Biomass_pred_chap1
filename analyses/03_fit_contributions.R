# This script run the six biomass contribution models (glm, gam, rf, sprf, spamm and brt)

source("R/03_glm_cont_var.R")


# load fish biomass data and covariates
load("data/new_derived_data/biomass_contribution.RData")
load("data/new_derived_data/rls_covariates.RData")

base_dir <- "outputs/biomass_contribution/"

# run glm 
print("glm biomass prediction")
glm_function(biomass = biomass_scv,
             covariates = rls_covariates,
             species_name = colnames(biomass_scv[[1]]$fitting)[!colnames(biomass_scv[[1]]$fitting) %in% c("survey_id", "latitude", "longitude")],
             base_dir = base_dir)