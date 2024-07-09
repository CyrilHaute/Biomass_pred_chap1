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

species_name <- colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")]

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "latitude", "longitude", "site_code", species_name[i])]
  
  spatialrf_function(biomass = rls_biomass_i,
                     covariates = rls_covariates,
                     base_dir = base_dir)
  
}, mc.cores = 1)


########### Species-specific biomass models per realm ###########

rls_biomass_realm <- rls_biomass |> 
  dplyr::inner_join(rls_surveys[, c("survey_id", "realm")]) |> 
  dplyr::group_split(realm)

# Remove ecoregion with less than 50 transects
nrow_realm <- sapply(1:length(rls_biomass_realm), function(i) {

  nrow(rls_biomass_realm[[i]]) < 100

})

rls_biomass_realm <- rls_biomass_realm[which(nrow_realm == FALSE)]

base_dir <- "outputs/biomass_prediction_ecoregion/"

# run glm per realm
print("glm per realm")

pbmcapply::pbmclapply(1:length(rls_biomass_realm), function(i) {
  
  realm <- rls_biomass_realm[[i]]
  
  species_name <- colnames(realm)[!colnames(realm) %in% c("survey_id", "latitude", "longitude", "site_code", "realm")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(realm[, species_name[j]])) < round(nrow(realm) / 10) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  realm <- realm[, c("survey_id", "latitude", "longitude", "site_code", "realm", new_species_name)]
  
  dir.create(base_dir)
  
  eco_base_dir <- paste0(base_dir, stringr::str_replace_all(unique(realm$realm), " ", "_"), "/")
  
  glm_function(biomass = realm,
               covariates = rls_covariates,
               species_name = new_species_name,
               base_dir = eco_base_dir)
  
})


# run random forest per realm
print("rf per realm")

pbmcapply::pbmclapply(1:length(rls_biomass_realm), function(i) {
  
  realm <- rls_biomass_realm[[i]]
  
  species_name <- colnames(realm)[!colnames(realm) %in% c("survey_id", "latitude", "longitude", "site_code", "realm")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(realm[, species_name[j]])) < round(nrow(realm) / 10) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  realm <- realm[, c("survey_id", "latitude", "longitude", "site_code", "realm", new_species_name)]
  
  dir.create(base_dir)
  
  eco_base_dir <- paste0(base_dir, stringr::str_replace_all(unique(realm$realm), " ", "_"), "/")
  
  rf_function(biomass = realm,
              covariates = rls_covariates,
              species_name = new_species_name,
              base_dir = eco_base_dir)
  
})


# run gbm per realm
print("gbm per realm")

pbmcapply::pbmclapply(1:length(rls_biomass_realm), function(i) {
  
  realm <- rls_biomass_realm[[i]]
  
  species_name <- colnames(realm)[!colnames(realm) %in% c("survey_id", "latitude", "longitude", "site_code", "realm")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(realm[, species_name[j]])) < round(nrow(realm) / 10) })

  new_species_name <- species_name[which(new_sp == FALSE)]
  
  realm <- realm[, c("survey_id", "latitude", "longitude", "site_code", "realm", new_species_name)]
  
  dir.create(base_dir)
  
  eco_base_dir <- paste0(base_dir, stringr::str_replace_all(unique(realm$realm), " ", "_"), "/")
  
  brt_function(biomass = realm,
               covariates = rls_covariates,
               species_name = new_species_name,
               base_dir = eco_base_dir)
  
})


# run gam per realm
print("gam per realm")

pbmcapply::pbmclapply(1:length(rls_biomass_realm), function(i) {
  
  realm <- rls_biomass_realm[[i]]
  
  species_name <- colnames(realm)[!colnames(realm) %in% c("survey_id", "latitude", "longitude", "site_code", "realm")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(realm[, species_name[j]])) < round(nrow(realm) / 10) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  realm <- realm[, c("survey_id", "latitude", "longitude", "site_code", "realm", new_species_name)]
  
  dir.create(base_dir)
  
  eco_base_dir <- paste0(base_dir, stringr::str_replace_all(unique(realm$realm), " ", "_"), "/")
  
  gam_function(biomass = realm,
               covariates = rls_covariates,
               species_name = new_species_name,
               base_dir = eco_base_dir)
  
})


# run spamm per realm
print("spamm per realm")

pbmcapply::pbmclapply(1:length(rls_biomass_realm), function(i) {
  
  realm <- rls_biomass_realm[[i]]
  
  species_name <- colnames(realm)[!colnames(realm) %in% c("survey_id", "latitude", "longitude", "site_code", "realm")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(realm[, species_name[j]])) < round(nrow(realm) / 10) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  realm <- realm[, c("survey_id", "latitude", "longitude", "site_code", "realm", new_species_name)]
  
  rls_covariates_eco <- rls_covariates |> 
    dplyr::filter(survey_id %in% realm$survey_id)
  
  dir.create(base_dir)
  
  eco_base_dir <- paste0(base_dir, stringr::str_replace_all(unique(realm$realm), " ", "_"), "/")
  
  spamm_function(biomass = realm,
                 covariates = rls_covariates_eco,
                 species_name = new_species_name,
                 base_dir = eco_base_dir)
  
})


# run sprf per realm
print("sprf per realm")

pbmcapply::pbmclapply(1:length(rls_biomass_realm), function(i) {
  
  realm <- rls_biomass_realm[[i]]
  
  species_name <- colnames(realm)[!colnames(realm) %in% c("survey_id", "latitude", "longitude", "site_code", "realm")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(realm[, species_name[j]])) < round(nrow(realm) / 10) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  realm <- realm[, c("survey_id", "latitude", "longitude", "site_code", "realm", new_species_name)]
  
  dir.create(base_dir)
  
  eco_base_dir <- paste0(base_dir, stringr::str_replace_all(unique(realm$realm), " ", "_"), "/")
  
  pbmcapply::pbmclapply(1:length(new_species_name), function(k) {
    
    rls_biomass_i <- rls_biomass[, c("survey_id", "latitude", "longitude", "site_code", new_species_name[k])]
    
    spatialrf_function(biomass = rls_biomass_i,
                       covariates = rls_covariates,
                       base_dir = eco_base_dir)
    
  }, mc.cores = 1)
  
})
