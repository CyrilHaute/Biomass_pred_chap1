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

rls_surveys$survey_id <- as.character(rls_surveys$survey_id)

species_name <- colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")]

# run glm for covariates contribution
print("glm biomass contribution")

base_dir <- "outputs/glm_contribution/"

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "longitude", "latitude", "site_code", species_name[i])]
  
  glm_function_cont(biomass = rls_biomass_i,
                    covariates = rls_covariates,
                    species_name = species_name[i],
                    base_dir_cont = base_dir)
  
}, mc.cores = 1)

# run random Forest for covariates contribution
print("rf biomass contribution")

base_dir <- "outputs/rf_contribution/"

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "longitude", "latitude", "site_code", species_name[i])]
  
  rf_function_cont(biomass = rls_biomass_i,
                   covariates = rls_covariates,
                   species_name = species_name[i],
                   base_dir_cont = base_dir)
  
}, mc.cores = 1)

# run gbm for covariates contribution
print("gbm biomass contribution")

base_dir <- "outputs/brt_contribution/"

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "longitude", "latitude", "site_code", species_name[i])]
  
  brt_function_cont(biomass = rls_biomass_i,
                    covariates = rls_covariates,
                    species_name = species_name[i],
                    base_dir_cont = base_dir)
  
}, mc.cores = 1)

# run gam for covariates contribution
print("gam biomass contribution")

base_dir <- "outputs/gam_contribution/"

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "longitude", "latitude", "site_code", species_name[i])]
  
  gam_function_cont(biomass = rls_biomass_i,
                    covariates = rls_covariates,
                    species_name = species_name[i],
                    base_dir_cont = base_dir)
  
}, mc.cores = 1)

# run spamm for covariates contribution
print("spamm biomass contribution")

base_dir <- "outputs/spamm_contribution/"

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "longitude", "latitude", "site_code", species_name[i])]
  
  spamm_function_cont(biomass = rls_biomass_i,
                      covariates = rls_covariates,
                      species_name = species_name[i],
                      base_dir_cont = base_dir)
  
}, mc.cores = 1)

# run spatial Random Forest for covariates contribution
print("sprf biomass contribution")

base_dir <- "outputs/sprf_contribution/"

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  rls_biomass_i <- rls_biomass[, c("survey_id", "latitude", "longitude", "site_code", species_name[i])]
  
  spatialrf_function_cont(biomass = rls_biomass_i,
                          covariates = rls_covariates,
                          species_name = species_name[i],
                          base_dir_cont = base_dir)
  
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

base_dir_contribution <- "outputs/biomass_contribution_realm/"

# run glm per realm
print("glm per realm")

pbmcapply::pbmclapply(1:length(rls_biomass_realm), function(i) {
  
  realm <- rls_biomass_realm[[i]]
  
  species_name <- colnames(realm)[!colnames(realm) %in% c("survey_id", "latitude", "longitude", "site_code", "realm")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(realm[, species_name[j]])) < round(nrow(realm) / 10) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  realm <- realm[, c("survey_id", "latitude", "longitude", "site_code", "realm", new_species_name)]
  
  dir.create(base_dir_contribution)
  
  eco_base_dir <- paste0(base_dir_contribution, stringr::str_replace_all(unique(realm$realm), " ", "_"), "/")
  
  glm_function_cont(biomass = realm,
                    covariates = rls_covariates,
                    species_name = new_species_name,
                    base_dir_cont = eco_base_dir)
  
})


# run random forest per realm
print("rf per realm")

pbmcapply::pbmclapply(1:length(rls_biomass_realm), function(i) {
  
  realm <- rls_biomass_realm[[i]]
  
  species_name <- colnames(realm)[!colnames(realm) %in% c("survey_id", "latitude", "longitude", "site_code", "realm")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(realm[, species_name[j]])) < round(nrow(realm) / 10) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  realm <- realm[, c("survey_id", "latitude", "longitude", "site_code", "realm", new_species_name)]
  
  dir.create(base_dir_contribution)
  
  eco_base_dir <- paste0(base_dir_contribution, stringr::str_replace_all(unique(realm$realm), " ", "_"), "/")
  
  rf_function_cont(biomass = realm,
                   covariates = rls_covariates,
                   species_name = new_species_name,
                   base_dir_cont = eco_base_dir)
  
})


# run gbm per realm
print("gbm per realm")

pbmcapply::pbmclapply(1:length(rls_biomass_realm), function(i) {
  
  realm <- rls_biomass_realm[[i]]
  
  species_name <- colnames(realm)[!colnames(realm) %in% c("survey_id", "latitude", "longitude", "site_code", "realm")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(realm[, species_name[j]])) < round(nrow(realm) / 10) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  realm <- realm[, c("survey_id", "latitude", "longitude", "site_code", "realm", new_species_name)]
  
  dir.create(base_dir_contribution)
  
  eco_base_dir <- paste0(base_dir_contribution, stringr::str_replace_all(unique(realm$realm), " ", "_"), "/")
  
  brt_function_cont(biomass = realm,
                    covariates = rls_covariates,
                    species_name = new_species_name,
                    base_dir_cont = eco_base_dir)
  
})


# run gam per realm
print("gam per realm")

pbmcapply::pbmclapply(1:length(rls_biomass_realm), function(i) {
  
  realm <- rls_biomass_realm[[i]]
  
  species_name <- colnames(realm)[!colnames(realm) %in% c("survey_id", "latitude", "longitude", "site_code", "realm")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(realm[, species_name[j]])) < round(nrow(realm) / 10) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  realm <- realm[, c("survey_id", "latitude", "longitude", "site_code", "realm", new_species_name)]
  
  dir.create(base_dir_contribution)
  
  eco_base_dir <- paste0(base_dir_contribution, stringr::str_replace_all(unique(realm$realm), " ", "_"), "/")
  
  gam_function_cont(biomass = realm,
                    covariates = rls_covariates,
                    species_name = new_species_name,
                    base_dir_cont = eco_base_dir)
  
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
  
  dir.create(base_dir_contribution)
  
  eco_base_dir <- paste0(base_dir_contribution, stringr::str_replace_all(unique(realm$realm), " ", "_"), "/")
  
  spamm_function_cont(biomass = realm,
                      covariates = rls_covariates_eco,
                      species_name = new_species_name,
                      base_dir = eco_base_dir)
  
})


# run sprf per realm
print("sprf per realm")

for(i in 1:length(rls_biomass_realm)) {
  
  realm <- rls_biomass_realm[[i]]
  
  species_name <- colnames(realm)[!colnames(realm) %in% c("survey_id", "latitude", "longitude", "site_code", "realm")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(realm[, species_name[j]])) < round(nrow(realm) / 10) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  realm <- realm[, c("survey_id", "latitude", "longitude", "site_code", "realm", new_species_name)]
  
  dir.create(base_dir_contribution)
  
  eco_base_dir <- paste0(base_dir_contribution, stringr::str_replace_all(unique(realm$realm), " ", "_"), "/")
  
  pbmcapply::pbmclapply(1:length(new_species_name), function(k) {
    
    rls_biomass_i <- rls_biomass[, c("survey_id", "latitude", "longitude", "site_code", new_species_name[k])]
    
    spatialrf_function_cont(biomass = rls_biomass_i,
                            covariates = rls_covariates,
                            base_dir_cont = eco_base_dir)
    
  }, mc.cores = 1)
  
}
