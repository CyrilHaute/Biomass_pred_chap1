# This script run the six biomass prediction models (glm, gam, rf, sprf, spamm and brt)

# General functions to run models
source("R/01_noise_function.R")

# Model functions
source("R/09_partial_function.R")

# load fish biomass data and covariates
load("data/new_derived_data/rls_biomass.RData")
load("data/new_derived_data/rls_covariates.RData")
load("data/new_raw_data/00_rls_surveys.Rdata")
load("outputs/best_models.Rdata")

rls_surveys$survey_id <- as.character(rls_surveys$survey_id)

species_name <- best_models[best_models$best_model == "rf",]$species_name

########### Species-specific biomass models per realm ###########

rls_biomass_realm <- rls_biomass |> 
  dplyr::inner_join(rls_surveys[, c("survey_id", "realm")]) |> 
  dplyr::group_split(realm)

# Remove ecoregion with less than 50 transects
nrow_realm <- sapply(1:length(rls_biomass_realm), function(i) {
  
  nrow(rls_biomass_realm[[i]]) < 100
  
})

rls_biomass_realm <- rls_biomass_realm[which(nrow_realm == FALSE)]

base_dir <- "outputs/realm_partial/"

for(i in 1:length(rls_biomass_realm)) {

  realm <- rls_biomass_realm[[i]]
  
  species_name_realm <- colnames(realm)[!colnames(realm) %in% c("survey_id", "latitude", "longitude", "site_code", "realm")]

  new_sp <- sapply(1:length(species_name_realm), function(j) { nrow(unique(realm[, species_name_realm[j]])) < round(nrow(realm) / 10) })
  
  new_species_name <- species_name_realm[which(new_sp == FALSE)]
  new_species_name <- new_species_name[which(new_species_name %in% species_name)]
  
  realm <- realm[, c("survey_id", "latitude", "longitude", "site_code", "realm", new_species_name)]
  
  dir.create(base_dir)
  
  eco_base_dir <- paste0(base_dir, stringr::str_replace_all(unique(realm$realm), " ", "_"), "/")
  
  pbmcapply::pbmclapply(1:length(new_species_name), function(j) {
    
    realm_j <- realm[, c("survey_id", "longitude", "latitude", "site_code", new_species_name[j])]
    
    partial_function(biomass = realm_j,
                     covariates = rls_covariates,
                     species_name = new_species_name[j],
                     base_dir = eco_base_dir)
    
  }, mc.cores = parallel::detectCores() - 1)
  
}
