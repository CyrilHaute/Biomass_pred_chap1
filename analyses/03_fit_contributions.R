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


########### Variable contribution per ecoregion ###########

rls_biomass_ecoregion <- rls_biomass |> 
  dplyr::inner_join(rls_surveys[, c("survey_id", "ecoregion")]) |> 
  dplyr::group_split(ecoregion)

# Remove ecoregion with less than 50 transects
nrow_ecoregion <- sapply(1:length(rls_biomass_ecoregion), function(i) {
  
  nrow(rls_biomass_ecoregion[[i]]) < 50
  
})

rls_biomass_ecoregion <- rls_biomass_ecoregion[which(nrow_ecoregion == FALSE)]

base_dir_contribution <- "outputs/biomass_contribution_ecoregion/"

# run glm for covariates contribution per ecoregion
print("glm biomass contribution per ecoregion")

pbmcapply::pbmclapply(1:length(rls_biomass_ecoregion), function(i) {
  
  ecoregion <- rls_biomass_ecoregion[[i]]
  
  species_name <- colnames(ecoregion)[!colnames(ecoregion) %in% c("survey_id", "latitude", "longitude", "site_code", "ecoregion")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(ecoregion[, species_name[j]])) < round(nrow(ecoregion) / 5) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  ecoregion <- ecoregion[, c("survey_id", "latitude", "longitude", "site_code", "ecoregion", new_species_name)]
  
  dir.create(base_dir_contribution)
  
  eco_base_dir_contribution <- paste0(base_dir_contribution, stringr::str_replace_all(unique(ecoregion$ecoregion), " ", "_"), "/")
  
  glm_function_cont(biomass = ecoregion,
                    covariates = rls_covariates,
                    species_name = new_species_name,
                    base_dir_cont = eco_base_dir_contribution)
  
})


# run random Forest for covariates contribution per ecoregion
print("rf biomass contribution per ecoregion")

pbmcapply::pbmclapply(1:length(rls_biomass_ecoregion), function(i) {
  
  ecoregion <- rls_biomass_ecoregion[[i]]
  
  species_name <- colnames(ecoregion)[!colnames(ecoregion) %in% c("survey_id", "latitude", "longitude", "site_code", "ecoregion")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(ecoregion[, species_name[j]])) < round(nrow(ecoregion) / 5) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  ecoregion <- ecoregion[, c("survey_id", "latitude", "longitude", "site_code", "ecoregion", new_species_name)]
  
  dir.create(base_dir_contribution)
  
  eco_base_dir_contribution <- paste0(base_dir_contribution, stringr::str_replace_all(unique(ecoregion$ecoregion), " ", "_"), "/")
  
  rf_function_cont(biomass = ecoregion,
                   covariates = rls_covariates,
                   species_name = new_species_name,
                   base_dir_cont = eco_base_dir_contribution)
  
})


# run gam for covariates contribution per ecoregion
print("gam biomass contribution per ecoregion")

pbmcapply::pbmclapply(1:length(rls_biomass_ecoregion), function(i) {
  
  ecoregion <- rls_biomass_ecoregion[[i]]
  
  species_name <- colnames(ecoregion)[!colnames(ecoregion) %in% c("survey_id", "latitude", "longitude", "site_code", "ecoregion")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(ecoregion[, species_name[j]])) < round(nrow(ecoregion) / 5) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  ecoregion <- ecoregion[, c("survey_id", "latitude", "longitude", "site_code", "ecoregion", new_species_name)]
  
  dir.create(base_dir_contribution)
  
  eco_base_dir_contribution <- paste0(base_dir_contribution, stringr::str_replace_all(unique(ecoregion$ecoregion), " ", "_"), "/")
  
  gam_function_cont(biomass = ecoregion,
                    covariates = rls_covariates,
                    species_name = new_species_name,
                    base_dir_cont = eco_base_dir_contribution)
  
})


# run spamm for covariates contribution per ecoregion
print("spamm biomass contribution per ecoregion")

pbmcapply::pbmclapply(1:length(rls_biomass_ecoregion), function(i) {
  
  ecoregion <- rls_biomass_ecoregion[[i]]
  
  species_name <- colnames(ecoregion)[!colnames(ecoregion) %in% c("survey_id", "latitude", "longitude", "site_code", "ecoregion")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(ecoregion[, species_name[j]])) < round(nrow(ecoregion) / 5) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  ecoregion <- ecoregion[, c("survey_id", "latitude", "longitude", "site_code", "ecoregion", new_species_name)]
  
  dir.create(base_dir_contribution)
  
  eco_base_dir_contribution <- paste0(base_dir_contribution, stringr::str_replace_all(unique(ecoregion$ecoregion), " ", "_"), "/")
  
  spamm_function_cont(biomass = ecoregion,
                      covariates = rls_covariates,
                      species_name = new_species_name,
                      base_dir_cont = eco_base_dir_contribution)
  
})


# run gbm for covariates contribution per ecoregion
print("gbm biomass contribution per ecoregion")

pbmcapply::pbmclapply(1:length(rls_biomass_ecoregion), function(i) {
  
  ecoregion <- rls_biomass_ecoregion[[i]]
  
  species_name <- colnames(ecoregion)[!colnames(ecoregion) %in% c("survey_id", "latitude", "longitude", "site_code", "ecoregion")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(ecoregion[, species_name[j]])) < round(nrow(ecoregion) / 5) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  ecoregion <- ecoregion[, c("survey_id", "latitude", "longitude", "site_code", "ecoregion", new_species_name)]
  
  dir.create(base_dir_contribution)
  
  eco_base_dir_contribution <- paste0(base_dir_contribution, stringr::str_replace_all(unique(ecoregion$ecoregion), " ", "_"), "/")

  brt_function_cont(biomass = ecoregion,
                    covariates = rls_covariates,
                    species_name = new_species_name,
                    base_dir_cont = eco_base_dir_contribution)
  
})


# run sprf for covariates contribution per ecoregion
print("sprf biomass contribution per ecoregion")

pbmcapply::pbmclapply(1:length(rls_biomass_ecoregion), function(i) {
  
  ecoregion <- rls_biomass_ecoregion[[i]]
  
  species_name <- colnames(ecoregion)[!colnames(ecoregion) %in% c("survey_id", "latitude", "longitude", "site_code", "ecoregion")]
  
  new_sp <- sapply(1:length(species_name), function(j) { nrow(unique(ecoregion[, species_name[j]])) < round(nrow(ecoregion) / 5) })
  
  new_species_name <- species_name[which(new_sp == FALSE)]
  
  ecoregion <- ecoregion[, c("survey_id", "latitude", "longitude", "site_code", "ecoregion", new_species_name)]
  
  dir.create(base_dir_contribution)
  
  eco_base_dir_contribution <- paste0(base_dir_contribution, stringr::str_replace_all(unique(ecoregion$ecoregion), " ", "_"), "/")
  
  pbmcapply::pbmclapply(1:length(new_species_name), function(i) {
    
    rls_biomass_i <- rls_biomass[, c("survey_id", "latitude", "longitude", "site_code", new_species_name[i])]
    
    spatialrf_function_cont(biomass = rls_biomass_i,
                            covariates = rls_covariates,
                            base_dir_cont = eco_base_dir_contribution)
    
  }, mc.cores = 1)
  
})


