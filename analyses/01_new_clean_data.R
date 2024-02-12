# This script load all the needed data, (covariates and fish biomass data), 
# select the covariates used in the models and create a spatial cross validation dataset

source("R/01_cross_validation_function.R")

#########################################################
############# Load environmental covariates #############
#########################################################

# List environmental covariates
environment_file <- list.files("data/new_raw_data/environmental_covariates", full.names = TRUE)

# Load environmental covariates within a list
environment_cov_list <- lapply(1:length(environment_file), function(i) {
  
  cov_i <- environment_file[i]
  
  environment_file_i <- list.files(cov_i, full.names = TRUE)
  environment_file_i_small <- list.files(cov_i, full.names = FALSE)
  
  cov_i <- environment_file_i[which(grepl("final", environment_file_i) == TRUE)]
  
  load(cov_i)
  
  final_derived_data_joined <- final_derived_data_joined[,colnames(final_derived_data_joined)[!colnames(final_derived_data_joined) %in% c("x", "y", "name_latitude", "name_longitude", "name_var_date")]]
  
  assign(gsub(".Rdata", "", gsub("final_", "", environment_file_i_small[which(grepl("final", environment_file_i_small) == TRUE)])), final_derived_data_joined)
  
})

# Bind all covariates in one dataframe
rls_env <- purrr::reduce(environment_cov_list, dplyr::inner_join)

rls_env <- rls_env |> 
  dplyr::rename(survey_id = name_data_id)

###################################################
############# Load habitat covariates #############
###################################################

load("data/new_raw_data/habitat_covariates/final_habitat.Rdata")
load("data/new_raw_data/habitat_covariates/Benthic_composition_inferred_tropical.Rdata")

# As the sum of habitats within a buffer do not always fit 100%, we converted percentage to relative in order to remove deep sea and land cover
# and acount only for coral habitats. First we assessed the sum of all habitat (either benthic and geomorphologic), for both buffers (500m, 10km),
# to keep the raw information of the sum of all habitat.
###################################################################################### changer chiffre par nom colonne
rls_habitat <- habitat |>  
  dplyr::rowwise() |> 
  dplyr::mutate(sum_geo_10 = sum(dplyr::c_across(colnames(habitat)[25:35])), # Assess reef extent from reef geomorphology
                sum_bent_500 = sum(dplyr::c_across(colnames(habitat)[2:7])), # Assess sum of habitat either benthic and geomorphomogic for both buffers 500m and 10km
                sum_geo_500 = sum(dplyr::c_across(colnames(habitat)[8:18])),
                sum_bent_10 = sum(dplyr::c_across(colnames(habitat)[19:24])))

# Convert to relatif percentage

rls_habitat <- rls_habitat |> 
  dplyr::rowwise() |> 
  dplyr::mutate((dplyr::across(colnames(habitat)[2:7]) * 100) / sum_bent_500,
                (dplyr::across(colnames(habitat)[8:18]) * 100) / sum_geo_500,
                (dplyr::across(colnames(habitat)[19:24]) * 100) / sum_bent_10,
                (dplyr::across(colnames(habitat)[25:35]) * 100) / sum_geo_10)

rls_habitat <- rls_habitat |> 
  dplyr::rename(reef_extent = sum_geo_10)

rls_habitat$survey_id <- as.character(rls_habitat$survey_id)
rls_habitat <- rls_habitat |> 
  dplyr::inner_join(inferred_benthos)

rls_habitat <- rls_habitat |> 
  dplyr::select(-c(longitude, latitude))

##################################################
############# Load social covariates #############
##################################################

load("data/new_raw_data/social_covariates/corruption.RData")
load("data/new_raw_data/social_covariates/gdp.RData")
load("data/new_raw_data/social_covariates/gravity.RData")
load("data/new_raw_data/social_covariates/hdi.RData")
load("data/new_raw_data/social_covariates/marine_ecosystem_dependency.RData")
load("data/new_raw_data/social_covariates/mpa.Rdata")
load("data/new_raw_data/social_covariates/natural_ressource_rent.RData")
load("data/new_raw_data/social_covariates/neartt.RData")
load("data/new_raw_data/social_covariates/ngo.RData")
load("data/new_raw_data/social_covariates/fishing_boat.RData")

ngo <- ngo[which(is.na(ngo$ngo) == FALSE),]

rls_soc <- corruption |> 
  dplyr::inner_join(gdp) |> 
  dplyr::inner_join(gravity) |> 
  dplyr::inner_join(hdi) |> 
  dplyr::inner_join(med) |> 
  dplyr::inner_join(mpa) |> 
  dplyr::inner_join(natural_ressource_rent) |> 
  dplyr::inner_join(neartt) |> 
  dplyr::inner_join(ngo) |> 
  dplyr::inner_join(fishing_boat)

rls_soc$effectiveness <- as.factor(rls_soc$effectiveness)

rls_habitat$survey_id <- as.integer(rls_habitat$survey_id)

rls_covariates <- rls_env |> 
  dplyr::inner_join(rls_soc) |> 
  dplyr::inner_join(rls_habitat)

rls_covariates$survey_id <- as.character(rls_covariates$survey_id)

##################################################
############# Load fish biomass data #############
##################################################

load("data/new_raw_data/RLS_actinopterygii_data.Rdata")
load("data/new_raw_data/00_rls_surveys.Rdata")

rls_surveys$survey_id <- as.character(rls_surveys$survey_id)
rls_fish_data <- dplyr::inner_join(RLS_actinopterygii_data, rls_surveys)

rls_coral_fish <- rls_fish_data |>
  dplyr::filter(survey_id %in% habitat$survey_id)

rls_coral_fish_mean_biomass <- rls_coral_fish |> 
  dplyr::group_by(survey_id, site_code, species_name, latitude, longitude, survey_date, depth) |> 
  dplyr::summarise(biomass = mean(biomass))

rls_coral_fish_mean_biomass_count <- rls_coral_fish_mean_biomass |> 
  dplyr::group_by(species_name, ) |> 
  dplyr::mutate(count = dplyr::n()) |> 
  dplyr::filter(count >= 50) |>
  dplyr::select(survey_id, species_name, biomass, latitude, longitude)

length(unique(rls_coral_fish_mean_biomass_count$survey_id))

rls_spread_coral_reef <- rls_coral_fish_mean_biomass_count |>
  tidyr::spread(species_name, biomass, fill = 0)

rls_fish_cov <- rls_spread_coral_reef |>
  dplyr::inner_join(rls_covariates)

##############################################
############# Variable selection #############
##############################################

# Select environmental covariates
rls_env_selec <- rls_fish_cov[,colnames(rls_env)]
rls_env_selec <- rls_env_selec[,which(grepl(pattern = paste0(c("max", "min", "mean"), collapse = "|"), x = colnames(rls_env_selec)) == TRUE)]

# Scale with a mean of 0 and a standard deviation of 1
rls_env_selec <- scale(rls_env_selec, center = TRUE, scale = TRUE)

cor_env <- stats::cor(rls_env_selec) #Look at correlation between covariates
corrplot::corrplot(cor_env, type = "upper") 

rls_env_selec2 <- rls_env_selec[,c("min_5year_ph", "mean_1year_chl", "mean_1year_so_mean" , "min_1year_analysed_sst", "max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_7days_chl", "mean_7days_analysed_sst")]
cor_env2 <- stats::cor(rls_env_selec2)
corrplot::corrplot(cor_env2, type = "upper")

rls_env_final <- rls_fish_cov[,c("survey_id", "min_5year_ph", "mean_1year_chl", "mean_1year_so_mean" , "min_1year_analysed_sst", "max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_7days_chl", "mean_7days_analysed_sst")]


# Select social covariates
rls_soc_selec <- rls_fish_cov[,colnames(rls_soc)]

# Transform gdp and gravity and scale with a mean of 0 and a standard deviation of 1
rls_soc_selec$gdp <- log10(rls_soc_selec$gdp + 1)
rls_soc_selec$gravtot2 <- log10(rls_soc_selec$gravtot2 + mean(rls_soc_selec$gravtot2))
rls_soc_selec[,-c(1:3, 9)] <- scale(rls_soc_selec[,-c(1:3, 9)], center = TRUE, scale = TRUE)

cor_soc <- stats::cor(rls_soc_selec[,-c(1:3, 9)])
corrplot::corrplot(cor_soc)

rls_soc_selec2 <- rls_soc_selec[,c("gravtot2", "neartt", "gdp", "hdi", "natural_ressource_rent", "n_fishing_vessels", "ngo")]

cor_soc2 <- stats::cor(rls_soc_selec2)
corrplot::corrplot(cor_soc2)

rls_soc_final <- rls_fish_cov[,c("survey_id", "gravtot2", "neartt", "gdp", "hdi", "natural_ressource_rent", "n_fishing_vessels", "ngo", "effectiveness")]


# Select habitat covariates
rls_hab_selec <- rls_fish_cov[, colnames(rls_habitat)]
rls_hab_selec <- rls_hab_selec |> 
  dplyr::inner_join(rls_surveys[,c("survey_id", "depth")])

cor_hab <- stats::cor(rls_hab_selec[,-c(1)])
corrplot::corrplot(cor_hab, type = "upper")

colnames(inferred_benthos)[-c(1:3)]
rls_hab_selec2 <- rls_hab_selec[,c("depth", "reef_extent", "coral_algae_500m", "Sand_500m", "Rock_500m", "Rubble_500m", "coral", "coralline algae")]
cor_hab2 <- stats::cor(rls_hab_selec2)
corrplot::corrplot(cor_hab2, type = "upper")

rls_hab_final <- rls_hab_selec[,c("survey_id", "depth", "reef_extent", "coral_algae_500m", "Sand_500m", "Rock_500m", "Rubble_500m", "coral", "coralline algae")]


#####################################################################
############# Create a spatial cross validation dataset #############
#####################################################################

species_name <- colnames(rls_spread_coral_reef)[!colnames(rls_spread_coral_reef) %in% c("survey_id", "latitude", "longitude")]

# species_name <- unique(rls_fish_cov$species_name)

rls_biomass <- rls_fish_cov |> 
  dplyr::inner_join(rls_surveys[,c("survey_id", "site_code", "depth")])

rls_biomass <- rls_biomass |>
  dplyr::select(survey_id, 
                latitude, 
                longitude, 
                site_code, 
                species_name, 
                colnames(rls_env_final)[!colnames(rls_env_final) %in% "survey_id"], 
                colnames(rls_soc_final)[!colnames(rls_soc_final) %in% "survey_id"], 
                colnames(rls_hab_final)[!colnames(rls_hab_final) %in% "survey_id"])

# biomass_scv <- pbmcapply::pbmclapply(1:length(species_name), function(i) {
#   
#   sp_i <- rls_biomass |> 
#     dplyr::select(survey_id, latitude, longitude, site_code, species_name[i])
#   
#   biomass_scv <- scv_function(dats = sp_i, 
#                               n.folds = 10)
#   
# }, mc.cores = parallel::detectCores() - 1)

biomass_scv <- scv_function(dats = rls_biomass,
                            n.folds = 10)

names(biomass_scv) <- sapply(1:length(biomass_scv), function(i) { paste0("cv_", i)})

# save derived data

rls_covariates <- rls_env_final |> 
  dplyr::inner_join(rls_soc_final) |> 
  dplyr::inner_join(rls_hab_final)

save(rls_covariates, file = "new_data/new_derived_data/rls_covariates.RData")
save(biomass_scv, file = "new_data/new_derived_data/biomass_scv.RData")
