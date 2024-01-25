load("outputs/RLS_actinopterygii_data.Rdata")
load("data/derived_data/00_rls_surveys.Rdata")

rls_surveys$survey_id <- as.character(rls_surveys$survey_id)
rls_fish_data <- dplyr::inner_join(RLS_actinopterygii_data, rls_surveys)

length(unique(rls_fish_data$survey_id))

rls_fish_mean_biomass <- rls_fish_data |> 
  dplyr::group_by(survey_id, site_code, species_name, latitude, longitude, survey_date, depth) |> 
  dplyr::summarise(biomass = mean(biomass))

rls_fish_mean_biomass_count <- rls_fish_mean_biomass |> 
  dplyr::group_by(species_name) |> 
  dplyr::mutate(count = dplyr::n()) |> 
  dplyr::filter(count >= 50)

length(unique(rls_fish_mean_biomass_count$survey_id))

rls_spread <- rls_fish_mean_biomass_count |> 
  tidyr::spread(species_name, biomass, fill = 0)


################################ Coral Reef

rls_coral_fish <- rls_fish_data |> 
  dplyr::filter(latitude < 30 & latitude > -30)

rls_coral_fish_mean_biomass <- rls_coral_fish |> 
  dplyr::group_by(survey_id, site_code, species_name, latitude, longitude, survey_date, depth) |> 
  dplyr::summarise(biomass = mean(biomass))

rls_coral_fish_mean_biomass_count <- rls_coral_fish_mean_biomass |> 
  dplyr::group_by(species_name) |> 
  dplyr::mutate(count = dplyr::n()) |> 
  dplyr::filter(count >= 50)

length(unique(rls_coral_fish_mean_biomass_count$survey_id))

rls_spread_coral_reef <- rls_coral_fish_mean_biomass_count |> 
  tidyr::spread(species_name, biomass, fill = 0)


load("new_data/new_raw_data/environmental_covariates/CHL/final_chl.Rdata")
chl <- final_derived_data_joined[,colnames(final_derived_data_joined)[!colnames(final_derived_data_joined) %in% c("x", "y", "name_latitude", "name_longitude", "name_var_date")]]
load("new_data/new_raw_data/environmental_covariates/DHW/final_dhw.Rdata")
dhw <- final_derived_data_joined[,colnames(final_derived_data_joined)[!colnames(final_derived_data_joined) %in% c("x", "y", "name_latitude", "name_longitude", "name_var_date")]]
load("new_data/new_raw_data/environmental_covariates/NPP/final_npp.Rdata")
npp <- final_derived_data_joined[,colnames(final_derived_data_joined)[!colnames(final_derived_data_joined) %in% c("x", "y", "name_latitude", "name_longitude", "name_var_date")]]
load("new_data/new_raw_data/environmental_covariates/O2/final_o2.Rdata")
o2 <- final_derived_data_joined[,colnames(final_derived_data_joined)[!colnames(final_derived_data_joined) %in% c("x", "y", "name_latitude", "name_longitude", "name_var_date")]]
load("new_data/new_raw_data/environmental_covariates/pH/final_ph.RData")
pH <- final_derived_data_joined[,colnames(final_derived_data_joined)[!colnames(final_derived_data_joined) %in% c("x", "y", "name_latitude", "name_longitude", "name_var_date")]]
load("new_data/new_raw_data/environmental_covariates/SSS/final_sss.Rdata")
sss <- final_derived_data_joined[,colnames(final_derived_data_joined)[!colnames(final_derived_data_joined) %in% c("x", "y", "name_latitude", "name_longitude", "name_var_date")]]
load("new_data/new_raw_data/environmental_covariates/SST/final_sst.Rdata")
sst <- final_derived_data_joined[,colnames(final_derived_data_joined)[!colnames(final_derived_data_joined) %in% c("x", "y", "name_latitude", "name_longitude", "name_var_date")]]

load("new_data/new_raw_data/mpa.Rdata")

rls_env <- chl |> 
  dplyr::inner_join(dhw) |> 
  dplyr::inner_join(npp) |> 
  dplyr::inner_join(o2) |> 
  dplyr::inner_join(pH) |> 
  dplyr::inner_join(sss) |> 
  dplyr::inner_join(sst) |> 
  dplyr::rename(survey_id = name_data_id)

rls_env <- rls_env[,c("survey_id", "longitude", "latitude", colnames(rls_env)[which(grepl("mean_7days", colnames(rls_env)) == TRUE)], 
                      colnames(rls_env)[which(grepl("min_1year", colnames(rls_env)) == TRUE)], colnames(rls_env)[which(grepl("max_1year", colnames(rls_env)) == TRUE)],
                      "max_5year_degree_heating_week")]
rls_env <- rls_env[!colnames(rls_env) %in% "min_1year_degree_heating_week"]

load("new_data/new_raw_data/final_habitat.Rdata")

load("new_data/new_raw_data/00_rls_surveys.Rdata")

habitat <- habitat |> 
  dplyr::inner_join(rls_surveys)

rls_env <- rls_env |> 
  dplyr::filter(latitude < max(habitat$latitude) & latitude > min(habitat$latitude))

cor_env <- stats::cor(rls_env[,-c(1:3)]) #Look at correlation between covariates
corrplot::corrplot(cor_env, type = "upper") 
