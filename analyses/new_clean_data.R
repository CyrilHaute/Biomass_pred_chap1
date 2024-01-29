############# Load environmental covariates #############

environment_file <- list.files("new_data/new_raw_data/environmental_covariates", full.names = TRUE)

environment_cov_list <- lapply(1:length(environment_file), function(i) {
  
  cov_i <- environment_file[i]
  
  environment_file_i <- list.files(cov_i, full.names = TRUE)
  environment_file_i_small <- list.files(cov_i, full.names = FALSE)
  
  cov_i <- environment_file_i[which(grepl("final", environment_file_i) == TRUE)]
  
  load(cov_i)
  
  final_derived_data_joined <- final_derived_data_joined[,colnames(final_derived_data_joined)[!colnames(final_derived_data_joined) %in% c("x", "y", "name_latitude", "name_longitude", "name_var_date")]]
  
  assign(gsub(".Rdata", "", gsub("final_", "", environment_file_i_small[which(grepl("final", environment_file_i_small) == TRUE)])), final_derived_data_joined)
  
})

rls_env <- purrr::reduce(environment_cov_list, dplyr::inner_join)

rls_env <- rls_env |> 
  dplyr::rename(survey_id = name_data_id)

############# Load habitat covariates #############

load("new_data/new_raw_data/final_habitat.Rdata")

# Assess reef extent from reef geomorphology

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

############# Load social covariates #############

load("new_data/new_raw_data/social_covariates/corruption.RData")
load("new_data/new_raw_data/social_covariates/gdp.RData")
load("new_data/new_raw_data/social_covariates/gravity.RData")
load("new_data/new_raw_data/social_covariates/hdi.RData")
load("new_data/new_raw_data/social_covariates/marine_ecosystem_dependency.RData")
load("new_data/new_raw_data/social_covariates/mpa.Rdata")
load("new_data/new_raw_data/social_covariates/natural_ressource_rent.RData")
load("new_data/new_raw_data/social_covariates/neartt.RData")
load("new_data/new_raw_data/social_covariates/ngo.RData")

ngo <- ngo[which(is.na(ngo$ngo) == FALSE),]

rls_soc <- corruption |> 
  dplyr::inner_join(gdp) |> 
  dplyr::inner_join(gravity) |> 
  dplyr::inner_join(hdi) |> 
  dplyr::inner_join(med) |> 
  dplyr::inner_join(mpa) |> 
  dplyr::inner_join(natural_ressource_rent) |> 
  dplyr::inner_join(neartt) |> 
  dplyr::inner_join(ngo)

rls_soc$effectiveness <- as.factor(rls_soc$effectiveness)

rls_covariates <- rls_env |> 
  dplyr::inner_join(rls_soc) |> 
  dplyr::inner_join(habitat)

rls_covariates$survey_id <- as.character(rls_covariates$survey_id)

############################################################################################

load("new_data/new_raw_data/RLS_actinopterygii_data.Rdata")
load("new_data/new_raw_data/00_rls_surveys.Rdata")

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

# diff_survey <- habitat[!habitat$survey_id %in% rls_spread_coral_reef$survey_id,]

############# Variable selection #############

rls_env_selec <- rls_fish_cov[,colnames(rls_env)]
rls_env_selec <- rls_env_selec[,which(grepl(pattern = paste0(c("max", "min", "mean"), collapse = "|"), x = colnames(rls_env_selec)) == TRUE)]

cor_env <- stats::cor(rls_env_selec) #Look at correlation between covariates
corrplot::corrplot(cor_env, type = "upper") 

rls_env_selec2 <- rls_env_selec[,c("min_5year_o2", "mean_1year_chl", "mean_1year_so_mean", "mean_1year_ph" , "min_1year_analysed_sst", "min_1year_degree_heating_week", "max_5year_degree_heating_week")]
cor_env2 <- stats::cor(rls_env_selec2)
corrplot::corrplot(cor_env2, type = "upper")


rls_soc_selec <- rls_fish_cov[,colnames(rls_soc)]

cor_soc <- stats::cor(rls_soc_selec[,-c(1:3, 9)])

rls_soc_selec2 <- rls_soc_selec[,c("gravtot2", "neartt", "gdp", "hdi", "natural_ressource_rent")]
cor_soc2 <- stats::cor(rls_soc_selec2)
