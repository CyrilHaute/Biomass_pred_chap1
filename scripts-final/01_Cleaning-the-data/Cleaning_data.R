libs <- c('tidyverse', 'parallel', 'pbmcapply', 'corrplot')
lapply(libs, library, character.only = T, lib.loc = '/home/marbec/R/x86_64-pc-linux-gnu-library/4.1')

# check all packages are loaded
if(sum(libs %in% (.packages())) != length(libs)){
  stop('packages not loaded correctly')}

source("scripts-final/00_functions/CV_cleaning_data_function.R")

######Load raw covariates data

rls_habitat <- readRDS("data/Cyril_data/RLS_habitat.rds")
rls_env <- readRDS("data/Cyril_data/RLS_env_spatio_temporal.rds")
rls_socio <- readRDS("data/Cyril_data/rls_socio_withoutNA.rds")
mpa <- readRDS("data/Cyril_data/RLS_mpa.rds")
rls_coral_cover <- read.csv("data/Cyril_data/RLS_benthic_data_imputed.txt", header = TRUE, sep = "") #coral cover from RLS
rls_sitesInfos <- readRDS("data/Cyril_data/RLS_sitesInfos.rds")
test <- dplyr::inner_join(rls_socio, rls_sitesInfos)

######Select environmental covariates

#As correlation between covariates change with the sites you consider, keep only the final number of sites (covariates dataset with the fewer number of sites)
rls_env <- rls_env[rls_env$SurveyID %in% rls_socio$SurveyID,]
cor_env <- cor(rls_env[,-1]) #Look at correlation between covariates
corrplot(cor_env, type = "upper") 
rls_env <- rls_env[,c(1,2,14,20,38,50,66,67)] #Select covariates

summary(rls_env) #Any transformation ?

rls_env[,-1] <- scale(rls_env[,-1], center = TRUE, scale = TRUE) #Scale with a mean of 0 and a standard deviation of 1

######Select social covariates

cor_soc <- cor(rls_socio[,-1])
corrplot(cor_soc, type = "upper")
rls_socio <- rls_socio[,c(1,2:4,10,13,15)]

summary(rls_socio)

rls_socio$gdp <- log10(rls_socio$gdp + 1)
rls_socio$gravtot2 <- log10(rls_socio$gravtot2 + mean(rls_socio$gravtot2)) #log10(x+mean(x)) transform gravtot2

rls_socio[,-1] <- scale(rls_socio[,-1], center = TRUE, scale = TRUE)

######MPA data

mpa$Effectiveness <- as.character(mpa$Effectiveness)
mpa <- mpa %>% replace_na(list(Effectiveness = "out"))
mpa <- mpa[,c(1,3)]
mpa$Effectiveness <- as.factor(mpa$Effectiveness)

######Select habitat covariates

rls_coral_cover <- rls_coral_cover[,c(1,5)]
rls_habitat <- dplyr::inner_join(rls_habitat, rls_coral_cover, by = "SurveyID") #join RLS coral cover and Allen Coral Atlas habitat data

rls_habitat <- dplyr::inner_join(rls_habitat, rls_sitesInfos[,c(1,13)], by = "SurveyID") #add depth

names(rls_habitat) <- c("SurveyID", "sand_500m", "coral_algae_500m", "rock_500m", "rubble_500m", "seagrass_500m",
                        "microalgal_mats_500m", "coral_algae_10km", "rock_10km", "rubble_10km", "sand_10km", 
                        "seagrass_10km", "microalgal_mats_10km", "reef_slope_500m", "outer_reef.flat_500m", 
                        "back_reef.flat_500m", "deep_lagoon_500m", "plateau_500m", "inner_reef.flat_500m", "reef_crest_500m",
                        "shallow_lagoon_500m", "sheltered_reef.slope_500m", "terrestrial_reef.flat_500m", "patch_reefs_500m",
                        "back_reef.slope_10km", "inner_reef.flat_10km", "outer_reef.flat_10km", "reef_crest_10km", "reef_slope_10km",
                        "shallow_lagoon_10km", "deep_lagoon_10km", "plateau_10km", "sheltered_reef.slope_10km", "terrestrial_reef.flat_10km",
                        "patch_reefs_10km", "coral", "depth")

rls_habitat <- rls_habitat |>  
  dplyr::rowwise() |> 
  dplyr::mutate(Sum_geo_10 = sum(dplyr::c_across(colnames(rls_habitat)[25:35]))) #assess reef extent from reef geomorphology

cor_hab <- cor(rls_habitat[,-1])
corrplot(cor_hab, type = "upper")
rls_habitat <- rls_habitat[,c(1,36,2:5,38,37)]

rls_habitat[,-1] <- scale(rls_habitat[,-1], center = TRUE, scale = TRUE)

#Save selected and scaled covariates

saveRDS(rls_env, "data/Cyril_data/RLS_env.rds")
saveRDS(rls_socio, "data/Cyril_data/RLS_soc.rds")
saveRDS(rls_habitat, "data/Cyril_data/RLS_hab.rds")
saveRDS(mpa, "data/Cyril_data/RLS_mpa2.rds")

#Load fish data

RLS_species_list <- read.csv("data/Cyril_data/Species_list_RLS_Aug2021.csv", header = TRUE)
RLS_fishdata <- readRDS("data/Cyril_data/RLS_fishdata.rds")

RLS_data <- RLS_fishdata |> 
  dplyr::left_join(rls_sitesInfos[, 1:2]) |>  # add sites
  dplyr::left_join(RLS_species_list[, c(2:4, 6:9)]) |>  # add species infos
  dplyr::mutate(species = dplyr::if_else(Level == "species", valid_name_FishBase, TAXONOMIC_NAME), # new column with unique species names
                species = dplyr::recode(species, "Ophieleotris spp." = "Giuris spp.")) # rename "Ophieleotris spp." as "Giuris spp." 

# There are 14 repeated observations
RLS_data |> 
  dplyr::group_by(dplyr::across()) |> 
  dplyr::count() |>  
  dplyr::filter(n > 1)
# Mainly observations from a species of wrasse (Ophthalmolepis lineolata)
# Remove these
RLS_data <- unique(RLS_data)
# Remove fish identified at a higher taxonomic level than "genus" (i.e., family, class)
RLS_data <- RLS_data |>  dplyr::filter(Level != "higher")
Fish_RLS <- RLS_data

Fish_RLS <- Fish_RLS[-which(grepl('spp.', Fish_RLS$TAXONOMIC_NAME, fixed = T)),]
Fish_RLS <- Fish_RLS[-which(grepl('sp.', Fish_RLS$TAXONOMIC_NAME, fixed = T)),]

#Keep only needed column

Fish_RLS <- Fish_RLS[,-c(2,4,5,7:13)]
Fish_RLS <- Fish_RLS |>  dplyr::inner_join(rls_sitesInfos, by = "SurveyID")

rls_cov_hab <- readRDS("data/Cyril_data/RLS_hab.rds")
rls_cov_env <- readRDS("data/Cyril_data/RLS_env.rds")
rls_cov_soc <- readRDS("data/Cyril_data/RLS_soc.rds")
rls_cov_mpa <- readRDS("data/Cyril_data/RLS_mpa2.rds")
covariates <- rls_cov_env %>% 
  inner_join(rls_cov_soc, by = "SurveyID") %>%
  inner_join(rls_cov_hab, by = "SurveyID") %>% 
  inner_join(rls_cov_mpa, by = "SurveyID")

Fish_RLS$Biomass <- as.numeric(levels(Fish_RLS$Biomass)[Fish_RLS$Biomass])
Fish_RLS <- Fish_RLS[-which(is.na(Fish_RLS$Biomass)),]
Fish_RLS <- Fish_RLS %>% filter(SurveyID %in% unique(covariates$SurveyID))
Fish_RLS <- Fish_RLS %>% group_by(TAXONOMIC_NAME, SiteCode, SiteLatitude, SiteLongitude, SurveyID, SurveyDate, Depth) %>% summarise(Biomass = mean(Biomass))
Fish_RLS <- Fish_RLS %>% group_by(TAXONOMIC_NAME) %>% mutate(Taille = n())
Fish_RLS <- Fish_RLS %>% filter(Taille >= 50)
Fish_RLS <- Fish_RLS[,-c(3,4,6,7,9)]
Fish_RLS$TAXONOMIC_NAME <- as.factor(Fish_RLS$TAXONOMIC_NAME) #Save fish species name in a vector
nm_sp <- levels(Fish_RLS$TAXONOMIC_NAME)


######For relative variable contribution, create a species matrix biomass

Fish_RLS_contribution <- Fish_RLS %>% spread(TAXONOMIC_NAME, Biomass)

SurveyID <- Fish_RLS_contribution$SurveyID

Fish_RLS_contribution <- mclapply(2:ncol(Fish_RLS_contribution), function(j){ #replace NA by 0 = pseudo-absence
  
  species_j <- Fish_RLS_contribution[,j]
  
  species_j[is.na(species_j)] <- 0
  
  species_j
  
}, mc.cores = 1)

Fish_RLS_contribution <- do.call(cbind, Fish_RLS_contribution)
Fish_RLS_contribution <- cbind(SurveyID, Fish_RLS_contribution)
saveRDS(Fish_RLS_contribution, file = "data/Cyril_data/Fish_RLS_cont_var.rds")

##################################################################################################################
########SPATIAL CROSS VALIDATION

Fish_RLS <- inner_join(Fish_RLS, rls_sitesInfos)

Fish_RLS <- Fish_RLS[,-c(5:15)]

Fish_RLS <- Fish_RLS %>% spread(TAXONOMIC_NAME,Biomass)

species <- mclapply(3:ncol(Fish_RLS), function(i){
  
  species_i <- Fish_RLS[,i] #select the i species
  
  species_i[is.na(species_i)] <- 0 #replace NA by 0
  
  species_i
  
}, mc.cores = 1)

species <- do.call(cbind, species)
species <- cbind(Fish_RLS[,c(1:2)], species)

biomass_scv <- CV(species, 10)

names(biomass_scv) <- sapply(1:length(biomass_scv), function(i) { paste0("cv_", i)})

saveRDS(biomass_scv, file = "data/Cyril_data/rls_biomass_SCV.rds")
