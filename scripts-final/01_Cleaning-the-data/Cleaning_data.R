libs <- c('tidyverse', 'parallel', 'pbmcapply', 'corrplot')
lapply(libs, library, character.only = T, lib.loc = '/home/marbec/R/x86_64-pc-linux-gnu-library/4.2')

# check all packages are loaded
if(sum(libs %in% (.packages())) != length(libs)){
  stop('packages not loaded correctly')}

######Load raw covariates data

rls_habitat <- readRDS("data/Cyril_data/RLS_habitat.rds")
rls_env <- readRDS("data/Cyril_data/RLS_env_spatio_temporal.rds")
rls_socio <- readRDS("data/Cyril_data/rls_socio_withoutNA.rds")
mpa <- readRDS("data/Cyril_data/RLS_mpa.rds")
rls_coral_cover <- read.csv("data/Cyril_data/RLS_benthic_data_imputed.txt", header = TRUE, sep = "") #coral cover from RLS
rls_sitesInfos <- readRDS("data/Cyril_data/RLS_sitesInfos.rds")

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
rls_habitat <- inner_join(rls_habitat, rls_coral_cover, by = "SurveyID") #join RLS coral cover and Allen Coral Atlas habitat data

rls_habitat <- inner_join(rls_habitat, rls_sitesInfos[,c(1,13)], by = "SurveyID") #add depth

names(rls_habitat) <- c("SurveyID", "sand_500m", "coral_algae_500m", "rock_500m", "rubble_500m", "seagrass_500m",
                        "microalgal_mats_500m", "coral_algae_10km", "rock_10km", "rubble_10km", "sand_10km", 
                        "seagrass_10km", "microalgal_mats_10km", "reef_slope_500m", "outer_reef.flat_500m", 
                        "back_reef.flat_500m", "deep_lagoon_500m", "plateau_500m", "inner_reef.flat_500m", "reef_crest_500m",
                        "shallow_lagoon_500m", "sheltered_reef.slope_500m", "terrestrial_reef.flat_500m", "patch_reefs_500m",
                        "back_reef.slope_10km", "inner_reef.flat_10km", "outer_reef.flat_10km", "reef_crest_10km", "reef_slope_10km",
                        "shallow_lagoon_10km", "deep_lagoon_10km", "plateau_10km", "sheltered_reef.slope_10km", "terrestrial_reef.flat_10km",
                        "patch_reefs_10km", "coral", "depth")

rls_habitat <- rls_habitat %>% 
  mutate(Sum_geo_10 = rowSums(.[,c(25:35)])) #assess reef extent from reef geomorphology

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

Fish_RLS <- readRDS("data/Cyril_data/RLS_fishdata.rds")

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

#create a list of 10 Fish_RLS dataset, for now they are all the same
biomass_scv <- list(Fish_RLS, Fish_RLS, Fish_RLS, Fish_RLS, Fish_RLS, Fish_RLS, Fish_RLS, Fish_RLS, Fish_RLS, Fish_RLS)

biomass_scv <- mclapply(1:length(biomass_scv), function(i){ #for each dataset, split by species
  
  cv_i <- biomass_scv[[i]]
  
  split_species <- cv_i %>% group_split(TAXONOMIC_NAME)
  
}, mc.cores=1)

#for each dataset and species, split the data into two subset : fitting (80% of the dataset) and validation (20% of the dataset)
biomass_scv <- pbmclapply(1:length(biomass_scv), function(i){
  
  cv_i <- biomass_scv[[i]]
  
  Fish_RLS <- mclapply(1:length(cv_i), function(j){
    
    species_j <- cv_i[[j]]
    
    smp_size <- floor(0.8 * nrow(species_j))
    
    train_ind <- sample(seq_len(nrow(species_j)), size = smp_size)
    
    train_test <- list(species_j[train_ind,], species_j[-train_ind,])
    
    train_test[[1]] <- train_test[[1]] %>%            #delete from the train set transects that appears in the same site than the point in the validation set
      filter(!SiteCode %in% train_test[[2]]$SiteCode) 
    
    names(train_test) <- c("fitting", "validation")
    
    train_test$fitting$set <- rep("fitting", nrow(train_test$fitting))
    
    train_test$validation$set <- rep("validation", nrow(train_test$validation))
    
    train_test <- do.call(rbind, train_test)
    
    train_test
    
  }, mc.cores = detectCores() - 1)
   
}, mc.cores = 1)


biomass_scv <- pbmclapply(1:length(biomass_scv), function(i){
  
  cv_i <- biomass_scv[[i]]
  
  cv_i <- do.call(rbind, cv_i) #for each cross validation dataset, bind all the species together
  
  cv_i <- cv_i %>% 
    group_split(set) #then split it by set, either fitting or validation

  cv_i_matrix <- mclapply(1:length(cv_i), function(j){ #Create a species matrix (each column is a species) => create 0
  
    cv_j <- cv_i[[j]]
    
    cv_j <- cv_j[,-c(4,5)] #delete the set column
    
    cv_j <- cv_j %>% spread(TAXONOMIC_NAME,Biomass) #spread function to create your species matrix
    
  }, mc.cores = detectCores() - 1)
  
}, mc.cores = 1)


biomass_scv <- pbmclapply(1:length(biomass_scv), function(i){ #replace NA by 0
  
  cv_i <- biomass_scv[[i]] #select the i cross validation dataset (between 1 and 10)
  
  fit_val <- mclapply(1:length(cv_i), function(j){
    
    cv_j <- cv_i[[j]] #select the j set dataset (between fitting and validation)
    
    species <- mclapply(2:ncol(cv_j), function(k){
      
      species_k <- cv_j[,k] #select the k species
      
      species_k[is.na(species_k)] <- 0 #replace NA by 0
      
      species_k
      
    }, mc.cores = 1)
    
    species <- do.call(cbind, species)
    
    species <- cbind(cv_j[,1], species)
    
  }, mc.cores = 1)
  
  names(fit_val) <- c("fitting", "validation")
  
  fit_val
  
}, mc.cores = detectCores() - 1)

names(biomass_scv) <- c("cv1", "cv2", "cv3", "cv4", "cv5", "cv6", "cv7", "cv8", "cv9", "cv10")

saveRDS(biomass_scv, file = "data/Cyril_data/rls_biomass_SCV.rds")
