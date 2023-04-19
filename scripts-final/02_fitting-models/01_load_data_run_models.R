# script to fit each of the functions to species matrix

# load in packages ---- 

libs <- c('tidyverse', 'parallel', 'itertools', 'matrixStats', 'MASS', 'statmod', 'buildmer')
lapply(libs, library, character.only = T, lib.loc = '/home/marbec/R/x86_64-pc-linux-gnu-library/4.1')

# check all packages are loaded
if(sum(libs %in% (.packages())) != length(libs)){
  stop('packages not loaded correctly')}
  
####################### Biomass predictions #######################

# load in data ----
    
# load in biomass data for biomass predictions in spatial cross validation
  
rls_biomass_SCV <- readRDS("data/Cyril_data/rls_biomass_SCV.rds")
rls_biomass_SCV <- mclapply(1:length(rls_biomass_SCV), function(i){
  cv <- rls_biomass_SCV[[i]]
  cv$fitting <- cv$fitting[,c(1:38)]
  cv$validation <- cv$validation[,c(1:38)]
  cv
},mc.cores = 1)
  
# load in covariates
  
rls_cov_hab <- readRDS("data/Cyril_data/RLS_hab.rds")
rls_cov_env <- readRDS("data/Cyril_data/RLS_env.rds")
rls_cov_soc <- readRDS("data/Cyril_data/RLS_soc.rds")
rls_cov_mpa <- readRDS("data/Cyril_data/RLS_mpa2.rds")
covariates <- rls_cov_env %>% 
  inner_join(rls_cov_soc, by = "SurveyID") %>%
  inner_join(rls_cov_hab, by = "SurveyID") %>% 
  inner_join(rls_cov_mpa, by = "SurveyID")
  
# add coordonates (ONLY for spatial model : spaMM and SpatialRF)

rls_sitesInfos <- readRDS("data/Cyril_data/RLS_sitesInfos.rds")
spatial_covariates <- inner_join(covariates, rls_sitesInfos, by = "SurveyID")
spatial_covariates <- spatial_covariates %>%
  dplyr::rename(X = "SiteLongitude", Y = "SiteLatitude") %>%
  dplyr::select(names(covariates), X, Y, Depth)
spatial_covariates <- spatial_covariates[,c(1,23,24,2:22)]

# set up base-directory for file saves

base_dir <- 'results/rls_basic_all_R2/SCV'

# run biomass models

source('scripts-final/02_fitting-models/02_fit-models/fit_glm_SCV.R')
source('scripts-final/02_fitting-models/02_fit-models/fit_gam_SCV.R')
source('scripts-final/02_fitting-models/02_fit-models/fit_rf_SCV.R')
source('scripts-final/02_fitting-models/02_fit-models/fit_sprf_SCV.R')
source('scripts-final/02_fitting-models/02_fit-models/fit_spamm_SCV.R')
source('scripts-final/02_fitting-models/02_fit-models/fit_brt_SCV.R')




####################### Variable contribution #######################

# load in biomass data for biomass variable contribution estimation

rls_biomass_cont <- readRDS("data/Cyril_data/Fish_RLS_cont_var.rds")
rls_biomass_cont <- rls_biomass_cont[,1:16]

# load in covariates

rls_cov_hab <- readRDS("../Biomass_prediction/data/Cyril_data/RLS_hab.rds")
rls_cov_env <- readRDS("../Biomass_prediction/data/Cyril_data/RLS_env.rds")
rls_cov_soc <- readRDS("../Biomass_prediction/data/Cyril_data/RLS_soc.rds")
rls_cov_mpa <- readRDS("../Biomass_prediction/data/Cyril_data/RLS_mpa2.rds")
covariates <- rls_cov_env %>% 
  inner_join(rls_cov_soc, by = "SurveyID") %>%
  inner_join(rls_cov_hab, by = "SurveyID") %>% 
  inner_join(rls_cov_mpa, by = "SurveyID")
covariates_cont <- covariates

# add coordonates (ONLY for spatial model : spaMM and SpatialRF)

RLS_sitesInfos <- readRDS("../Biomass_prediction/data/Cyril_data/RLS_sitesInfos.rds")
spatial_covariates_cont <- inner_join(covariates_cont, RLS_sitesInfos, by = "SurveyID")
spatial_covariates_cont <- spatial_covariates_cont %>% 
  dplyr::rename(X = "SiteLongitude", Y = "SiteLatitude") %>% 
  dplyr::select(names(covariates), X, Y)
spatial_covariates_cont <- spatial_covariates_cont[,c(1,23,24,2:22)]

# set up base-directory for file saves

base_dir_cont <- 'results/contributions'

# run contribution models

source('scripts-final/02_fitting-models/03_fit-contributions/fit_glm_cont_var.R')
source('scripts-final/02_fitting-models/03_fit-contributions/fit_gam_cont_var.R')
source('scripts-final/02_fitting-models/03_fit-contributions/fit_rf_cont_var.R')
source('scripts-final/02_fitting-models/03_fit-contributions/fit_sprf_cont_var.R')
source('scripts-final/02_fitting-models/03_fit-contributions/fit_spamm_cont_var.R')
source('scripts-final/02_fitting-models/03_fit-contributions/fit_brt_cont_var.R')
