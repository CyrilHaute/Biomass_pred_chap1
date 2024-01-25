# R version 4.1.1

# script to fit each of the functions to species matrix

# load in packages ---- 

libs <- c('tidyverse', 'parallel', 'itertools', 'matrixStats', 'MASS', 'statmod', 'buildmer')
lapply(libs, library, character.only = T, lib.loc = '/home/marbec/R/x86_64-pc-linux-gnu-library/4.1')

# check all packages are loaded
if(sum(libs %in% (.packages())) != length(libs)){
  stop('packages not loaded correctly')}

rls_biomass_SCV <- readRDS("data/Cyril_data/rls_biomass_SCV.rds")
rls_biomass_SCV <- rls_biomass_SCV[[1]]
rls_biomass_SCV <- do.call(rbind, rls_biomass_SCV)

# load in covariates

rls_cov_hab <- readRDS("data/Cyril_data/RLS_hab.rds")
rls_cov_env <- readRDS("data/Cyril_data/RLS_env.rds")
rls_cov_soc <- readRDS("data/Cyril_data/RLS_soc.rds")
rls_cov_mpa <- readRDS("data/Cyril_data/RLS_mpa2.rds")
covariates <- rls_cov_env %>% 
  inner_join(rls_cov_soc, by = "SurveyID") %>%
  inner_join(rls_cov_hab, by = "SurveyID") %>% 
  inner_join(rls_cov_mpa, by = "SurveyID")

all_assessments_SCV <- readRDS("results/model_assessment_validation/validation.rds")
all_assessments_SCV <- do.call(rbind, all_assessments_SCV)

metrics = c('Intercept', 'Slope', 'Pearson', 'Spearman')

# select only the columns to be used later 
all_assessments_SCV <- all_assessments_SCV %>% 
  
  dplyr::select(fitted_model, species_name, metrics)

all_assessments_SCV <- all_assessments_SCV[all_assessments_SCV$Pearson > 0.5 & all_assessments_SCV$Spearman > 0.5,]
bestsp <- unique(all_assessments_SCV$species_name)

data_test <- inner_join(rls_biomass_SCV[c("SurveyID", bestsp[15])], covariates)
names(data_test)[2] <- "Biomass"
data_test <- data_test[data_test$Biomass > 0,]
ggplot(data_test, aes(x = log10(Biomass+1), y = mean_chl_1year)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = mean_DHW_1year)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = mean_DHW_5year)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = mean_pH_1year)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = mean_sss_1year)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = min_sst_1year)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = max_sst_1year)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = conflicts)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = Naturalresourcesrents)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = gravtot2)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = gdp)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = coral)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = sand_500m)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = coral_algae_500m)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = rock_500m)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = rubble_500m)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = Sum_geo_10)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = depth)) +
  geom_point()
ggplot(data_test, aes(x = log10(Biomass+1), y = Effectiveness)) +
  geom_point()



