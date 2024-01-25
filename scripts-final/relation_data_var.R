# script to produce figures that evaluate model performances

# load in packages ---- 

libs <- c('tidyverse', 'gridExtra', 'ggplot2', 'patchwork', 'matrixStats', 'parallel','PNWColors')
lapply(libs, library, character.only = T, lib.loc = '/home/marbec/R/x86_64-pc-linux-gnu-library/4.1')

# check all packages are loaded
if(sum(libs %in% (.packages())) != length(libs)){
  stop('packages not loaded correctly')}

# Set palette colors for performance figures

pal_best = pnw_palette("Bay", 6 , type = "continuous")
pal_perf = pnw_palette("Bay",6, type = "continuous")

# select best fitted model for each model type based on a concensus metrics ----

metrics = c('Intercept', 'Slope', 'Pearson', 'Spearman')

# read data

all_assessments_SCV <- readRDS("results/model_assessment_validation/validation.rds")
all_assessments_SCV <- do.call(rbind, all_assessments_SCV)
all_assessments_SCV <- all_assessments_SCV[,c(1:2,12:15)]

sp_car <- readRDS("data/Cyril_data/RLS_species_traits.rds")
sp_car <- sp_car[,c(1,2,4,6,11)]

all_assessments_SCV <- inner_join(all_assessments_SCV, sp_car, by = "species_name")

ggplot(all_assessments_SCV, aes(x = log10(MaxLength), y = Pearson)) +
  facet_wrap(~fitted_model) +
  geom_point()
ggplot(all_assessments_SCV, aes(x = log10(MaxLength), y = Spearman)) +
  facet_wrap(~fitted_model) +
  geom_point()
ggplot(all_assessments_SCV, aes(x = Water.column, y = Pearson)) +
  facet_wrap(~fitted_model) +
  geom_boxplot()
ggplot(all_assessments_SCV, aes(x = Water.column, y = Spearman)) +
  facet_wrap(~fitted_model) +
  geom_boxplot()
ggplot(all_assessments_SCV, aes(x = Habitat, y = Pearson)) +
  facet_wrap(~fitted_model) +
  geom_boxplot()
ggplot(all_assessments_SCV, aes(x = Habitat, y = Spearman)) +
  facet_wrap(~fitted_model) +
  geom_boxplot()
all_assessments_SCV %>% group_by(Habitat) %>% summarise(n = n())
ggplot(all_assessments_SCV, aes(x = Trophic_guild_name, y = Pearson)) +
  facet_wrap(~fitted_model) +
  geom_boxplot()
ggplot(all_assessments_SCV, aes(x = Trophic_guild_name, y = Spearman)) +
  facet_wrap(~fitted_model) +
  geom_boxplot()

best_models <- readRDS("results/overall_best_models.rds")

best_assessments_SCV <- inner_join(all_assessments_SCV, best_models, by = "species_name")
best_assessments_SCV <- best_assessments_SCV[best_assessments_SCV$fitted_model == best_assessments_SCV$best_model,]


ggplot(best_assessments_SCV, aes(x = log10(MaxLength), y = Pearson)) +
  facet_wrap(~fitted_model) +
  geom_point()
ggplot(best_assessments_SCV, aes(x = log10(MaxLength), y = Spearman)) +
  facet_wrap(~fitted_model) +
  geom_point()
ggplot(best_assessments_SCV, aes(x = Water.column, y = Pearson)) +
  facet_wrap(~fitted_model) +
  geom_boxplot()
ggplot(best_assessments_SCV, aes(x = Water.column, y = Spearman)) +
  facet_wrap(~fitted_model) +
  geom_boxplot()
ggplot(best_assessments_SCV, aes(x = Habitat, y = Pearson)) +
  facet_wrap(~fitted_model) +
  geom_boxplot()
ggplot(best_assessments_SCV, aes(x = Habitat, y = Spearman)) +
  facet_wrap(~fitted_model) +
  geom_boxplot()
ggplot(best_assessments_SCV, aes(x = Trophic_guild_name, y = Pearson)) +
  facet_wrap(~fitted_model) +
  geom_boxplot()
ggplot(best_assessments_SCV, aes(x = Trophic_guild_name, y = Spearman)) +
  facet_wrap(~fitted_model) +
  geom_boxplot()

######################## tester relation entre occurence/biomasse et pearson/spearman #######################

rls_biomass_SCV <- readRDS("data/Cyril_data/rls_biomass_SCV.rds")
rls_biomass_SCV <- rls_biomass_SCV[[1]]
rls_biomass_SCV <- do.call(rbind, rls_biomass_SCV)
test <- mclapply(2:ncol(rls_biomass_SCV), function(i) {
  
  sp_i <- rls_biomass_SCV[,c(1,i)]
  names(sp_i)[2] <- "Biomass"
  df_sp_i <- tibble(species_name = names(rls_biomass_SCV)[i],
                    Occurence = nrow(sp_i[sp_i$Biomass > 0,]),
                    Biomass = sum(sp_i$Biomass))
}, mc.cores = 1)
test <- do.call(rbind, test)
all_assessments_SCV <- inner_join(all_assessments_SCV, test, by = "species_name")

ggplot(all_assessments_SCV, aes(x = Occurence, y = Pearson)) +
  facet_wrap(~fitted_model) +
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(all_assessments_SCV, aes(x = Occurence, y = Spearman)) +
  facet_wrap(~fitted_model) +
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(all_assessments_SCV, aes(x = Biomass, y = Pearson)) +
  facet_wrap(~fitted_model) +
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(all_assessments_SCV, aes(x = Biomass, y = Spearman)) +
  facet_wrap(~fitted_model) +
  geom_point() + 
  geom_smooth(method = "lm")


####################################### tester relation entre importance variables et traits spécifique #####################################

bind_files <- list.files('results/contributions_bind/', full.names = T)

##### For bind_files

plot_data <- readRDS(bind_files)

plot_level <- unique(plot_data$fitted_model)

ext_cont_biomass <- mclapply(1:length(plot_level), function(i) {
  
  only_model <- plot_data %>% filter(fitted_model == plot_level[i])
  only_model <- only_model %>% filter(species_name %in% best_models[best_models$best_model == plot_level[i],1]$species_name)
  
  ENV <- lapply(1:nrow(only_model), function(i) { only_model$contributions[[i]][c(11:17)]})
  ENV <- do.call(rbind, ENV)
  colnames(ENV) <- c("max_sst_1year", "mean_DHW_1year", "mean_DHW_5year", "mean_chl_1year", "mean_pH_1year", "mean_sss_1year", "min_sst_1year")
  ENV <- as.data.frame(ENV)
  ENV$species_name <- only_model$species_name
  ENV$cat <- rep("ENV", nrow(ENV))
  ENV <- ENV %>% inner_join(sp_car, by = "species_name")
  
  SOC <- lapply(1:nrow(only_model), function(i) { only_model$contributions[[i]][c(1:3,5,9,10,18)]})
  SOC <- do.call(rbind, SOC)
  colnames(SOC) <- c("MPA Effectiveness",  "Human Development Index", "Natural Ressource", "Conflicts", "GDP", "Gravity", "Neartt")
  SOC <- as.data.frame(SOC)
  SOC$species_name <- only_model$species_name
  SOC$cat <- rep("SOC", nrow(SOC))
  SOC <- SOC %>% inner_join(sp_car, "species_name")
  
  HAB <- lapply(1:nrow(only_model), function(i) { only_model$contributions[[i]][c(4,6:8,19:21)]})
  HAB <- do.call(rbind, HAB)
  colnames(HAB) <- c("Reef extent(Allen)" ,"Coral", "Coral_Algea(Allen)", "Depth", "Rock(Allen)", "Rubble(Allen)", "Sand(Allen)")
  HAB <- as.data.frame(HAB)
  HAB$species_name <- only_model$species_name
  HAB$cat <- rep("HAB", nrow(HAB))
  HAB <- HAB %>% inner_join(sp_car, "species_name")
  
  final_object <- list(ENV, SOC, HAB)
  names(final_object) <- c("ENV", "SOC", "HAB")
  final_object
  
}, mc.cores = 1)
names(ext_cont_biomass) <- plot_level

ggplot(ext_cont_biomass[[6]]$ENV, aes(x = log10(MaxLength+1), y = mean_chl_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[6]]$ENV, aes(x = Water.column, y = mean_chl_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[6]]$ENV, aes(x = Habitat, y = mean_chl_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[6]]$ENV, aes(x = Trophic_guild_name, y = mean_chl_1year)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[6]]$ENV, aes(x = log10(MaxLength+1), y = min_sst_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[6]]$ENV, aes(x = Water.column, y = min_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[6]]$ENV, aes(x = Habitat, y = min_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[6]]$ENV, aes(x = Trophic_guild_name, y = min_sst_1year)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[4]]$ENV, aes(x = log10(MaxLength+1), y = mean_chl_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[4]]$ENV, aes(x = Water.column, y = mean_chl_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[4]]$ENV, aes(x = Habitat, y = mean_chl_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[4]]$ENV, aes(x = Trophic_guild_name, y = mean_chl_1year)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[4]]$ENV, aes(x = log10(MaxLength+1), y = min_sst_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[4]]$ENV, aes(x = Water.column, y = min_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[4]]$ENV, aes(x = Habitat, y = min_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[4]]$ENV, aes(x = Trophic_guild_name, y = min_sst_1year)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[1]]$ENV, aes(x = log10(MaxLength+1), y = mean_chl_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[1]]$ENV, aes(x = Water.column, y = mean_chl_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[1]]$ENV, aes(x = Habitat, y = mean_chl_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[1]]$ENV, aes(x = Trophic_guild_name, y = mean_chl_1year)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[1]]$ENV, aes(x = log10(MaxLength+1), y = min_sst_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[1]]$ENV, aes(x = Water.column, y = min_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[1]]$ENV, aes(x = Habitat, y = min_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[1]]$ENV, aes(x = Trophic_guild_name, y = min_sst_1year)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[3]]$ENV, aes(x = log10(MaxLength+1), y = mean_sss_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[3]]$ENV, aes(x = Water.column, y = mean_sss_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[3]]$ENV, aes(x = Habitat, y = mean_sss_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[3]]$ENV, aes(x = Trophic_guild_name, y = mean_sss_1year)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[3]]$ENV, aes(x = log10(MaxLength+1), y = min_sst_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[3]]$ENV, aes(x = Water.column, y = min_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[3]]$ENV, aes(x = Habitat, y = min_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[3]]$ENV, aes(x = Trophic_guild_name, y = min_sst_1year)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[2]]$ENV, aes(x = log10(MaxLength+1), y = mean_sss_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[2]]$ENV, aes(x = Water.column, y = mean_sss_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[2]]$ENV, aes(x = Habitat, y = mean_sss_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[2]]$ENV, aes(x = Trophic_guild_name, y = mean_sss_1year)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[2]]$ENV, aes(x = log10(MaxLength+1), y = min_sst_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[2]]$ENV, aes(x = Water.column, y = min_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[2]]$ENV, aes(x = Habitat, y = min_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[2]]$ENV, aes(x = Trophic_guild_name, y = min_sst_1year)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[2]]$ENV, aes(x = log10(MaxLength+1), y = max_sst_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[2]]$ENV, aes(x = Water.column, y = max_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[2]]$ENV, aes(x = Habitat, y = max_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[2]]$ENV, aes(x = Trophic_guild_name, y = max_sst_1year)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[2]]$ENV, aes(x = log10(MaxLength+1), y = mean_pH_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[2]]$ENV, aes(x = Water.column, y = mean_pH_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[2]]$ENV, aes(x = Habitat, y = mean_pH_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[2]]$ENV, aes(x = Trophic_guild_name, y = mean_pH_1year)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[5]]$ENV, aes(x = log10(MaxLength+1), y = mean_sss_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[5]]$ENV, aes(x = Water.column, y = mean_sss_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[5]]$ENV, aes(x = Habitat, y = mean_sss_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[5]]$ENV, aes(x = Trophic_guild_name, y = mean_sss_1year)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[5]]$ENV, aes(x = log10(MaxLength+1), y = min_sst_1year)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[5]]$ENV, aes(x = Water.column, y = min_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[5]]$ENV, aes(x = Habitat, y = min_sst_1year)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[5]]$ENV, aes(x = Trophic_guild_name, y = min_sst_1year)) +
  geom_boxplot()





ggplot(ext_cont_biomass[[6]]$HAB, aes(x = log10(MaxLength+1), y = Depth)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[6]]$HAB, aes(x = Water.column, y = Depth)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[6]]$HAB, aes(x = Habitat, y = Depth)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[6]]$HAB, aes(x = Trophic_guild_name, y = Depth)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[6]]$HAB, aes(x = log10(MaxLength+1), y = Coral)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[6]]$HAB, aes(x = Water.column, y = Coral)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[6]]$HAB, aes(x = Habitat, y = Coral)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[6]]$HAB, aes(x = Trophic_guild_name, y = Coral)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[6]]$HAB, aes(x = log10(MaxLength+1), y = `Reef extent(Allen)`)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[6]]$HAB, aes(x = Water.column, y = `Reef extent(Allen)`)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[6]]$HAB, aes(x = Habitat, y = `Reef extent(Allen)`)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[6]]$HAB, aes(x = Trophic_guild_name, y = `Reef extent(Allen)`)) +
  geom_boxplot()


ggplot(ext_cont_biomass[[4]]$HAB, aes(x = log10(MaxLength+1), y = Depth)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[4]]$HAB, aes(x = Water.column, y = Depth)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[4]]$HAB, aes(x = Habitat, y = Depth)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[4]]$HAB, aes(x = Trophic_guild_name, y = Depth)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[4]]$HAB, aes(x = log10(MaxLength+1), y = Coral)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[4]]$HAB, aes(x = Water.column, y = Coral)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[4]]$HAB, aes(x = Habitat, y = Coral)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[4]]$HAB, aes(x = Trophic_guild_name, y = Coral)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[4]]$HAB, aes(x = log10(MaxLength+1), y = `Reef extent(Allen)`)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[4]]$HAB, aes(x = Water.column, y = `Reef extent(Allen)`)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[4]]$HAB, aes(x = Habitat, y = `Reef extent(Allen)`)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[4]]$HAB, aes(x = Trophic_guild_name, y = `Reef extent(Allen)`)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[1]]$HAB, aes(x = log10(MaxLength+1), y = Coral)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[1]]$HAB, aes(x = Water.column, y = Coral)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[1]]$HAB, aes(x = Habitat, y = Coral)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[1]]$HAB, aes(x = Trophic_guild_name, y = Coral)) +
  geom_boxplot()


ggplot(ext_cont_biomass[[6]]$SOC, aes(x = log10(MaxLength+1), y = Gravity)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[6]]$SOC, aes(x = Water.column, y = Gravity)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[6]]$SOC, aes(x = Habitat, y = Gravity)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[6]]$SOC, aes(x = Trophic_guild_name, y = Gravity)) +
  geom_boxplot()


ggplot(ext_cont_biomass[[4]]$SOC, aes(x = log10(MaxLength+1), y = Gravity)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[4]]$SOC, aes(x = Water.column, y = Gravity)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[4]]$SOC, aes(x = Habitat, y = Gravity)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[4]]$SOC, aes(x = Trophic_guild_name, y = Gravity)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[3]]$SOC, aes(x = log10(MaxLength+1), y = `Natural Ressource`)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[3]]$SOC, aes(x = Water.column, y = `Natural Ressource`)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[3]]$SOC, aes(x = Habitat, y = `Natural Ressource`)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[3]]$SOC, aes(x = Trophic_guild_name, y = `Natural Ressource`)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[2]]$SOC, aes(x = log10(MaxLength+1), y = `MPA Effectiveness`)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[2]]$SOC, aes(x = Water.column, y = `MPA Effectiveness`)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[2]]$SOC, aes(x = Habitat, y = `MPA Effectiveness`)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[2]]$SOC, aes(x = Trophic_guild_name, y = `MPA Effectiveness`)) +
  geom_boxplot()

ggplot(ext_cont_biomass[[5]]$SOC, aes(x = log10(MaxLength+1), y = `MPA Effectiveness`)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(ext_cont_biomass[[5]]$SOC, aes(x = Water.column, y = `MPA Effectiveness`)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[5]]$SOC, aes(x = Habitat, y = `MPA Effectiveness`)) +
  geom_boxplot()
ggplot(ext_cont_biomass[[5]]$SOC, aes(x = Trophic_guild_name, y = `MPA Effectiveness`)) +
  geom_boxplot()


