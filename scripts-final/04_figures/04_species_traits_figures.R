# load in packages ---- 

libs <- c('tidyverse', 'gridExtra', 'ggplot2', 'patchwork', 'matrixStats', 'parallel','PNWColors', 'agricolae')
lapply(libs, library, character.only = T, lib.loc = '/home/marbec/R/x86_64-pc-linux-gnu-library/4.1')

# check all packages are loaded
if(sum(libs %in% (.packages())) != length(libs)){
  stop('packages not loaded correctly')}

# source functions ----

source("scripts-final/00_functions/species_traits_figures_function.R")

pal_sp_trait = pnw_palette("Bay",3, type = "discrete")

best_models <- readRDS("results/overall_best_models_R3.rds")

#### Covariates contribution plot ####

# Merge all files together by species

path <- "results/contributions/contributions_biomass"

# load in abundance data
# Save directory as character object
directory <- path

# Extract names all elements in folder
all_files <- list.files(path = directory)

Contributions_biomass <- mclapply(1:length(all_files), function(i) {
  all_files_full <- list()
  all_files_full[i] <- paste0(directory, '/', all_files[i])
  readRDS(all_files_full[[i]])
},mc.cores = 1)

Contributions_biomass <- do.call(rbind, Contributions_biomass)

fitted_model <- unique(Contributions_biomass$fitted_model)

sp_car <- readRDS("data/Cyril_data/RLS_species_traits_new.rds")

plot_max.length <- species_traits_function(plot_data = Contributions_biomass,
                                           trait = "ML_cat",
                                           color = pal_sp_trait,
                                           labs_title = "A")

plot_water.column <- species_traits_function(plot_data = Contributions_biomass,
                                             trait = "Water.column",
                                             color = pal_sp_trait,
                                             labs_title = "B")

plot_habitat <- species_traits_function(plot_data = Contributions_biomass,
                                        trait = "Habitat",
                                        color = pal_sp_trait,
                                        labs_title = "C")

plot_trophic <- species_traits_function(plot_data = Contributions_biomass,
                                        trait = "Trophic_guild_name",
                                        color = pal_sp_trait,
                                        labs_title = "D")

plot_species_traits <- plot_max.length / plot_water.column / plot_habitat / plot_trophic

ggsave("figures-R3/plot_species_traits.png", plot_species_traits, height = 25, width = 19)
