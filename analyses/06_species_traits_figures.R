# source functions ----

source("R/06_species_traits_figures_function.R")

pal_sp_trait <- PNWColors::pnw_palette("Bay",3, type = "discrete")

load("outputs/best_models.Rdata")

#### Covariates contribution plot ####

# Merge all files together by species

list_files_path <- list.files("outputs/biomass_contribution", full.names = T)
bind_files <- lapply(1:length(list_files_path), function(i) {
  
  load(list_files_path[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
bind_files <- do.call(rbind, bind_files)

fitted_model <- unique(bind_files$fitted_model)

# Load species traits
sp_car <- read.csv("data/new_raw_data/Traits_tropical_spp_1906.csv", header = TRUE)

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
