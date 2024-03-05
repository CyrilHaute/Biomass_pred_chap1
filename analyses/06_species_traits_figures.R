# source functions ----

library(patchwork)

source("R/06_species_traits_figures_function.R")

pal_sp_trait <- PNWColors::pnw_palette("Bay", 3, type = "discrete")

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
sp_car <- read.csv("data/new_raw_data/Traits_tropical_spp_1906.csv", header = TRUE) |> 
  dplyr::rename(species_name = Species) |> 
  dplyr::filter(species_name %in% unique(bind_files$species_name))
sp_car$ML_cat <- NA
sp_car <- sp_car[which(is.na(sp_car$MaxLength) == FALSE),]
sp_car <- sp_car[which(is.na(sp_car$Trophic_guild_name) == FALSE),]

sp_car[sp_car$MaxLength > 0 & sp_car$MaxLength <= 20,]$ML_cat <- "0-20 cm"
sp_car[sp_car$MaxLength > 20 & sp_car$MaxLength <= 40,]$ML_cat <- "20-40 cm"
sp_car[sp_car$MaxLength > 40 & sp_car$MaxLength <= 60,]$ML_cat <- "40-60 cm"
sp_car[sp_car$MaxLength > 60 & sp_car$MaxLength <= 80,]$ML_cat <- "60-80 cm"
sp_car[sp_car$MaxLength > 80 & sp_car$MaxLength <= 300,]$ML_cat <- "80-300 cm"

sp_car[sp_car$Water.column == "Demersal",]$Water.column <- "demersal"
sp_car[sp_car$Water.column == "pelagic non-site attached",]$Water.column <- "pelagic"
sp_car[sp_car$Water.column == "pelagic site attached",]$Water.column <- "pelagic"

sp_car[sp_car$Habitat == "Coral",]$Habitat <- "coral"

sp_car[sp_car$Trophic_guild_name == "Herbivores Microvores Detritivores",]$Trophic_guild_name <- "herbivores"

plot_max.length <- species_traits_function(plot_data = bind_files,
                                           trait = "ML_cat",
                                           color = pal_sp_trait,
                                           labs_title = "A",
                                           aes_string_x = c(3.1, 4.85, 4.1, 2.9, 5.12, 2.12, 3.9, 1.9, 0.85, 1.12, 5.37, 4.37, 3.35, 2.3, 1.3),
                                           x_text_angle = NULL)

plot_water.column <- species_traits_function(plot_data = bind_files,
                                             trait = "Water.column",
                                             color = pal_sp_trait,
                                             labs_title = "B",
                                             aes_string_x = c(3.05, 2.8, 2.05, 1.8, 1.05, 3.3, 0.8, 2.3, 1.3),
                                             x_text_angle = NULL)

plot_habitat <- species_traits_function(plot_data = bind_files,
                                        trait = "Habitat",
                                        color = pal_sp_trait,
                                        labs_title = "C",
                                        aes_string_x = c(4.05, 3.85, 2.1, 3.1, 1.825, 0.85, 1.1, 2.85, 1.3, 4.3, 3.3, 2.3),
                                        x_text_angle = NULL)

plot_trophic <- species_traits_function(plot_data = bind_files,
                                        trait = "Trophic_guild_name",
                                        color = pal_sp_trait,
                                        labs_title = "D",
                                        aes_string_x = c(5.82, 2.82, 6.1, 3.1, 2.1, 4.05, 8.15, 3.83, 1.85, 6.85, 1.1, 7.87, 4.87, 5.13, 0.88, 7.15, 3.4, 7.45, 4.45, 6.45, 2.35, 8.4, 5.4, 1.3),
                                        x_text_angle = 65)

plot_species_traits <- plot_max.length / plot_water.column / plot_habitat / plot_trophic

ggsave("figures/plot_species_traits.pdf", plot_species_traits, height = 25, width = 19)
ggsave("figures/plot_species_traits.png", plot_species_traits, height = 25, width = 19)
