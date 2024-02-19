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

sp_car[sp_car$MaxLength > 0 & sp_car$MaxLength <= 10,]$ML_cat <- "0-10 cm"
sp_car[sp_car$MaxLength > 10 & sp_car$MaxLength <= 20,]$ML_cat <- "10-20 cm"
sp_car[sp_car$MaxLength > 20 & sp_car$MaxLength <= 30,]$ML_cat <- "20-30 cm"
sp_car[sp_car$MaxLength > 30 & sp_car$MaxLength <= 40,]$ML_cat <- "30-40 cm"
sp_car[sp_car$MaxLength > 40 & sp_car$MaxLength <= 50,]$ML_cat <- "40-50 cm"
sp_car[sp_car$MaxLength > 50 & sp_car$MaxLength <= 60,]$ML_cat <- "50-60 cm"
sp_car[sp_car$MaxLength > 60 & sp_car$MaxLength <= 70,]$ML_cat <- "60-70 cm"
sp_car[sp_car$MaxLength > 70 & sp_car$MaxLength <= 80,]$ML_cat <- "70-80 cm"
sp_car[sp_car$MaxLength > 80 & sp_car$MaxLength <= 300,]$ML_cat <- "80-300 cm"

sp_car[sp_car$Water.column == "Demersal",]$Water.column <- "demersal"
sp_car[sp_car$Water.column == "pelagic non-site attached",]$Water.column <- "pelagic"
sp_car[sp_car$Water.column == "pelagic site attached",]$Water.column <- "pelagic"

sp_car[sp_car$Habitat == "Coral",]$Habitat <- "coral"

sp_car[sp_car$Trophic_guild_name == "Herbivores Microvores Detritivores",]$Trophic_guild_name <- "herbivores"

plot_max.length <- species_traits_function(plot_data = bind_files,
                                           trait = "ML_cat",
                                           color = pal_sp_trait,
                                           labs_title = "A")

plot_water.column <- species_traits_function(plot_data = bind_files,
                                             trait = "Water.column",
                                             color = pal_sp_trait,
                                             labs_title = "B")

plot_habitat <- species_traits_function(plot_data = bind_files,
                                        trait = "Habitat",
                                        color = pal_sp_trait,
                                        labs_title = "C")

plot_trophic <- species_traits_function(plot_data = bind_files,
                                        trait = "Trophic_guild_name",
                                        color = pal_sp_trait,
                                        labs_title = "D")

plot_species_traits <- plot_max.length / plot_water.column / plot_habitat / plot_trophic

ggsave("figures/plot_species_traits.pdf", plot_species_traits, height = 25, width = 19)
