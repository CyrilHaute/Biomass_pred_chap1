# source functions ----

library(patchwork)
library(ggplot2)

source("R/05_contributions_figures_functions.R")
source("R/05_load_realm_contribution_function.R")

pal_contribution <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1"), PNWColors::pnw_palette("Bay", 6, type = "continuous"))

load("outputs/best_models.Rdata")

env_var <- c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_nppv", "mean_1year_so_mean", "min_1year_analysed_sst", "min_5year_ph")
hum_var <- c("protection_status2", "gdp", "gravtot2", "n_fishing_vessels", "neartt", "marine_ecosystem_dependency")
hab_var <- c("Rock_500m", "Sand_500m", "coral", "coral_algae_500m", "depth", "reef_extent")

#################### Global contribution ####################


#### Covariates contribution plot ####

glm <- list.files("outputs/glm_contribution", full.names = TRUE)
glm <- lapply(1:length(glm), function(i) {
  
  load(glm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
glm <- do.call(rbind, glm)
glm_inf <- list.files("outputs/glm_contribution_benthos_infered", full.names = TRUE)
glm_inf <- lapply(1:length(glm_inf), function(i) {
  
  load(glm_inf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
glm_inf <- do.call(rbind, glm_inf)

gam <- list.files("outputs/gam_contribution", full.names = TRUE)
gam <- lapply(1:length(gam), function(i) {
  
  load(gam[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gam <- do.call(rbind, gam)
gam_inf <- list.files("outputs/gam_contribution_benthos_infered", full.names = TRUE)
gam_inf <- lapply(1:length(gam_inf), function(i) {
  
  load(gam_inf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gam_inf <- do.call(rbind, gam_inf)

spamm <- list.files("outputs/spamm_contribution", full.names = TRUE)
spamm <- lapply(1:length(spamm), function(i) {
  
  load(spamm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
spamm <- do.call(rbind, spamm)
spamm_inf <- list.files("outputs/spamm_contribution_benthos_infered", full.names = TRUE)
spamm_inf <- lapply(1:length(spamm_inf), function(i) {
  
  load(spamm_inf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
spamm_inf <- do.call(rbind, spamm_inf)

rf <- list.files("outputs/rf_contribution", full.names = TRUE)
rf <- lapply(1:length(rf), function(i) {
  
  load(rf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
rf <- do.call(rbind, rf)
rf_inf <- list.files("outputs/rf_contribution_benthos_infered", full.names = TRUE)
rf_inf <- lapply(1:length(rf_inf), function(i) {
  
  load(rf_inf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
rf_inf <- do.call(rbind, rf_inf)

sprf <- list.files("outputs/sprf_contribution", full.names = TRUE)
sprf <- lapply(1:length(sprf), function(i) {
  
  load(sprf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
sprf <- do.call(rbind, sprf)
colnames(sprf)[colnames(sprf) %in% c("global_dropout_loss", "global_sd_dropout_loss")] <- c("Dropout_loss", "sd_dropout_loss")
sprf_inf <- list.files("outputs/sprf_contribution_benthos_infered", full.names = TRUE)
sprf_inf <- lapply(1:length(sprf_inf), function(i) {
  
  load(sprf_inf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
sprf_inf <- do.call(rbind, sprf_inf)
colnames(sprf_inf)[colnames(sprf_inf) %in% c("global_dropout_loss", "global_sd_dropout_loss")] <- c("Dropout_loss", "sd_dropout_loss")

gbm <- list.files("outputs/brt_contribution", full.names = TRUE)
gbm <- lapply(1:length(gbm), function(i) {
  
  load(gbm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gbm <- do.call(rbind, gbm)
gbm_inf <- list.files("outputs/brt_contribution_benthos_infered", full.names = TRUE)
gbm_inf <- lapply(1:length(gbm_inf), function(i) {
  
  load(gbm_inf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gbm_inf <- do.call(rbind, gbm_inf)

bind_files <- list(glm, gam, spamm, sprf, rf, gbm)
bind_files <- purrr::reduce(bind_files, dplyr::full_join)

bind_files_inf <- list(glm_inf, gam_inf, spamm_inf, sprf_inf, rf_inf, gbm_inf)
bind_files_inf <- purrr::reduce(bind_files_inf, dplyr::full_join)

##################### Contribution per species #####################

# Load species traits
sp_car <- read.csv("data/new_raw_data/Traits_tropical_spp_1906.csv", header = TRUE) |> 
  dplyr::rename(species_name = Species) |> 
  dplyr::filter(species_name %in% unique(bind_files$species_name))

sp_car <- sp_car[which(is.na(sp_car$Trophic_guild_name) == FALSE),]

sp_car[sp_car$Trophic_guild_name == "Herbivores Microvores Detritivores",]$Trophic_guild_name <- "Herbivores"
sp_car[sp_car$Trophic_guild_name == "planktivore",]$Trophic_guild_name <- "Planktivores"
sp_car[sp_car$Trophic_guild_name == "piscivore",]$Trophic_guild_name <- "Piscivores"
sp_car[sp_car$Trophic_guild_name == "microinvertivore",]$Trophic_guild_name <- "Microinvertivores"
sp_car[sp_car$Trophic_guild_name == "macroinvertivore",]$Trophic_guild_name <- "Macroinvertivores"
sp_car[sp_car$Trophic_guild_name == "sessile invertivores",]$Trophic_guild_name <- "Sessile invertivores"
sp_car[sp_car$Trophic_guild_name == "corallivore",]$Trophic_guild_name <- "Corallivores"
sp_car[sp_car$Trophic_guild_name == "crustacivore",]$Trophic_guild_name <- "Crustacivores"

trophic_group <- bind_files |> 
  dplyr::inner_join(sp_car[colnames(sp_car) %in% c("species_name", "Trophic_guild_name")])

trophic_group <- trophic_group |>
  dplyr::inner_join(best_models)

trophic_group <- trophic_group[trophic_group$best_model == trophic_group$fitted_model,]

trophic_group <- trophic_group |> 
  dplyr::inner_join(sp_car[colnames(sp_car) %in% c("species_name", "Trophic_guild_name")]) |> 
  dplyr::mutate(VAR = dplyr::case_when(variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week",
                                                       "mean_1year_nppv", "mean_1year_so_mean", "min_1year_analysed_sst",
                                                       "min_5year_ph") ~ "ENV",
                                       variable %in% c("Rock_500m", "Sand_500m", "coral", "coral_algae_500m", "depth",
                                                       "reef_extent") ~ "HAB",
                                       variable %in% c("gdp", "gravtot2", "marine_ecosystem_dependency", "n_fishing_vessels",            
                                                       "neartt", "protection_status2") ~ "HUM"))
trophic_group <- trophic_group |> 
  dplyr::group_by(Trophic_guild_name, variable, VAR) |> 
  dplyr::summarise(value = median(Dropout_loss),
                   sd = median(sd_dropout_loss))

trophic_group[trophic_group$variable == "max_1year_analysed_sst",]$variable <- "Sea Surface Temperature (max 1 year)"
trophic_group[trophic_group$variable == "min_1year_analysed_sst",]$variable <- "Sea Surface Temperature (min 1 year)"
trophic_group[trophic_group$variable == "max_5year_degree_heating_week",]$variable <- "Degree Heating Week (max 5 year)"
trophic_group[trophic_group$variable == "mean_1year_nppv",]$variable <- "Net Primary Productivity (mean 1 year)"
trophic_group[trophic_group$variable == "mean_1year_so_mean",]$variable <- "Sea Surface Salinity (mean 1 year)"
trophic_group[trophic_group$variable == "min_5year_ph",]$variable <- "pH (min 5 year)"
trophic_group[trophic_group$variable == "protection_status2",]$variable <- "MPA status"
trophic_group[trophic_group$variable == "gdp",]$variable <- "Gross Domestic Product"
trophic_group[trophic_group$variable == "gravtot2",]$variable <- "Human Gravity"
trophic_group[trophic_group$variable == "n_fishing_vessels",]$variable <- "Fishing Vessels Density"
trophic_group[trophic_group$variable == "neartt",]$variable <- "Nearest Population"
trophic_group[trophic_group$variable == "marine_ecosystem_dependency",]$variable <- "Marine Ecosystem Dependency"
trophic_group[trophic_group$variable == "Rock_500m",]$variable <- "Rock (%)"
trophic_group[trophic_group$variable == "Sand_500m",]$variable <- "Sand (%)"
trophic_group[trophic_group$variable == "coral_algae_500m",]$variable <- "Coral/Algae (%)"
trophic_group[trophic_group$variable == "coral",]$variable <- "Coral (RLS)"
trophic_group[trophic_group$variable == "depth",]$variable <- "Depth"
trophic_group[trophic_group$variable == "reef_extent",]$variable <- "Reef extent"

trophic_group$variable <- as.factor(trophic_group$variable)

trophic_group_plot <- trophic_group |> 
  dplyr::mutate(variable = forcats::fct_relevel(variable, "MPA status",
                                                "Fishing Vessels Density",
                                                "Marine Ecosystem Dependency",
                                                "Gross Domestic Product",
                                                "Degree Heating Week (max 5 year)",
                                                "Sea Surface Temperature (max 1 year)",
                                                "pH (min 5 year)",
                                                "Sand (%)",
                                                "Coral/Algae (%)",
                                                "Rock (%)",
                                                "Nearest Population",
                                                "Sea Surface Temperature (min 1 year)",
                                                "Coral (RLS)",
                                                "Reef extent",
                                                "Human Gravity",
                                                "Depth",
                                                "Sea Surface Salinity (mean 1 year)",
                                                "Net Primary Productivity (mean 1 year)")) |>
  ggplot(aes(x = Trophic_guild_name, y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "Reds", direction = 1) +
  labs(y = "", x = "", fill = "Relative importance (RMSE)") +
  theme(
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white", color = NA),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = "grey50",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

ggsave("figures/trophic_group_plot.png", trophic_group_plot, height = 10, width = 16)
