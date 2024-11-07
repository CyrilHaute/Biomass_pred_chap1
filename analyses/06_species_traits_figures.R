# source functions ----

library(patchwork)

source("R/06_species_traits_figures_function.R")

pal_sp_trait <- PNWColors::pnw_palette("Bay", 3, type = "discrete")

load("outputs/best_models.Rdata")

#### Covariates contribution plot ####

# Merge all files together by species

glm <- list.files("outputs/glm_contribution", full.names = TRUE)
glm <- lapply(1:length(glm), function(i) {
  
  load(glm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
glm <- do.call(rbind, glm)

gam <- list.files("outputs/gam_contribution", full.names = TRUE)
gam <- lapply(1:length(gam), function(i) {
  
  load(gam[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gam <- do.call(rbind, gam)

spamm <- list.files("outputs/spamm_contribution", full.names = TRUE)
spamm <- lapply(1:length(spamm), function(i) {
  
  load(spamm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
spamm <- do.call(rbind, spamm)

rf <- list.files("outputs/rf_contribution", full.names = TRUE)
rf <- lapply(1:length(rf), function(i) {
  
  load(rf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
rf <- do.call(rbind, rf)

sprf <- list.files("outputs/sprf_contribution", full.names = TRUE)
sprf <- lapply(1:length(sprf), function(i) {
  
  load(sprf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
sprf <- do.call(rbind, sprf)
colnames(sprf)[colnames(sprf) %in% c("global_dropout_loss", "global_sd_dropout_loss")] <- c("Dropout_loss", "sd_dropout_loss")

gbm <- list.files("outputs/brt_contribution", full.names = TRUE)
gbm <- lapply(1:length(gbm), function(i) {
  
  load(gbm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gbm <- do.call(rbind, gbm)

bind_files <- list(glm, gam, spamm, rf, sprf, gbm)
bind_files <- purrr::reduce(bind_files, dplyr::full_join)

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
                                           labs_title = "A. Maximum length (cm)",
                                           aes_string_x = c(3.1, 4.85, 4.1, 2.9, 5.12, 2.12, 3.9, 1.9, 0.85, 1.12, 5.37, 4.37, 3.35, 2.3, 1.3),
                                           x_text_angle = NULL,
                                           legend.position = "none")

plot_water.column <- species_traits_function(plot_data = bind_files,
                                             trait = "Water.column",
                                             color = pal_sp_trait,
                                             labs_title = "B. Water column position ",
                                             aes_string_x = c(3.05, 2.8, 2.05, 1.8, 1.05, 3.3, 0.8, 2.3, 1.3),
                                             x_text_angle = NULL,
                                             legend.position = "none")

plot_habitat <- species_traits_function(plot_data = bind_files,
                                        trait = "Habitat",
                                        color = pal_sp_trait,
                                        labs_title = "C. Habitat",
                                        aes_string_x = c(4.05, 3.85, 2.1, 3.1, 1.825, 0.85, 1.1, 2.85, 1.3, 4.3, 3.3, 2.3),
                                        x_text_angle = NULL,
                                        legend.position = "none")

plot_trophic <- species_traits_function(plot_data = bind_files,
                                        trait = "Trophic_guild_name",
                                        color = pal_sp_trait,
                                        labs_title = "D. Trophic classes",
                                        aes_string_x = c(5.82, 2.82, 6.1, 3.1, 2.1, 4.05, 8.15, 3.83, 1.85, 6.85, 1.1, 7.87, 4.87, 5.13, 0.88, 7.15, 3.4, 7.45, 4.45, 6.45, 2.35, 8.4, 5.4, 1.3),
                                        x_text_angle = 65,
                                        legend.position = c(0.85, 0.8))

plot_species_traits <- plot_max.length / plot_water.column / plot_habitat / plot_trophic

ggsave("figures/plot_species_traits.pdf", plot_species_traits, height = 25, width = 19)
ggsave("figures/plot_species_traits_presentation.png", plot_species_traits, height = 25, width = 19)

plot_species_traits <- (plot_max.length / plot_water.column)
ggsave("figures/plot_species_traits_presentation.png", plot_species_traits, height = 13, width = 16)
plot_species_traits2 <- (plot_habitat / plot_trophic)
ggsave("figures/plot_species_traits_presentation2.png", plot_species_traits2, height = 13, width = 16)




cont <- lapply(1:length(fitted_model), function(i) {
  
  plot_level <- fitted_model[i]
  plot_data <- plot_data |> 
    dplyr::filter(fitted_model == plot_level)
  plot_data <- plot_data |>
    dplyr::filter(species_name %in% best_models[best_models$best_model == plot_level,1]$species_name)
  
  # ENV <- lapply(1:nrow(only_model), function(i) { only_model$contributions_and_sd[[i]][only_model$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$Dropout_loss})
  # ENV_sd <- lapply(1:nrow(only_model), function(i) { only_model$contributions_and_sd[[i]][only_model$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$sd_dropout_loss})
  # ENV <- do.call(rbind, ENV)
  # ENV_sd <- do.call(rbind, ENV_sd)
  # ENV <- dplyr::tibble(species_name = only_model$species_name,
  #                      value = matrixStats::rowMedians(ENV),
  #                      sd = matrixStats::rowMedians(ENV_sd),
  #                      var = rep("ENV",nrow(ENV)),
  #                      plot_level = rep(plot_level, nrow(ENV)))
  ENV <- plot_data |> 
    dplyr::filter(variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "max_5year_nppv", "min_1year_analysed_sst", "min_5year_ph", "min_7days_o2")) |> 
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "ENV",
                     plot_level = plot_level)
  
  # SOC <- lapply(1:nrow(only_model), function(i) { only_model$contributions_and_sd[[i]][only_model$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$Dropout_loss})
  # SOC_sd <- lapply(1:nrow(only_model), function(i) { only_model$contributions_and_sd[[i]][only_model$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$sd_dropout_loss})
  # SOC <- do.call(rbind, SOC)
  # SOC_sd <- do.call(rbind, SOC_sd)
  # SOC <- dplyr::tibble(species_name = only_model$species_name,
  #                      value = matrixStats::rowMedians(SOC),
  #                      sd = matrixStats::rowMedians(SOC_sd),
  #                      var = rep("HUM",nrow(SOC)),
  #                      plot_level = rep(plot_level, nrow(SOC)))
  
  SOC <- plot_data |> 
    dplyr::filter(variable %in% c("effectiveness", "gdp", "gravtot2", "no_violence", "n_fishing_vessels", "natural_ressource_rent", "neartt", "marine_ecosystem_dependency")) |> 
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "HUM",
                     plot_level = plot_level)
  
  # HAB <- lapply(1:nrow(only_model), function(i) { only_model$contributions_and_sd[[i]][only_model$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$Dropout_loss})
  # HAB_sd <- lapply(1:nrow(only_model), function(i) { only_model$contributions_and_sd[[i]][only_model$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$sd_dropout_loss})
  # HAB <- do.call(rbind, HAB)
  # HAB_sd <- do.call(rbind, HAB_sd)
  # HAB <- dplyr::tibble(species_name = only_model$species_name,
  #                      value = matrixStats::rowMedians(HAB),
  #                      sd = matrixStats::rowMedians(HAB_sd),
  #                      var = rep("HAB",nrow(HAB)),
  #                      plot_level = rep(plot_level,nrow(HAB)))
  
  HAB <- plot_data |> 
    dplyr::filter(variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "seagrass", "depth", "reef_extent")) |> 
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "HAB",
                     plot_level = plot_level)
  
  cont <- ENV |> 
    dplyr::full_join(HAB) |> 
    dplyr::full_join(SOC)
  cont$plot_level <- rep(plot_level, nrow(cont))
  cont
  
})

cont <- do.call(rbind, cont)

load("data/new_raw_data/RLS_actinopterygii_data.Rdata")
phylo <- RLS_actinopterygii_data |> 
  dplyr::select(species_name, order, family)
phylo <- unique(phylo)

cont <- cont |> 
  dplyr::inner_join(phylo, by = "species_name")

n_trait <- cont |>
  dplyr::group_by(family) |> 
  dplyr::summarise(n = dplyr::n() / 3)

cont <- cont |> 
  dplyr::inner_join(n_trait, by = "family") |> 
  dplyr::filter(n >= 10) |> 
  dplyr::select(-n)

library(ggplot2)

plot_trait <- cont |> 
  dplyr::mutate(VAR = forcats::fct_relevel(VAR, "ENV", "HAB", "HUM")) |>
  ggplot(aes_string(x = "family", y = "value", fill = "VAR")) +
  geom_boxplot(aes(fill = factor(VAR)), width=0.6, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("ENV" = pal_sp_trait[1],
                               "HUM" = pal_sp_trait[3],
                               "HAB" = pal_sp_trait[2])) +
  theme_bw() +
  coord_cartesian(ylim = c(0,0.25)) +
  labs(y = "Relative importance (RMSE)", x = "", fill = "") +
  theme(
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white", color = NA),
    title = element_text(size = 20),
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 20),
    strip.text.x = element_text(size = 20),
    strip.text.y = element_text(size = 20),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = "grey50",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

ggsave("figures/plot_family.png", plot_trait, height = 10, width = 20)

scaridae <- cont |> 
  dplyr::filter(family == "Scaridae")

scaridae <- scaridae |> 
  dplyr::mutate(var_reordered = tidytext::reorder_within(VAR, value, species_name))

plot_scaridae <- ggplot(scaridae) +
  geom_col(aes(x = var_reordered, y = value, fill = VAR)) +
  geom_errorbar(aes(x = var_reordered, y = value, ymin = value - sd, ymax = value + sd), width = .1, position = position_dodge(.9)) +
  theme_bw() +
  coord_flip() +
  facet_wrap(~species_name, scales = "free_y") +
  tidytext::scale_x_reordered() +
  scale_fill_manual(values = c("ENV" = pal_sp_trait[1], 
                               "HUM" = pal_sp_trait[3], 
                               "HAB" = pal_sp_trait[2])) +
  labs(y = "Relative importance (RMSE)", x = "", fill = "") +
  theme(
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white", color = NA),
    title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = "grey50",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

ggsave("figures/plot_scaridae.png", plot_scaridae, height = 10, width = 20)
