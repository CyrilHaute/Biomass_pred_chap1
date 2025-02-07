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

covariates_importance_all_function <- function(plot_data,
                                               title,
                                               legend.position,
                                               title.size,
                                               axis.text.x,
                                               axis.text.y,
                                               axis.title,
                                               legend.text,
                                               strip.text.x,
                                               strip.text.y,
                                               fill)
  
{
  
  require(ggplot2)
  
  # covariates relative importance by mean
  
  only_model_best <- best_models |>
    dplyr::inner_join(bind_files)
  
  only_model_best <- only_model_best[only_model_best$best_model == only_model_best$fitted_model,]
  
  ENV <- only_model_best |> 
    dplyr::filter(variable %in% env_var) |> 
    dplyr::mutate(VAR = "ENV")
  ENV[ENV$variable == "max_1year_analysed_sst",]$variable <- "Sea Surface Temperature (max 1 year)"
  ENV[ENV$variable == "min_1year_analysed_sst",]$variable <- "Sea Surface Temperature (min 1 year)"
  ENV[ENV$variable == "max_5year_degree_heating_week",]$variable <- "Degree Heating Week (max 5 year)"
  ENV[ENV$variable == "mean_1year_nppv",]$variable <- "Net Primary Productivity (mean 1 year)"
  ENV[ENV$variable == "mean_1year_so_mean",]$variable <- "Sea Surface Salinity (mean 1 year)"
  ENV[ENV$variable == "min_5year_ph",]$variable <- "pH (min 5 year)"
  
  SOC <- only_model_best |> 
    dplyr::filter(variable %in% hum_var) |> 
    dplyr::mutate(VAR = "HUM")
  SOC[SOC$variable == "protection_status2",]$variable <- "MPA status"
  SOC[SOC$variable == "gdp",]$variable <- "Gross Domestic Product"
  SOC[SOC$variable == "gravtot2",]$variable <- "Human Gravity"
  SOC[SOC$variable == "n_fishing_vessels",]$variable <- "Fishing Vessels Density"
  SOC[SOC$variable == "neartt",]$variable <- "Nearest Population"
  SOC[SOC$variable == "marine_ecosystem_dependency",]$variable <- "Marine Ecosystem Dependency"
  
  HAB <- only_model_best |> 
    dplyr::filter(variable %in% hab_var) |> 
    dplyr::mutate(VAR = "HAB")
  HAB[HAB$variable == "Rock_500m",]$variable <- "Rock (%)"
  HAB[HAB$variable == "Sand_500m",]$variable <- "Sand (%)"
  HAB[HAB$variable == "coral_algae_500m",]$variable <- "Coral/Algae (%)"
  HAB[HAB$variable == "coral",]$variable <- "Coral (RLS)"
  HAB[HAB$variable == "depth",]$variable <- "Depth"
  HAB[HAB$variable == "reef_extent",]$variable <- "Reef extent"
  
  cont <- ENV |>  
    dplyr::full_join(HAB) |> 
    dplyr::full_join(SOC)
  
  importance_plot <- ggplot(cont) +
    geom_point(aes(x = reorder(variable, Dropout_loss, FUN = median), y = Dropout_loss, fill = VAR),
               position = "jitter", 
               size = 0.2) +
    geom_boxplot(aes(x = reorder(variable, Dropout_loss, FUN = median), y = Dropout_loss, fill = VAR),
                 outlier.shape = NA,
                 alpha = 0.9,
                 size = 0.3) +
    scale_y_continuous(limits = c(0, quantile(cont[cont$variable == "Net Primary Productivity (mean 1 year)",]$Dropout_loss, probs = 0.95))) + 
    scale_fill_manual(values = c("ENV" = pal_contribution[2],
                                 "HUM" = pal_contribution[1],
                                 "HAB" = pal_contribution[13])) +
    theme_minimal() +
    coord_flip() +
    # facet_wrap(~ realm) +
    labs(y = "Relative importance (RMSE)", x = "", fill = fill, title = title) +
    theme(legend.position = legend.position) +
    theme(title = element_text(size = title.size),
          axis.text.x = element_text(size = axis.text.x),
          axis.text.y = element_text(size = axis.text.y),
          axis.title = element_text(axis.title),
          legend.text = element_text(size = legend.text),
          strip.text.x = element_text(size = strip.text.x),
          strip.text.y = element_text(size = strip.text.y))
  
}

merged_covariates_importance_all_function <- function(plot_data,
                                                      title,
                                                      legend.position,
                                                      title.size,
                                                      axis.text.x,
                                                      axis.text.y,
                                                      axis.title,
                                                      legend.text,
                                                      strip.text.x,
                                                      strip.text.y,
                                                      geom.text.size,
                                                      fill,
                                                      ylim){
  
  require(ggplot2)
  
  # covariates relative importance by median
  
  only_model_best_best <- best_models |>
    dplyr::inner_join(bind_files)
  
  only_model_best_best <- only_model_best_best[only_model_best_best$best_model == only_model_best_best$fitted_model,]
  
  ENV <- only_model_best_best |> 
    dplyr::filter(variable %in% env_var) |> 
    dplyr::mutate(VAR = "ENV")
  
  SOC <- only_model_best_best |> 
    dplyr::filter(variable %in% hum_var) |> 
    dplyr::mutate(VAR = "HUM")
  
  HAB <- only_model_best_best |> 
    dplyr::filter(variable %in% hab_var) |> 
    dplyr::mutate(VAR = "HAB")
  
  cont <- ENV |>  
    dplyr::full_join(HAB) |> 
    dplyr::full_join(SOC)
  
  var_max <- var_max_all_function(plot_data = bind_files)
  
  y <- cont |> 
    dplyr::group_by(VAR) |> 
    dplyr::summarise(y = quantile(Dropout_loss, probs = 0.96)) 
  y <- y |> 
    dplyr::mutate(y = max(y))
  var_max <- var_max |> 
    dplyr::inner_join(y)
  #merge contribution per var and model
  
  cont_merge <- cont |>  
    dplyr::group_by(species_name, VAR) |>
    dplyr::summarise(value = mean(Dropout_loss),
                     sd = mean(sd_dropout_loss))
  
  cont_merge <- dplyr::full_join(cont_merge, var_max, by = c("VAR"))
  cont_merge[which(is.na(cont_merge$n)),]$n <- 0
  
  merged_importance_plot <- ggplot(cont_merge) +
    geom_point(aes(x = reorder(VAR, value, FUN = mean), y = value, fill = VAR),
               size = 0.2,
               position = "jitter") +
    geom_boxplot(aes(x = reorder(VAR, value, FUN = mean), y = value, fill = VAR),
                 outlier.shape = NA,
                 alpha = 0.9,
                 size = 0.3) +
    geom_text(data = var_max, aes(x = VAR, y = y, label = n), size = geom.text.size) +
    scale_y_continuous(limits = c(0, 0.3)) + 
    scale_fill_manual(values = c("ENV" = pal_contribution[2],
                                 "HUM" = pal_contribution[1],
                                 "HAB" = pal_contribution[13])) +
    theme_minimal() +
    coord_flip() +
    labs(y = "Relative importance (RMSE)", x = "", fill = fill, title = title) +
    theme(legend.position = legend.position) +
    theme(title = element_text(size = title.size),
          axis.text.x = element_text(size = axis.text.x),
          axis.text.y = element_text(size = axis.text.y),
          axis.title = element_text(size = axis.title),
          legend.text = element_text(size = legend.text),
          strip.text.x = element_text(size = strip.text.x),
          strip.text.y = element_text(size = strip.text.y))
  
}


