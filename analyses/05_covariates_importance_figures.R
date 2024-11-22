# source functions ----

library(patchwork)

source("R/05_contributions_figures_functions.R")
source("R/05_load_realm_contribution_function.R")

# pal_contribution <- PNWColors::pnw_palette("Bay", type = "discrete")
pal_contribution <- RColorBrewer::brewer.pal(n = 9, name = "Set1")

load("outputs/best_models.Rdata")

#################### Global contribution ####################


#### Covariates contribution plot ####

glm <- list.files("outputs/glm_contribution2", full.names = TRUE)
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

##### For bind_files

covariates_importance_all <- covariates_importance_all_function(plot_data = bind_files,
                                                                title = "A.",
                                                                legend.position = "none",
                                                                title.size = 18,
                                                                axis.text.x = 15,
                                                                axis.text.y = 17,
                                                                axis.title = 21,
                                                                legend.text = 15,
                                                                strip.text.x = 20,
                                                                strip.text.y = 20,
                                                                fill = "")

merged_covariates_importance_all <- merged_covariates_importance_all_function(plot_data = bind_files,
                                                                              title = "B.",
                                                                              legend.position = c(0.85, 0.16),
                                                                              title.size = 18,
                                                                              axis.text.x = 15,
                                                                              axis.text.y = 17,
                                                                              axis.title = 21,
                                                                              legend.text = 15,
                                                                              strip.text.x = 20,
                                                                              strip.text.y = 20,
                                                                              geom.text.size = 5,
                                                                              fill = "")

covariates_importance_all_and_merged <- covariates_importance_all + merged_covariates_importance_all

ggsave("figures/covariates_importance_all_and_merged.png", covariates_importance_all_and_merged, height = 7, width = 15)

covariates_importance_GLM <- covariates_importance_function(plot_data = glm,
                                                            fitted_model = "glm",
                                                            color = pal_contribution,
                                                            labs_y = "",
                                                            labs_fill = "",
                                                            ylim = c(0, 2),
                                                            legend.position = "none")

covariates_importance_GAM <- covariates_importance_function(plot_data = bind_files,
                                                            fitted_model = "gam",
                                                            color = pal_contribution,
                                                            labs_y = "",
                                                            labs_fill = "",
                                                            ylim = c(0, 0.35),
                                                            legend.position = "none")

covariates_importance_SPAMM <- covariates_importance_function(plot_data = bind_files,
                                                              fitted_model = "spamm",
                                                              color = pal_contribution,
                                                              labs_y = "",
                                                              labs_fill = "",
                                                              ylim = c(0, 0.06),
                                                              legend.position = "none")

covariates_importance_RF <- covariates_importance_function(plot_data = rf,
                                                           fitted_model = "rf",
                                                           color = pal_contribution,
                                                           labs_y = "",
                                                           labs_fill = "",
                                                           ylim = c(0, 0.13),
                                                           legend.position = "none")

covariates_importance_GBM <- covariates_importance_function(plot_data = bind_files,
                                                            fitted_model = "gbm",
                                                            color = pal_contribution,
                                                            labs_y = "Relative importance (RMSE)",
                                                            labs_fill = "",
                                                            ylim = c(0, 0.1),
                                                            legend.position = "none")

covariates_importance_SPRF <- covariates_importance_function(plot_data = bind_files,
                                                             fitted_model = "sprf",
                                                             color = pal_contribution,
                                                             labs_y = "Relative importance (RMSE)",
                                                             labs_fill = "",
                                                             ylim = c(0, 0.1),
                                                             legend.position = c(0.75, 0.16))

covariates_importance_all <- (covariates_importance_GLM + covariates_importance_GAM) / (covariates_importance_SPAMM + covariates_importance_RF) / (covariates_importance_GBM + covariates_importance_SPRF)

ggsave("figures/covariates_importance_all.png", covariates_importance_all, height = 16, width = 13)

merged_covariates_importance_GLM <- plot_merged_covariates_importance_function(plot_data =  merged_covariates_importance_function(plot_data = bind_files,
                                                                                                                                  fitted_model = "glm"),
                                                                               color = pal_contribution,
                                                                               labs_y = "",
                                                                               labs_fill = "",
                                                                               legend.position = "none",
                                                                               mul = 6)

merged_covariates_importance_GAM <- plot_merged_covariates_importance_function(plot_data =  merged_covariates_importance_function(plot_data = bind_files,
                                                                                                                                  fitted_model = "gam"),
                                                                               color = pal_contribution,
                                                                               labs_y = "",
                                                                               labs_fill = "",
                                                                               legend.position = "none",
                                                                               mul = 2)

merged_covariates_importance_SPAMM <- plot_merged_covariates_importance_function(merged_covariates_importance_function(plot_data = bind_files,
                                                                                                                       fitted_model = "spamm"),
                                                                                 color = pal_contribution,
                                                                                 labs_y = "",
                                                                                 labs_fill = "",
                                                                                 legend.position = "none",
                                                                                 mul = 2)

merged_covariates_importance_RF <- plot_merged_covariates_importance_function(merged_covariates_importance_function(plot_data = rf,
                                                                                                                    fitted_model = "rf"),
                                                                              color = pal_contribution,
                                                                              # labs_y = "",
                                                                              labs_y = "Relative importance (RMSE)",
                                                                              labs_fill = "",
                                                                              # legend.position = "none",
                                                                              legend.position = c(0.8, 0.18),
                                                                              mul = 3)
ggsave("figures/merged_covariates_importance_rf_biot.pdf", merged_covariates_importance_RF)

merged_covariates_importance_GBM <- plot_merged_covariates_importance_function(merged_covariates_importance_function(plot_data = bind_files,
                                                                                                                     fitted_model = "gbm"),
                                                                               color = pal_contribution,
                                                                               labs_y = "Relative importance (RMSE)",
                                                                               labs_fill = "",
                                                                               legend.position = "none",
                                                                               mul = 3)

merged_covariates_importance_SPRF <- plot_merged_covariates_importance_function(merged_covariates_importance_function(plot_data = bind_files,
                                                                                                                      fitted_model = "sprf"),
                                                                                color = pal_contribution,
                                                                                labs_y = "Relative importance (RMSE)",
                                                                                labs_fill = "",
                                                                                legend.position = c(0.8, 0.18),
                                                                                mul = 3)

merged_covariates_importance <- (merged_covariates_importance_GLM + merged_covariates_importance_GAM) / (merged_covariates_importance_SPAMM + merged_covariates_importance_RF) / (merged_covariates_importance_GBM + merged_covariates_importance_SPRF)

ggsave("figures/merged_covariates_importance.pdf", merged_covariates_importance, height = 15, width = 11)



#################### Per realm contribution ####################

glm_realm <- load_realm_cont_function(files_path = "outputs/glm_biomass_contribution_realm")
rf_realm <- load_realm_cont_function(files_path = "outputs/rf_biomass_contribution_realm")
gam_realm <- load_realm_cont_function(files_path = "outputs/gam_biomass_contribution_realm")
gbm_realm <- load_realm_cont_function(files_path = "outputs/gbm_biomass_contribution_realm")
spamm_realm <- load_realm_cont_function(files_path = "outputs/spamm_biomass_contribution_realm")
sprf_realm <- load_realm_cont_function(files_path = "outputs/sprf_biomass_contribution_realm")
sprf_realm <- do.call(rbind, sprf_realm)
colnames(sprf_realm)[colnames(sprf_realm) %in% c("global_dropout_loss", "global_sd_dropout_loss")] <- c("Dropout_loss", "sd_dropout_loss")
sprf_realm <- sprf_realm |> 
  dplyr::group_split(realm)

contribution_realm_data <- list(glm_realm, gam_realm, spamm_realm, rf_realm, gbm_realm, sprf_realm)
contribution_realm_data <- list(do.call(rbind, lapply(contribution_realm_data, '[[', 1)),
                                do.call(rbind, lapply(contribution_realm_data, '[[', 2)),
                                do.call(rbind, lapply(contribution_realm_data, '[[', 3)),
                                do.call(rbind, lapply(contribution_realm_data, '[[', 4)),
                                do.call(rbind, lapply(contribution_realm_data, '[[', 5)))

merged_covariates_importance_all <- lapply(1:length(contribution_realm_data), function(i) {
  
  merged_covariates_importance_all_function(plot_data = contribution_realm_data[[i]],
                                            title = stringr::str_replace_all(unique(contribution_realm_data[[i]]$realm), c("_" = " ", "-" = " ")),
                                            legend.position = "none",
                                            title.size = 13,
                                            axis.text.x = 11,
                                            axis.text.y = 11,
                                            axis.title = 13,
                                            legend.text = 15,
                                            strip.text.x = 9,
                                            strip.text.y = 9,
                                            geom.text.size = 4,
                                            fill = "")

})

covariates_importance_all <- lapply(1:length(contribution_realm_data), function(i) {
  
  covariates_importance_all_function(plot_data = contribution_realm_data[[i]],
                                     title = stringr::str_replace_all(unique(contribution_realm_data[[i]]$realm), c("_" = " ", "-" = " ")),
                                     legend.position = "none",
                                     title.size = 13,
                                     axis.text.x = 11,
                                     axis.text.y = 11,
                                     axis.title = 13,
                                     legend.text = 15,
                                     strip.text.x = 9,
                                     strip.text.y = 9,
                                     fill = "Covariate categories")
  
})

world <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")

rls_map <- ggplot(data = world) +
  geom_sf(fill = "white", color = "grey", size = 0.03) +
  scale_x_continuous(limits = c(-180, 180)) +
  scale_y_continuous(limits = c(-60, 60)) +
  theme_classic() +
  coord_sf(expand = FALSE) +
  theme(legend.position = "none",
        axis.text = element_text(size = 15))

map_contribution_merged_realm <- rls_map + 
  patchwork::inset_element(merged_covariates_importance_all[[3]], left = 0.27, bottom = 0.53, right = 0.48, top = 1) +
  patchwork::inset_element(merged_covariates_importance_all[[4]], left = 0.30, bottom = 0.07, right = 0.53, top = 0.52) +
  patchwork::inset_element(merged_covariates_importance_all[[1]], left = 0.79, bottom = 0.26, right = 1, top = 0.74) +
  patchwork::inset_element(merged_covariates_importance_all[[2]], left = 0.60, bottom = 0.003, right = 0.81, top = 0.41) +
  patchwork::inset_element(merged_covariates_importance_all[[5]], left = 0.05, bottom = 0.14, right = 0.26, top = 0.61) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggplot2::ggsave("figures/map_contribution_merged_realm.pdf", map_contribution_merged_realm, height = 6, width = 15)

covariates_importance_all_bind <- (covariates_importance_all[[1]] + covariates_importance_all[[2]] + covariates_importance_all[[3]]) / (covariates_importance_all[[4]] + covariates_importance_all[[5]])

ggplot2::ggsave("figures/contribution_realm.png", covariates_importance_all_bind, height = 10, width = 20)
