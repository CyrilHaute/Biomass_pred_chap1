# source functions ----

library(patchwork)

source("R/05_contributions_figures_functions.R")

pal_contribution <- PNWColors::pnw_palette("Bay", 3, type = "discrete")

load("outputs/best_models.Rdata")

#################### Global contribution ####################


#### Covariates contribution plot ####

sprf <- list.files("outputs/biomass_contribution/sprf", full.names = T)

bind_files_sprf <- lapply(1:length(sprf), function(i) {
  
  load(sprf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
bind_files_sprf <- do.call(rbind, bind_files_sprf)
bind_files_sprf$contributions_and_sd <- lapply(1:nrow(bind_files_sprf), function(i) { bind_files_sprf$contributions_and_sd[[i]] |> 
    dplyr::rename(Dropout_loss = global_dropout_loss,
                  sd_dropout_loss = global_sd_dropout_loss)
  })
bind_files_sprf <- list(bind_files_sprf)

list_files_path <- list.files("outputs/biomass_contribution", full.names = T)
list_files_path <- list_files_path[which(grepl(pattern = "sprf", list_files_path) == FALSE)]
bind_files <- lapply(1:length(list_files_path), function(i) {
  
  load(list_files_path[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
bind_files <- c(bind_files, bind_files_sprf)
bind_files <- do.call(rbind, bind_files)
models <- unique(bind_files$fitted_model)

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
                                                                strip.text.y = 20)

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
                                                                              geom.text.size = 5)

covariates_importance_all_and_merged <- covariates_importance_all + merged_covariates_importance_all

ggsave("figures/covariates_importance_all_and_merged.pdf", covariates_importance_all_and_merged, height = 7, width = 14)
ggsave("figures/covariates_importance_all_and_merged.png", covariates_importance_all_and_merged, height = 7, width = 14)

covariates_importance_GLM <- covariates_importance_function(plot_data = bind_files,
                                                            fitted_model = "glm",
                                                            color = pal_contribution,
                                                            labs_y = "",
                                                            labs_fill = "",
                                                            ylim = c(0,1),
                                                            legend.position = "none")

covariates_importance_GAM <- covariates_importance_function(plot_data = bind_files,
                                                            fitted_model = "gam",
                                                            color = pal_contribution,
                                                            labs_y = "",
                                                            labs_fill = "",
                                                            ylim = c(0,0.7),
                                                            legend.position = "none")

covariates_importance_SPAMM <- covariates_importance_function(plot_data = bind_files,
                                                              fitted_model = "spamm",
                                                              color = pal_contribution,
                                                              labs_y = "",
                                                              labs_fill = "",
                                                              ylim = c(0,0.8),
                                                              legend.position = "none")

covariates_importance_RF <- covariates_importance_function(plot_data = bind_files,
                                                           fitted_model = "rf",
                                                           color = pal_contribution,
                                                           labs_y = "",
                                                           labs_fill = "",
                                                           ylim = c(0,0.15),
                                                           legend.position = "none")

covariates_importance_GBM <- covariates_importance_function(plot_data = bind_files,
                                                            fitted_model = "gbm",
                                                            color = pal_contribution,
                                                            labs_y = "Relative importance (RMSE)",
                                                            labs_fill = "",
                                                            ylim = c(0,0.15),
                                                            legend.position = "none")

covariates_importance_SPRF <- covariates_importance_function(plot_data = bind_files,
                                                             fitted_model = "sprf",
                                                             color = pal_contribution,
                                                             labs_y = "Relative importance (RMSE)",
                                                             labs_fill = "",
                                                             ylim = c(0,0.15),
                                                             legend.position = c(0.75, 0.16))

covariates_importance_all <- (covariates_importance_GLM + covariates_importance_GAM) / (covariates_importance_SPAMM + covariates_importance_RF) / (covariates_importance_GBM + covariates_importance_SPRF)

ggsave("figures/covariates_importance_all.pdf", covariates_importance_all, height = 15, width = 11)
ggsave("figures/covariates_importance_all.png", covariates_importance_all, height = 15, width = 11)

merged_covariates_importance_GLM <- merged_covariates_importance_function(plot_data = bind_files,
                                                                          fitted_model = "glm",
                                                                          color = pal_contribution,
                                                                          labs_y = "",
                                                                          labs_fill = "",
                                                                          legend.position = "none",
                                                                          mul = 2)

merged_covariates_importance_GAM <- merged_covariates_importance_function(plot_data = bind_files,
                                                                          fitted_model = "gam",
                                                                          color = pal_contribution,
                                                                          labs_y = "",
                                                                          labs_fill = "",
                                                                          legend.position = "none",
                                                                          mul = 2)

merged_covariates_importance_SPAMM <- merged_covariates_importance_function(plot_data = bind_files,
                                                                            fitted_model = "spamm",
                                                                            color = pal_contribution,
                                                                            labs_y = "",
                                                                            labs_fill = "",
                                                                            legend.position = "none",
                                                                            mul = 2)

merged_covariates_importance_RF <- merged_covariates_importance_function(plot_data = bind_files,
                                                                         fitted_model = "rf",
                                                                         color = pal_contribution,
                                                                         labs_y = "",
                                                                         labs_fill = "",
                                                                         legend.position = "none",
                                                                         mul = 3)

merged_covariates_importance_GBM <- merged_covariates_importance_function(plot_data = bind_files,
                                                                          fitted_model = "gbm",
                                                                          color = pal_contribution,
                                                                          labs_y = "Relative importance (RMSE)",
                                                                          labs_fill = "",
                                                                          legend.position = "none",
                                                                          mul = 3)

merged_covariates_importance_SPRF <- merged_covariates_importance_function(plot_data = bind_files,
                                                                           fitted_model = "sprf",
                                                                           color = pal_contribution,
                                                                           labs_y = "Relative importance (RMSE)",
                                                                           labs_fill = "",
                                                                           legend.position = c(0.8, 0.18),
                                                                           mul = 3)

merged_covariates_importance <- (merged_covariates_importance_GLM + merged_covariates_importance_GAM) / (merged_covariates_importance_SPAMM + merged_covariates_importance_RF) / (merged_covariates_importance_GBM + merged_covariates_importance_SPRF)

ggsave("figures/merged_covariates_importance.pdf", merged_covariates_importance, height = 15, width = 11)
ggsave("figures/merged_covariates_importance.png", merged_covariates_importance, height = 15, width = 11)



#################### Per realm contribution ####################


realm_full <- list.files("outputs/biomass_contribution_realm", full.names = T)
realm_small <- list.files("outputs/biomass_contribution_realm", full.names = F)

load_outputs <- lapply(1:length(realm_full), function(i) {
  
  realms_all <- list.files(realm_full[i], full.names = T)
  realms <- realms_all[which(grepl(pattern = "sprf", realms_all) == FALSE)]
  
  realm_j <- lapply(1:length(realms), function(j) {
    
    load(realms[j])
    assign(paste0("model_", j), extracted_contributions)
    
  })
  
  sprf <- list.files(realms_all[which(grepl(pattern = "sprf", realms_all) == TRUE)], full.names = T)

  load_outputs_sprf <- lapply(1:length(sprf), function(i) {

    load(sprf[i])
    assign(paste0("model_", i), extracted_contributions)

  })
  load_outputs_sprf <- do.call(rbind, load_outputs_sprf)
  load_outputs_sprf$contributions_and_sd <- lapply(1:nrow(load_outputs_sprf), function(i) { load_outputs_sprf$contributions_and_sd[[i]] |> 
      dplyr::rename(Dropout_loss = global_dropout_loss,
                    sd_dropout_loss = global_sd_dropout_loss)
  })
  load_outputs_sprf <- list(load_outputs_sprf)

  load_outputs <- c(realm_j, load_outputs_sprf)
  
})

names(load_outputs) <- realm_small

merged_covariates_importance_all <- lapply(1:length(load_outputs), function(i) {
  
  merged_covariates_importance_all_function(plot_data = do.call(rbind, load_outputs[[i]]),
                                            title = stringr::str_replace_all(names(load_outputs)[i], c("_" = " ", "-" = " ")),
                                            legend.position = "none",
                                            title.size = 12,
                                            axis.text.x = 11,
                                            axis.text.y = 0,
                                            axis.title = 12,
                                            legend.text = 9,
                                            strip.text.x = 9,
                                            strip.text.y = 9,
                                            geom.text.size = 4,
                                            fill = "Covariate categories")

})

covariates_importance_all <- lapply(1:length(load_outputs), function(i) {
  
  covariates_importance_all_function(plot_data = do.call(rbind, load_outputs[[i]]),
                                     title = stringr::str_replace_all(names(load_outputs)[i], c("_" = " ", "-" = " ")),
                                     legend.position = "none",
                                     title.size = 12,
                                     axis.text.x = 8,
                                     axis.text.y = 8,
                                     axis.title = 13,
                                     legend.text = 9,
                                     strip.text.x = 9,
                                     strip.text.y = 9,
                                     fill = "Covariate categories")
  
})

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

rls_map <- ggplot(data = world) +
  geom_sf() +
  scale_x_continuous(limits=c(-180, 180)) +
  scale_y_continuous(limits=c(-60, 60)) +
  geom_sf(data = sf::st_union(world), fill = "white", color = "gray90", size = 0.01) +
  theme_classic() +
  coord_sf(expand = FALSE) +
  theme(legend.position = "none",
        axis.text = element_text(size = 15))

map_contribution_merged_realm <- rls_map + 
  patchwork::inset_element(merged_covariates_importance_all[[3]], left = 0.29, bottom = 0.62, right = 0.48, top = 0.92) +
  patchwork::inset_element(merged_covariates_importance_all[[4]], left = 0.34, bottom = 0.30, right = 0.53, top = 0.60) +
  patchwork::inset_element(merged_covariates_importance_all[[1]], left = 0.81, bottom = 0.43, right = 1, top = 0.74) +
  patchwork::inset_element(merged_covariates_importance_all[[2]], left = 0.62, bottom = 0.20, right = 0.81, top = 0.50) +
  patchwork::inset_element(merged_covariates_importance_all[[5]], left = 0.07, bottom = 0.31, right = 0.26, top = 0.61) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")


ggplot2::ggsave("figures/map_contribution_merged_realm.pdf", map_contribution_merged_realm, height = 10, width = 20)
ggplot2::ggsave("figures/map_contribution_merged_realm.png", map_contribution_merged_realm, height = 10, width = 20)

covariates_importance_all_bind <- (covariates_importance_all[[1]] + covariates_importance_all[[2]] + covariates_importance_all[[3]]) / (covariates_importance_all[[4]] + covariates_importance_all[[5]])

# ggplot2::ggsave("figures/contribution_realm.pdf", covariates_importance_all_bind, height = 10, width = 20)
# ggplot2::ggsave("figures/contribution_realm.png", covariates_importance_all_bind, height = 10, width = 20)
