# load in packages ---- 

libs <- c('tidyverse', 'gridExtra', 'ggplot2', 'patchwork', 'matrixStats', 'parallel','PNWColors')
lapply(libs, library, character.only = T, lib.loc = .libPaths()[1])

# check all packages are loaded
if(sum(libs %in% (.packages())) != length(libs)){
  stop('packages not loaded correctly')}

# source functions ----

source("scripts-final/00_functions/contributions_figures_functions.R")

pal_contribution <- PNWColors::pnw_palette("Bay",3, type = "discrete")

best_models <- readRDS("results/overall_best_models.rds")

#### Covariates contribution plot ####

bind_files <- list.files('results/contributions_bind/', full.names = T)

##### For bind_files

Contributions_biomass <- readRDS(bind_files)

covariates_importance_GLM <- covariates_importance_function(plot_data = Contributions_biomass,
                                                           fitted_model = "GLM",
                                                           color = pal_contribution,
                                                           labs_y = "",
                                                           labs_fill = "",
                                                           ylim = c(0,0.36),
                                                           legend.position = "none")

covariates_importance_GAM <- covariates_importance_function(plot_data = Contributions_biomass,
                                                           fitted_model = "GAM",
                                                           color = pal_contribution,
                                                           labs_y = "",
                                                           labs_fill = "",
                                                           ylim = c(0,0.21),
                                                           legend.position = "none")

covariates_importance_SPAMM <- covariates_importance_function(plot_data = Contributions_biomass,
                                                           fitted_model = "SPAMM",
                                                           color = pal_contribution,
                                                           labs_y = "",
                                                           labs_fill = "",
                                                           ylim = c(0,0.21),
                                                           legend.position = "none")

covariates_importance_RF <- covariates_importance_function(plot_data = Contributions_biomass,
                                                           fitted_model = "RF",
                                                           color = pal_contribution,
                                                           labs_y = "",
                                                           labs_fill = "",
                                                           ylim = c(0,0.21),
                                                           legend.position = "none")

covariates_importance_GBM <- covariates_importance_function(plot_data = Contributions_biomass,
                                                           fitted_model = "GBM",
                                                           color = pal_contribution,
                                                           labs_y = "Change in RMSE",
                                                           labs_fill = "",
                                                           ylim = c(0,0.21),
                                                           legend.position = "none")

covariates_importance_SPRF <- covariates_importance_function(plot_data = Contributions_biomass,
                                                             fitted_model = "SPRF",
                                                             color = pal_contribution,
                                                             labs_y = "Change in RMSE",
                                                             labs_fill = "",
                                                             ylim = c(0,0.21),
                                                             legend.position = c(0.65, 0.2))

covariates_importance_all <- (covariates_importance_GLM + covariates_importance_GAM) / (covariates_importance_SPAMM + covariates_importance_RF) / (covariates_importance_GBM + covariates_importance_SPRF)

ggsave("figures-R3/covariates_importance_all.pdf", covariates_importance_all, height = 15, width = 11)

merged_covariates_importance_GLM <- merged_covariates_importance_function(plot_data = Contributions_biomass,
                                                                          fitted_model = "GLM",
                                                                          color = pal_contribution,
                                                                          labs_y = "",
                                                                          labs_fill = "",
                                                                          legend.position = "none",
                                                                          mul = 2)

merged_covariates_importance_GAM <- merged_covariates_importance_function(plot_data = Contributions_biomass,
                                                                          fitted_model = "GAM",
                                                                          color = pal_contribution,
                                                                          labs_y = "",
                                                                          labs_fill = "",
                                                                          legend.position = "none",
                                                                          mul = 2)

merged_covariates_importance_SPAMM <- merged_covariates_importance_function(plot_data = Contributions_biomass,
                                                                            fitted_model = "SPAMM",
                                                                            color = pal_contribution,
                                                                            labs_y = "",
                                                                            labs_fill = "",
                                                                            legend.position = "none",
                                                                            mul = 2)

merged_covariates_importance_RF <- merged_covariates_importance_function(plot_data = Contributions_biomass,
                                                                         fitted_model = "RF",
                                                                         color = pal_contribution,
                                                                         labs_y = "",
                                                                         labs_fill = "",
                                                                         legend.position = "none",
                                                                         mul = 3)

merged_covariates_importance_GBM <- merged_covariates_importance_function(plot_data = Contributions_biomass,
                                                                          fitted_model = "GBM",
                                                                          color = pal_contribution,
                                                                          labs_y = "Change in RMSE",
                                                                          labs_fill = "",
                                                                          legend.position = "none",
                                                                          mul = 3)

merged_covariates_importance_SPRF <- merged_covariates_importance_function(plot_data = Contributions_biomass,
                                                                           fitted_model = "SPRF",
                                                                           color = pal_contribution,
                                                                           labs_y = "Change in RMSE",
                                                                           labs_fill = "",
                                                                           legend.position = c(0.8, 0.19),
                                                                           mul = 3)

merged_covariates_importance <- (merged_covariates_importance_GLM + merged_covariates_importance_GAM) / (merged_covariates_importance_SPAMM + merged_covariates_importance_RF) / (merged_covariates_importance_GBM + merged_covariates_importance_SPRF)

ggsave("figures-R3/merged_covariates_importance.pdf", merged_covariates_importance, height = 15, width = 11)

