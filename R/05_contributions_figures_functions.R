# function for evaluating covariates importance

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
    dplyr::inner_join(plot_data)
  
  only_model_best <- only_model_best[only_model_best$best_model == only_model_best$fitted_model,]
  
  # only_model_best$contributions_and_sd[[]]
  
  ENV <- lapply(1:nrow(only_model_best), function(i) { only_model_best$contributions_and_sd[[i]][only_model_best$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$Dropout_loss})
  ENV_sd <- lapply(1:nrow(only_model_best), function(i) { only_model_best$contributions_and_sd[[i]][only_model_best$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$sd_dropout_loss})
  ENV <- do.call(rbind, ENV)
  ENV_sd <- do.call(rbind, ENV_sd)
  ENV <- dplyr::tibble(value = matrixStats::colMedians(ENV),
                       sd = matrixStats::colMedians(ENV_sd),
                       var = c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),
                       VAR = rep("ENV", 8))
  # ENV[ENV$var == "mean_1year_so_mean",3]$var <- "mean_1year_sss"
  # ENV[ENV$var == "min_1year_analysed_sst",3]$var <- "min_1year_sst"
  # ENV[ENV$var == "mean_7days_analysed_sst",3]$var <- "mean_7days_sst"
  # ENV[ENV$var == "max_1year_analysed_sst",3]$var <- "max_1year_sst"
  # ENV[ENV$var == "max_5year_degree_heating_week",3]$var <- "max_5year_dhw"
  
  ENV[ENV$var == "mean_1year_so_mean",3]$var <- "Sea Surface Salinity (mean 1 year)"
  ENV[ENV$var == "min_1year_analysed_sst",3]$var <- "Sea Surface Temperature (min 1 year)"
  ENV[ENV$var == "mean_7days_analysed_sst",3]$var <- "Sea Surface Temperature (mean 7 days) "
  ENV[ENV$var == "max_1year_analysed_sst",3]$var <- "Sea Surface Temperature (max 1 year)"
  ENV[ENV$var == "max_5year_degree_heating_week",3]$var <- "Degree Heating Weeks (max 5 year)"
  ENV[ENV$var == "mean_1year_chl",3]$var <- "Chlorophyll (mean 1 year)"
  ENV[ENV$var == "mean_7days_chl",3]$var <- "Chlorophyll (mean 7 days)"
  ENV[ENV$var == "min_5year_ph",3]$var <- "pH (min 5 year)"
  
  SOC <- lapply(1:nrow(only_model_best), function(i) { only_model_best$contributions_and_sd[[i]][only_model_best$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$Dropout_loss})
  SOC_sd <- lapply(1:nrow(only_model_best), function(i) { only_model_best$contributions_and_sd[[i]][only_model_best$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$sd_dropout_loss})
  SOC <- do.call(rbind, SOC)
  SOC_sd <- do.call(rbind, SOC_sd)
  SOC <- dplyr::tibble(value = matrixStats::colMedians(SOC),
                       sd = matrixStats::colMedians(SOC_sd),
                       var = c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),
                       VAR = rep("HUM", 8))
  # SOC[SOC$var == "n_fishing_vessels",3]$var <- "fishing_vessels_density"
  # SOC[SOC$var == "gravtot2",3]$var <- "gravity"
  
  SOC[SOC$var == "n_fishing_vessels",3]$var <- "Fishing vessels density"
  SOC[SOC$var == "gravtot2",3]$var <- "Gravity"
  SOC[SOC$var == "hdi",3]$var <- "Human development index"
  SOC[SOC$var == "natural_ressource_rent",3]$var <- "Natural Ressource Rent"
  SOC[SOC$var == "neartt",3]$var <- "Neartt"
  SOC[SOC$var == "ngo",3]$var <- "Non governmental organization"
  SOC[SOC$var == "effectiveness",3]$var <- "MPA effectiveness"
  SOC[SOC$var == "gdp",3]$var <- "Growth development product"
  
  
  
  HAB <- lapply(1:nrow(only_model_best), function(i) { only_model_best$contributions_and_sd[[i]][only_model_best$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$Dropout_loss})
  HAB_sd <- lapply(1:nrow(only_model_best), function(i) { only_model_best$contributions_and_sd[[i]][only_model_best$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$sd_dropout_loss})
  HAB <- do.call(rbind, HAB)
  HAB_sd <- do.call(rbind, HAB_sd)
  HAB <- dplyr::tibble(value = matrixStats::colMedians(HAB),
                       sd = matrixStats::colMedians(HAB_sd),
                       var = c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),
                       VAR = rep("HAB", 8))
  # HAB[HAB$var == "Rock_500m",3]$var <- "%_of_rock"
  # HAB[HAB$var == "Sand_500m",3]$var <- "%_of_sand"
  # HAB[HAB$var == "coral_algae_500m",3]$var <- "%_of_coral/algae"
  # HAB[HAB$var == "Rubble_500m",3]$var <- "%_of_rubble"
  # HAB[HAB$var == "coral",3]$var <- "RLS_coral"
  # HAB[HAB$var == "coralline_algae",3]$var <- "RLS_coralline_algae"
  
  HAB[HAB$var == "Rock_500m",3]$var <- "Rock (%)"
  HAB[HAB$var == "Sand_500m",3]$var <- "Sand (%)"
  HAB[HAB$var == "coral_algae_500m",3]$var <- "Coral/Algae (%)"
  HAB[HAB$var == "Rubble_500m",3]$var <- "Rubble (%)"
  HAB[HAB$var == "coral",3]$var <- "Coral (RLS)"
  HAB[HAB$var == "coralline_algae",3]$var <- "Coralline algae (RLS)"
  HAB[HAB$var == "depth",3]$var <- "Depth"
  HAB[HAB$var == "reef_extent",3]$var <- "Reef extent"
  
  cont <- ENV |>  
    dplyr::full_join(HAB) |>
    dplyr::full_join(SOC)
  
  importance_plot <- ggplot(cont) +
    geom_col(aes(x = reorder(var, value), y = value, fill = VAR)) +
    geom_errorbar(aes(x = var, y = value, ymin = value - sd, ymax = value + sd), width = .2,
                  position = position_dodge(.9)) +
    scale_fill_manual(values = c("ENV" = pal_contribution[1],
                                 "HUM" = pal_contribution[3],
                                 "HAB" = pal_contribution[2])) +
    theme_minimal() +
    coord_flip() +
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


var_max_all_function <- function(plot_data)
  
  {
  
  only_model_best <- best_models |>
    dplyr::inner_join(plot_data)
  
  only_model_best <- only_model_best[only_model_best$best_model == only_model_best$fitted_model,]
  
  ENV <- lapply(1:nrow(only_model_best), function(i) { only_model_best$contributions_and_sd[[i]][only_model_best$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$Dropout_loss})
  which_true <- which(lapply(1:length(ENV), function(i) { is.null(ENV[[i]])}) == TRUE)
  ENV[which_true] <- lapply(1:length(ENV[which_true]), function(i) {rep(0, 8)})
  ENV_sd <- lapply(1:nrow(only_model_best), function(i) { only_model_best$contributions_and_sd[[i]][only_model_best$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$sd_dropout_loss})
  ENV <- do.call(rbind, ENV)
  which_true <- which(lapply(1:length(ENV_sd), function(i) { is.null(ENV_sd[[i]])}) == TRUE)
  ENV_sd[which_true] <- lapply(1:length(ENV_sd[which_true]), function(i) {rep(0, 8)})
  ENV_sd <- do.call(rbind, ENV_sd)
  ENV <- dplyr::tibble(species_name = only_model_best$species_name,
                       value = rowMeans(ENV),
                       sd = rowMeans(ENV_sd),
                       var = rep("ENV",nrow(ENV)))
  
  SOC <- lapply(1:nrow(only_model_best), function(i) { only_model_best$contributions_and_sd[[i]][only_model_best$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$Dropout_loss})
  which_true <- which(lapply(1:length(SOC), function(i) { is.null(SOC[[i]])}) == TRUE)
  SOC[which_true] <- lapply(1:length(SOC[which_true]), function(i) {rep(0, 8)})
  SOC_sd <- lapply(1:nrow(only_model_best), function(i) { only_model_best$contributions_and_sd[[i]][only_model_best$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$sd_dropout_loss})
  SOC <- do.call(rbind, SOC)
  which_true <- which(lapply(1:length(SOC_sd), function(i) { is.null(SOC_sd[[i]])}) == TRUE)
  SOC_sd[which_true] <- lapply(1:length(SOC_sd[which_true]), function(i) {rep(0, 8)})
  SOC_sd <- do.call(rbind, SOC_sd)
  SOC <- dplyr::tibble(species_name = only_model_best$species_name,
                       value = rowMeans(SOC),
                       sd = rowMeans(SOC_sd),
                       var = rep("HUM",nrow(SOC)))
  
  HAB <- lapply(1:nrow(only_model_best), function(i) { only_model_best$contributions_and_sd[[i]][only_model_best$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$Dropout_loss})
  which_true <- which(lapply(1:length(HAB), function(i) { is.null(HAB[[i]])}) == TRUE)
  HAB[which_true] <- lapply(1:length(HAB[which_true]), function(i) {rep(0, 8)})
  HAB_sd <- lapply(1:nrow(only_model_best), function(i) { only_model_best$contributions_and_sd[[i]][only_model_best$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$sd_dropout_loss})
  HAB <- do.call(rbind, HAB)
  which_true <- which(lapply(1:length(HAB_sd), function(i) { is.null(HAB_sd[[i]])}) == TRUE)
  HAB_sd[which_true] <- lapply(1:length(HAB_sd[which_true]), function(i) {rep(0, 8)})
  HAB_sd <- do.call(rbind, HAB_sd)
  HAB <- dplyr::tibble(species_name = only_model_best$species_name,
                       value = rowMeans(HAB),
                       sd = rowMeans(HAB_sd),
                       var = rep("HAB",nrow(HAB)))
  
  cont <- ENV |>
    dplyr::full_join(HAB) |>
    dplyr::full_join(SOC)
  
  best <- dplyr::tibble(species_name = unique(cont$species_name),
                        ENV = cont[cont$var == "ENV",2],
                        SOC = cont[cont$var == "HUM",2],
                        HAB = cont[cont$var == "HAB",2])
  best <- dplyr::inner_join(best, best_models, by = "species_name")
  best <- as.matrix(best)
  
  best <- lapply(1:nrow(best), function(i) {
    test <- best[i,]
    testt <- which.max(test[2:4])
    testtt <- dplyr::tibble(species_name = test[1],
                            varmax = testt)
  })
  
  best <- do.call(rbind, best)
  best$varmax <- as.character(best$varmax)
  best[best$varmax == "1" ,2] <- "ENV"
  best[best$varmax == "2" ,2] <- "HUM"
  best[best$varmax == "3" ,2] <- "HAB"
  
  cont <- cont |>
    dplyr::inner_join(best, by = c("species_name"))
  
  var_max <- cont |>
    dplyr::group_by(varmax) |>
    dplyr::summarise(n = dplyr::n()/3)
  var_max <- var_max |>
    dplyr::rename(VAR = varmax)
  
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
                                                      fill){
  
  require(ggplot2)
  
  # covariates relative importance by median
  
  
  only_model_best_best <- best_models |>
    dplyr::inner_join(plot_data)
  
  only_model_best_best <- only_model_best_best[only_model_best_best$best_model == only_model_best_best$fitted_model,]
  
  ENV <- lapply(1:nrow(only_model_best_best), function(i) { only_model_best_best$contributions_and_sd[[i]][only_model_best_best$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$Dropout_loss})
  ENV_sd <- lapply(1:nrow(only_model_best_best), function(i) { only_model_best_best$contributions_and_sd[[i]][only_model_best_best$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$sd_dropout_loss})
  ENV <- do.call(rbind, ENV)
  ENV_sd <- do.call(rbind, ENV_sd)
  ENV <- dplyr::tibble(value = matrixStats::colMedians(ENV),
                       sd = matrixStats::colMedians(ENV_sd),
                       var = c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),
                       VAR = rep("ENV", 8))
  
  SOC <- lapply(1:nrow(only_model_best_best), function(i) { only_model_best_best$contributions_and_sd[[i]][only_model_best_best$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$Dropout_loss})
  SOC_sd <- lapply(1:nrow(only_model_best_best), function(i) { only_model_best_best$contributions_and_sd[[i]][only_model_best_best$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$sd_dropout_loss})
  SOC <- do.call(rbind, SOC)
  SOC_sd <- do.call(rbind, SOC_sd)
  SOC <- dplyr::tibble(value = matrixStats::colMedians(SOC),
                       sd = matrixStats::colMedians(SOC_sd),
                       var = c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),
                       VAR = rep("HUM", 8))
  
  HAB <- lapply(1:nrow(only_model_best_best), function(i) { only_model_best_best$contributions_and_sd[[i]][only_model_best_best$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$Dropout_loss})
  HAB_sd <- lapply(1:nrow(only_model_best_best), function(i) { only_model_best_best$contributions_and_sd[[i]][only_model_best_best$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$sd_dropout_loss})
  HAB <- do.call(rbind, HAB)
  HAB_sd <- do.call(rbind, HAB_sd)
  HAB <- dplyr::tibble(value = matrixStats::colMedians(HAB),
                       sd = matrixStats::colMedians(HAB_sd),
                       var = c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),
                       VAR = rep("HAB", 8))
  
  cont <- ENV |>  
    dplyr::full_join(HAB) |>
    dplyr::full_join(SOC)
  
  var_max <- var_max_all_function(plot_data = plot_data)
  
  #merge contribution per var and model
  
  cont_merge <- cont |>  
    dplyr::group_by(VAR) |>
    dplyr::summarise(value = mean(value),
                     sd = mean(sd))
  
  cont_merge <- dplyr::inner_join(cont_merge, var_max, by = c("VAR"))
  
  merged_importance_plot <- ggplot(cont_merge) +
    geom_col(aes(x = reorder(VAR, value), y = value, fill = VAR)) +
    geom_errorbar(aes(x = VAR, y = value, ymin=value-sd, ymax=value+sd), width=.1,
                  position=position_dodge(.9)) +
    geom_text(aes(x = VAR, y = value+3*sd, label = n), size = geom.text.size) +
    scale_fill_manual(values = c("ENV" = pal_contribution[1],
                                 "HUM" = pal_contribution[3],
                                 "HAB" = pal_contribution[2])) +
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
          strip.text.y = element_text(size = strip.text.y)) +
    scale_y_continuous(breaks = c(0, round(mean(cont_merge$value), 2), round(max(cont_merge$value) + (max(cont_merge$value) * 0.2), 2)))
    
          
}

plot_data = bind_files
fitted_model = "glm"
color = pal_contribution
labs_y = ""
labs_fill = ""
ylim = c(0,0.2)
legend.position = "none"



covariates_importance_function <- function(plot_data,
                                           fitted_model,
                                           color,
                                           labs_y,
                                           labs_fill,
                                           ylim,
                                           legend.position
                                           ){
  
  require(ggplot2)

  # covariates relative importance by mean

  plot_level <- fitted_model
  plot_data <- plot_data |> 
    dplyr::filter(fitted_model == plot_level)
  plot_data <- plot_data |>
    dplyr::filter(species_name %in% best_models[best_models$best_model == plot_level,1]$species_name)

  ENV <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$Dropout_loss})
  ENV_sd <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$sd_dropout_loss})
  ENV <- do.call(rbind, ENV)
  ENV_sd <- do.call(rbind, ENV_sd)
  ENV <- dplyr::tibble(value = matrixStats::colMedians(ENV),
                       sd = matrixStats::colMedians(ENV_sd),
                       var = c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),
                       VAR = rep("ENV", 8),
                       plot_level = rep(plot_level, 8))
    
  SOC <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$Dropout_loss})
  SOC_sd <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$sd_dropout_loss})
  SOC <- do.call(rbind, SOC)
  SOC_sd <- do.call(rbind, SOC_sd)
  SOC <- dplyr::tibble(value = matrixStats::colMedians(SOC),
                       sd = matrixStats::colMedians(SOC_sd),
                       var = c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),
                       VAR = rep("HUM", 8),
                       plot_level = rep(plot_level, 8))
  SOC[SOC$var == "n_fishing_vessels",3]$var <- "fishing_vessels_density"
    
  HAB <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$Dropout_loss})
  HAB_sd <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$sd_dropout_loss})
  HAB <- do.call(rbind, HAB)
  HAB_sd <- do.call(rbind, HAB_sd)
  HAB <- dplyr::tibble(value = matrixStats::colMedians(HAB),
                       sd = matrixStats::colMedians(HAB_sd),
                       var = c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),
                       VAR = rep("HAB", 8),
                       plot_level = rep(plot_level, 8))
    
  cont <- ENV |>  
    dplyr::full_join(HAB) |> 
    dplyr::full_join(SOC)

  importance_plot <- ggplot(cont) +
    geom_col(aes(x = reorder(var, value), y = value, fill = VAR)) +
    geom_errorbar(aes(x = var, y = value, ymin = value - sd, ymax = value + sd), width = .2,
                    position = position_dodge(.9)) +
    scale_fill_manual(values = c("ENV" = color[1],
                                 "HUM" = color [3],
                                 "HAB" = color [2])) +
    theme_minimal() +
    coord_flip(ylim = ylim) +
    facet_grid(~ plot_level) +
    labs(y = labs_y, x = "", fill = labs_fill) +
    theme(legend.position = legend.position) +
    theme(title = element_text(size = 18),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 10),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20))

}


var_max_function <- function(plot_data,
                             fitted_model
                             ){
  
  plot_level <- fitted_model
  plot_data <- plot_data |> 
    dplyr::filter(fitted_model == plot_level)
  plot_data <- plot_data |> 
    dplyr::filter(species_name %in% best_models[best_models$best_model == plot_level,1]$species_name)

  ENV <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$Dropout_loss})
  ENV_sd <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$sd_dropout_loss})
  ENV <- do.call(rbind, ENV)
  ENV_sd <- do.call(rbind, ENV_sd)
  ENV <- dplyr::tibble(species_name = plot_data$species_name,
                       value = rowMeans(ENV),
                       sd = rowMeans(ENV_sd),
                       var = rep("ENV",nrow(ENV)),
                       plot_level = rep(plot_level, nrow(ENV)))
  
  SOC <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$Dropout_loss})
  SOC_sd <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$sd_dropout_loss})
  SOC <- do.call(rbind, SOC)
  SOC_sd <- do.call(rbind, SOC_sd)
  SOC <- dplyr::tibble(species_name = plot_data$species_name,
                       value = rowMeans(SOC),
                       sd = rowMeans(SOC_sd),
                       var = rep("HUM",nrow(SOC)),
                       plot_level = rep(plot_level, nrow(SOC)))
  
  HAB <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$Dropout_loss})
  HAB_sd <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$sd_dropout_loss})
  HAB <- do.call(rbind, HAB)
  HAB_sd <- do.call(rbind, HAB_sd)
  HAB <- dplyr::tibble(species_name = plot_data$species_name,
                       value = rowMeans(HAB),
                       sd = rowMeans(HAB_sd),
                       var = rep("HAB",nrow(HAB)),
                       plot_level = rep(plot_level,nrow(HAB)))
  
  cont <- ENV |> 
    dplyr::full_join(HAB) |> 
    dplyr::full_join(SOC)

  best <- dplyr::tibble(species_name = unique(cont$species_name),
                        ENV = cont[cont$var == "ENV",2],
                        SOC = cont[cont$var == "HUM",2],
                        HAB = cont[cont$var == "HAB",2])
  best <- dplyr::inner_join(best, best_models, by = "species_name")
  best <- as.matrix(best)
  
  best <- lapply(1:nrow(best), function(i) {
    test <- best[i,]
    testt <- which.max(test[2:4])
    testtt <- dplyr::tibble(species_name = test[1],
                            varmax = testt,
                            plot_level = test[5])
  })
  
  best <- do.call(rbind, best)
  best$varmax <- as.character(best$varmax)
  best[best$varmax == "1" ,2] <- "ENV"
  best[best$varmax == "2" ,2] <- "HUM"
  best[best$varmax == "3" ,2] <- "HAB"
  
  cont <- cont |> 
    dplyr::inner_join(best, by = c("species_name", "plot_level"))
  
  var_max <- cont |> 
    dplyr::group_by(plot_level, varmax) |> 
    dplyr::summarise(n = dplyr::n()/3)
  var_max <- var_max |> 
    dplyr::rename(VAR = varmax)
  
}

merged_covariates_importance_function <- function(plot_data,
                                                  fitted_model,
                                                  color,
                                                  labs_y,
                                                  labs_fill,
                                                  legend.position,
                                                  mul
                                                  ){
  
  require(ggplot2)

  # covariates relative importance by median
  
  plot_level <- fitted_model
  plot_data <- plot_data |> 
    dplyr::filter(fitted_model == plot_level)
  plot_data <- plot_data |> 
    dplyr::filter(species_name %in% best_models[best_models$best_model == plot_level,1]$species_name)
  
  ENV <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$Dropout_loss})
  ENV_sd <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$sd_dropout_loss})
  ENV <- do.call(rbind, ENV)
  ENV_sd <- do.call(rbind, ENV_sd)
  ENV <- dplyr::tibble(value = matrixStats::colMedians(ENV),
                       sd = matrixStats::colMedians(ENV_sd),
                       var = c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),
                       VAR = rep("ENV", 8),
                       plot_level = rep(plot_level, 8))
  
  SOC <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$Dropout_loss})
  SOC_sd <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$sd_dropout_loss})
  SOC <- do.call(rbind, SOC)
  SOC_sd <- do.call(rbind, SOC_sd)
  SOC <- dplyr::tibble(value = matrixStats::colMedians(SOC),
                       sd = matrixStats::colMedians(SOC_sd),
                       var = c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),
                       VAR = rep("HUM", 8),
                       plot_level = rep(plot_level, 8))
  
  HAB <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$Dropout_loss})
  HAB_sd <- lapply(1:nrow(plot_data), function(i) { plot_data$contributions_and_sd[[i]][plot_data$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$sd_dropout_loss})
  HAB <- do.call(rbind, HAB)
  HAB_sd <- do.call(rbind, HAB_sd)
  HAB <- dplyr::tibble(value = matrixStats::colMedians(HAB),
                       sd = matrixStats::colMedians(HAB_sd),
                       var = c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),
                       VAR = rep("HAB", 8),
                       plot_level = rep(plot_level, 8))
  
  cont <- ENV |>  
    dplyr::full_join(HAB) |> 
    dplyr::full_join(SOC)

  var_max <- var_max_function(plot_data = bind_files,
                              fitted_model = plot_level)
  
  #merge contribution per var and model
  
  cont_merge <- cont |>  
    dplyr::group_by(VAR, plot_level) |> 
    dplyr::summarise(value = mean(value),
                     sd = mean(sd))
  
  cont_merge <- dplyr::inner_join(cont_merge, var_max, by = c("plot_level", "VAR"))
  
  merged_importance_plot <- ggplot(cont_merge) +
    geom_col(aes(x = reorder(VAR, value), y = value, fill = VAR)) +
    geom_errorbar(aes(x = VAR, y = value, ymin=value-sd, ymax=value+sd), width=.1,
                  position=position_dodge(.9)) +
    geom_text(aes(x = VAR, y = value+mul*sd, label = n), size = 5) +
    scale_fill_manual(values = c("ENV" = color[1],
                                 "HUM" = color[3],
                                 "HAB" = color[2])) +
    theme_minimal() +
    coord_flip() +
    facet_grid(~plot_level) +
    labs(y = labs_y, x = "", fill = labs_fill) +
    theme(legend.position = legend.position) +
    theme(title = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 10),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20))
  
}



