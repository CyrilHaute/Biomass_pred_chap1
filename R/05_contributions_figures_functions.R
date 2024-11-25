# function for evaluating covariates importance

# plot_data = bind_files
# title = "A."
# legend.position = "none"
# title.size = 18
# axis.text.x = 15
# axis.text.y = 17
# axis.title = 21
# legend.text = 15
# strip.text.x = 20
# strip.text.y = 20
# fill = ""

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

  ENV <- plot_data |> 
    dplyr::filter(variable %in% env_var) |>
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "ENV")
  ENV[ENV$variable == "max_1year_analysed_sst",]$variable <- "Sea Surface Temperature (max 1 year)"
  ENV[ENV$variable == "min_1year_analysed_sst",]$variable <- "Sea Surface Temperature (min 1 year)"
  ENV[ENV$variable == "max_5year_degree_heating_week",]$variable <- "Degree Heating Week (max 5 year)"
  ENV[ENV$variable == "mean_1year_nppv",]$variable <- "Net Primary Productivity (mean 1 year)"
  ENV[ENV$variable == "mean_1year_so_mean",]$variable <- "Sea Surface Salinity (mean 1 year)"
  ENV[ENV$variable == "min_5year_ph",]$variable <- "pH (min 5 year)"

  SOC <- plot_data |> 
    dplyr::filter(variable %in% hum_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "HUM")
  SOC[SOC$variable == "protection_status2",]$variable <- "MPA status"
  SOC[SOC$variable == "gdp",]$variable <- "Growth Development Product"
  SOC[SOC$variable == "gravtot2",]$variable <- "Human Gravity"
  SOC[SOC$variable == "n_fishing_vessels",]$variable <- "Fishing Vessels Density"
  SOC[SOC$variable == "neartt",]$variable <- "Nearest Population"
  SOC[SOC$variable == "marine_ecosystem_dependency",]$variable <- "Marine Ecosystem Dependency"
  
  HAB <- plot_data |> 
    dplyr::filter(variable %in% hab_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "HAB")
  HAB[HAB$variable == "Rock_500m",]$variable <- "Rock (%)"
  HAB[HAB$variable == "Sand_500m",]$variable <- "Sand (%)"
  HAB[HAB$variable == "coral_algae_500m",]$variable <- "Coral/Algae (%)"
  HAB[HAB$variable == "coral",]$variable <- "Coral (RLS)"
  HAB[HAB$variable == "depth",]$variable <- "Depth"
  HAB[HAB$variable == "reef_extent",]$variable <- "Reef extent"
  
  BIOT <- plot_data |> 
    dplyr::filter(variable %in% biot_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "BIOT")
  BIOT[BIOT$variable == "diversity",]$variable <- "Diversity"
  BIOT[BIOT$variable == "max_trophic",]$variable <- "Maximum Trophic Level"
  BIOT[BIOT$variable == "mean_trophic",]$variable <- "Trophic Level (mean)"
  BIOT[BIOT$variable == "n_trophic",]$variable <- "Number of Throphic Level"
  
  cont <- ENV |>  
    dplyr::full_join(HAB) |> 
    dplyr::full_join(SOC) |> 
    dplyr::full_join(BIOT)

  importance_plot <- ggplot(cont) +
    geom_col(aes(x = reorder(variable, value), y = value, fill = VAR)) +
    geom_errorbar(aes(x = variable, y = value, ymin = value - sd, ymax = value + sd), width = .2,
                  position = position_dodge(.9)) +
    scale_fill_manual(values = c("ENV" = pal_contribution[2],
                                 "HUM" = pal_contribution[1],
                                 "HAB" = pal_contribution[13],
                                 "BIOT" = pal_contribution[3])) +
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
  
  ENV <- only_model_best |> 
    dplyr::filter(variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "min_1year_analysed_sst", "min_5year_ph")) |> 
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = mean(Dropout_loss),
                     sd = mean(sd_dropout_loss),
                     VAR = "ENV")

  SOC <- only_model_best |> 
    dplyr::filter(variable %in% c("effectiveness", "gdp", "gravtot2", "no_violence", "n_fishing_vessels", "natural_ressource_rent", "neartt", "marine_ecosystem_dependency")) |> 
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = mean(Dropout_loss),
                     sd = mean(sd_dropout_loss),
                     VAR = "HUM")

  HAB <- only_model_best |> 
    dplyr::filter(variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "seagrass", "depth", "reef_extent")) |> 
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = mean(Dropout_loss),
                     sd = mean(sd_dropout_loss),
                     VAR = "HAB")
  
  BIOT <- only_model_best |> 
    dplyr::filter(variable %in% c("delta_biomass", "diversity", "max_trophic", "mean_biomass", "mean_trophic", "n_trophic")) |> 
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = mean(Dropout_loss),
                     sd = mean(sd_dropout_loss),
                     VAR = "BIOT")

  cont <- ENV |>
    dplyr::full_join(HAB) |>
    dplyr::full_join(SOC) |> 
    dplyr::full_join(BIOT)
  
  best <- dplyr::tibble(species_name = unique(cont$species_name),
                        ENV = cont[cont$VAR == "ENV",2],
                        SOC = cont[cont$VAR == "HUM",2],
                        HAB = cont[cont$VAR == "HAB",2],
                        BIOT = cont[cont$VAR == "BIOT",2])
  best <- dplyr::inner_join(best, best_models, by = "species_name")
  best <- as.matrix(best)
  
  best <- lapply(1:nrow(best), function(i) {
    test <- best[i,]
    # testt <- which.max(test[2:4])
    testt <- which.max(test[2:5])
    testtt <- dplyr::tibble(species_name = test[1],
                            varmax = testt,
                            # plot_level = test[5],
                            plot_level = test[6])
  })
  
  best <- do.call(rbind, best)
  best$varmax <- as.character(best$varmax)
  best[best$varmax == "1" ,2] <- "ENV"
  best[best$varmax == "2" ,2] <- "HUM"
  best[best$varmax == "3" ,2] <- "HAB"
  best[best$varmax == "4" ,2] <- "BIOT"
  
  cont <- cont |>
    dplyr::inner_join(best, by = c("species_name"))
  
  var_max <- cont |>
    dplyr::group_by(varmax) |>
    dplyr::summarise(n = dplyr::n()/4)
  var_max <- var_max |>
    dplyr::rename(VAR = varmax)
  
}

# plot_data = contribution_realm_data[[i]]
# title = stringr::str_replace_all(unique(contribution_realm_data[[i]]$realm), c("_" = " ", "-" = " "))
# legend.position = "none"
# title.size = 13
# axis.text.x = 11
# axis.text.y = 11
# axis.title = 13
# legend.text = 15
# strip.text.x = 9
# strip.text.y = 9
# geom.text.size = 4
# fill = ""

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
  
  ENV <- only_model_best_best |> 
    dplyr::filter(variable %in% env_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "ENV")
  
  SOC <- only_model_best_best |> 
    dplyr::filter(variable %in% hum_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "HUM")
  
  HAB <- only_model_best_best |> 
    dplyr::filter(variable %in% hab_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "HAB")
  
  BIOT <- only_model_best_best |> 
    dplyr::filter(variable %in% biot_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "BIOT")

  cont <- ENV |>  
    dplyr::full_join(HAB) |> 
    dplyr::full_join(SOC) |> 
    dplyr::full_join(BIOT)

  var_max <- var_max_all_function(plot_data = plot_data)
  
  #merge contribution per var and model
  
  cont_merge <- cont |>  
    dplyr::group_by(VAR) |>
    dplyr::summarise(value = mean(value),
                     sd = mean(sd))
  
  cont_merge <- dplyr::full_join(cont_merge, var_max, by = c("VAR"))
  cont_merge[which(is.na(cont_merge$n)),]$n <- 0
  
  merged_importance_plot <- ggplot(cont_merge) +
    geom_col(aes(x = reorder(VAR, value), y = value, fill = VAR)) +
    geom_errorbar(aes(x = VAR, y = value, ymin=value-sd, ymax=value+sd), width=.1,
                  position=position_dodge(.9)) +
    geom_text(aes(x = VAR, y = value+3*sd, label = n), size = geom.text.size) +
    scale_fill_manual(values = c("ENV" = pal_contribution[2],
                                 "HUM" = pal_contribution[1],
                                 "HAB" = pal_contribution[13],
                                 "BIOT" = pal_contribution[3])) +
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

# plot_data = rf
# fitted_model = "rf"
# color = pal_contribution
# labs_y = ""
# labs_fill = ""
# ylim = c(0,0.1)
# legend.position = "none"



covariates_importance_function <- function(plot_data
                                           ){
  
  require(ggplot2)

  # covariates relative importance by mean

  plot_level <- unique(plot_data$fitted_model)
  plot_data <- plot_data |>
    dplyr::filter(species_name %in% best_models[best_models$best_model == plot_level,1]$species_name)

  ENV <- plot_data |> 
    dplyr::filter(variable %in% env_var) |>
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "ENV",
                     plot_level = plot_level)
  ENV[ENV$variable == "max_1year_analysed_sst",]$variable <- "Sea Surface Temperature (max 1 year)"
  ENV[ENV$variable == "min_1year_analysed_sst",]$variable <- "Sea Surface Temperature (min 1 year)"
  ENV[ENV$variable == "max_5year_degree_heating_week",]$variable <- "Degree Heating Week (max 5 year)"
  ENV[ENV$variable == "mean_1year_nppv",]$variable <- "Net Primary Productivity (mean 1 year)"
  ENV[ENV$variable == "mean_1year_so_mean",]$variable <- "Sea Surface Salinity (mean 1 year)"
  ENV[ENV$variable == "min_5year_ph",]$variable <- "pH (min 5 year)"

  SOC <- plot_data |> 
    dplyr::filter(variable %in% hum_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "HUM",
                     plot_level = plot_level)
  SOC[SOC$variable == "protection_status2",]$variable <- "MPA status"
  SOC[SOC$variable == "gdp",]$variable <- "Growth Development Product"
  SOC[SOC$variable == "gravtot2",]$variable <- "Human Gravity"
  SOC[SOC$variable == "n_fishing_vessels",]$variable <- "Fishing Vessels Density"
  SOC[SOC$variable == "neartt",]$variable <- "Nearest Population"
  SOC[SOC$variable == "marine_ecosystem_dependency",]$variable <- "Marine Ecosystem Dependency"

  HAB <- plot_data |> 
    dplyr::filter(variable %in% hab_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "HAB",
                     plot_level = plot_level)
  HAB[HAB$variable == "Rock_500m",]$variable <- "Rock (%)"
  HAB[HAB$variable == "Sand_500m",]$variable <- "Sand (%)"
  HAB[HAB$variable == "coral_algae_500m",]$variable <- "Coral/Algae (%)"
  HAB[HAB$variable == "coral",]$variable <- "Coral (RLS)"
  HAB[HAB$variable == "depth",]$variable <- "Depth"
  HAB[HAB$variable == "reef_extent",]$variable <- "Reef extent"
  
  BIOT <- plot_data |> 
    dplyr::filter(variable %in% biot_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "BIOT",
                     plot_level = plot_level)
  BIOT[BIOT$variable == "diversity",]$variable <- "Diversity"
  BIOT[BIOT$variable == "max_trophic",]$variable <- "Trophic Level (max)"
  BIOT[BIOT$variable == "mean_trophic",]$variable <- "Trophic Level (mean)"
  BIOT[BIOT$variable == "n_trophic",]$variable <- "Number of Throphic Level"

  cont <- ENV |>  
    dplyr::full_join(HAB) |> 
    dplyr::full_join(SOC) |> 
    dplyr::full_join(BIOT)
  
}

plot_covariates_importance_function <- function(plot_data,
                                                color,
                                                labs_y,
                                                labs_fill,
                                                ylim,
                                                legend.position
){
  
  require(ggplot2)
  
  importance_plot <- ggplot(plot_data) +
    geom_col(aes(x = reorder(variable, value), y = value, fill = VAR)) +
    geom_errorbar(aes(x = variable, y = value, ymin = value - sd, ymax = value + sd), width = .2,
                    position = position_dodge(.9)) +
    scale_fill_manual(values = c("ENV" = color[2],
                                 "HUM" = color[1],
                                 "HAB" = color[13],
                                 "BIOT" = color[3])) +
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


# plot_data = plot_data
# fitted_model = plot_level

var_max_function <- function(plot_data,
                             fitted_model
){
  
  plot_level <- fitted_model
  plot_data <- plot_data |> 
    dplyr::filter(fitted_model == plot_level)
  plot_data <- plot_data |> 
    dplyr::filter(species_name %in% best_models[best_models$best_model == plot_level,1]$species_name)
  
  ENV <- plot_data |> 
    dplyr::filter(variable %in% env_var) |>
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = mean(Dropout_loss),
                     sd = mean(sd_dropout_loss),
                     VAR = "ENV",
                     plot_level = plot_level)
  
  SOC <- plot_data |> 
    dplyr::filter(variable %in% hum_var) |> 
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = mean(Dropout_loss),
                     sd = mean(sd_dropout_loss),
                     VAR = "HUM",
                     plot_level = plot_level)

  HAB <- plot_data |> 
    dplyr::filter(variable %in% hab_var) |> 
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = mean(Dropout_loss),
                     sd = mean(sd_dropout_loss),
                     VAR = "HAB",
                     plot_level = plot_level)
  
  BIOT <- plot_data |> 
    dplyr::filter(variable %in% biot_var) |> 
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = mean(Dropout_loss),
                     sd = mean(sd_dropout_loss),
                     VAR = "BIOT",
                     plot_level = plot_level)

  cont <- ENV |>
    dplyr::full_join(HAB) |>
    dplyr::full_join(SOC) |> 
    dplyr::full_join(BIOT)

  best <- dplyr::tibble(species_name = unique(cont$species_name),
                        ENV = cont[cont$VAR == "ENV",2],
                        SOC = cont[cont$VAR == "HUM",2],
                        HAB = cont[cont$VAR == "HAB",2],
                        BIOT = cont[cont$VAR == "BIOT",2])
  best <- dplyr::inner_join(best, best_models, by = "species_name")
  best <- as.matrix(best)
  
  best <- lapply(1:nrow(best), function(i) {
    test <- best[i,]
    testt <- which.max(test[2:5])
    testtt <- dplyr::tibble(species_name = test[1],
                            varmax = testt,
                            plot_level = test[6])
  })
  
  best <- do.call(rbind, best)
  best$varmax <- as.character(best$varmax)
  best[best$varmax == "1" ,2] <- "ENV"
  best[best$varmax == "2" ,2] <- "HUM"
  best[best$varmax == "3" ,2] <- "HAB"
  best[best$varmax == "4" ,2] <- "BIOT"
  
  if(is.na(fitted_model)){
    
    cont <- cont[which(colnames(cont) %in% "plot_level" == FALSE)] |> 
      dplyr::inner_join(best, by = c("species_name"))
    
  }else{
  
  cont <- cont |> 
    dplyr::inner_join(best, by = c("species_name", "plot_level"))
  
  }
  
  if(is.na(fitted_model)){
    
    var_max <- cont |> 
      dplyr::group_by(species_name, varmax) |> 
      dplyr::summarise(n = dplyr::n()/4)
    var_max <- var_max |> 
      dplyr::rename(VAR = varmax)
    
  }else{
    
    var_max <- cont |> 
      dplyr::group_by(plot_level, varmax) |> 
      dplyr::summarise(n = dplyr::n()/4)
    var_max <- var_max |> 
      dplyr::rename(VAR = varmax)
    
  }
  
}


# plot_data = scaridae
# fitted_model = ""
# color = pal_contribution
# labs_y = ""
# labs_fill = ""
# legend.position = "none"
# mul = 2

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
  
  ENV <- plot_data |> 
    dplyr::filter(variable %in% env_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "ENV",
                     plot_level = plot_level)

  SOC <- plot_data |> 
    dplyr::filter(variable %in% hum_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "HUM",
                     plot_level = plot_level)

  HAB <- plot_data |> 
    dplyr::filter(variable %in% hab_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "HAB",
                     plot_level = plot_level)
  
  BIOT <- plot_data |> 
    dplyr::filter(variable %in% biot_var) |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(value = median(Dropout_loss),
                     sd = median(sd_dropout_loss),
                     VAR = "BIOT",
                     plot_level = plot_level)
  
  cont <- ENV |>  
    dplyr::full_join(HAB) |> 
    dplyr::full_join(SOC) |> 
    dplyr::full_join(BIOT)

  var_max <- var_max_function(plot_data = plot_data,
                              fitted_model = plot_level)
  
  #merge contribution per var and model
  
  cont_merge <- cont |>  
    dplyr::group_by(VAR, plot_level) |> 
    dplyr::summarise(value = mean(value),
                     sd = mean(sd))
  
  cont_merge <- dplyr::inner_join(cont_merge, var_max, by = c("plot_level", "VAR"))
  
}

# plot_data =  merged_covariates_importance_function(plot_data = rf,
#                                                    fitted_model = "rf")
# color = pal_contribution
# labs_y = ""
# labs_fill = ""
# legend.position = "none"
# mul = 2

plot_merged_covariates_importance_function <- function(plot_data,
                                                       fitted_model,
                                                       color,
                                                       labs_y,
                                                       labs_fill,
                                                       legend.position,
                                                       mul
){
  
  merged_importance_plot <- ggplot(plot_data) +
    geom_col(aes(x = reorder(VAR, value), y = value, fill = VAR)) +
    geom_errorbar(aes(x = VAR, y = value, ymin=value-sd, ymax=value+sd), width=.1,
                  position=position_dodge(.9)) +
    geom_text(aes(x = VAR, y = value+mul*sd, label = n), size = 5) +
    scale_fill_manual(values = c("ENV" = color[2],
                                 "HUM" = color[1],
                                 "HAB" = color[13],
                                 "BIOT" = color[3])) +
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


species_covariates_importance_function <- function(plot_data){
  
  require(ggplot2)
  
  # covariates relative importance by median
  
  plot_data <- plot_data |> 
    dplyr::inner_join(best_models)
  plot_data <- plot_data |> 
    dplyr::filter(fitted_model == best_model)
  
  ENV <- plot_data |> 
    dplyr::filter(variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_nppv", "mean_1year_so_mean", "min_1year_analysed_sst", "min_5year_ph", "min_7days_o2")) |>
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = mean(Dropout_loss),
                     sd = mean(sd_dropout_loss),
                     VAR = "ENV")
  
  SOC <- plot_data |> 
    dplyr::filter(variable %in% c("protection_status2", "gdp", "gravtot2", "no_violence", "n_fishing_vessels", "neartt", "marine_ecosystem_dependency")) |> 
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = mean(Dropout_loss),
                     sd = mean(sd_dropout_loss),
                     VAR = "HUM")
  
  HAB <- plot_data |> 
    dplyr::filter(variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "depth", "reef_extent")) |> 
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = mean(Dropout_loss),
                     sd = mean(sd_dropout_loss),
                     VAR = "HAB")
  
  BIOT <- plot_data |> 
    dplyr::filter(variable %in% c("delta_biomass", "diversity", "max_trophic", "mean_biomass", "mean_trophic", "n_trophic")) |> 
    dplyr::group_by(species_name) |> 
    dplyr::summarise(value = mean(Dropout_loss),
                     sd = mean(sd_dropout_loss),
                     VAR = "BIOT")
  

  cont <- ENV |>
    dplyr::full_join(HAB) |>
    dplyr::full_join(SOC) |> 
    dplyr::full_join(BIOT)

}
