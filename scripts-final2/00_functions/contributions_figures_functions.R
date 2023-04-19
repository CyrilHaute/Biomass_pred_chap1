
covariates_importance_function <- function(plot_data,
                                           fitted_model,
                                           color,
                                           labs_y,
                                           labs_fill,
                                           legend.position
                                           ){

  # covariates relative importance by median

  plot_level <- fitted_model
  only_model <- plot_data %>% filter(fitted_model == plot_level)
  only_model <- only_model %>% filter(species_name %in% best_models[best_models$best_model == plot_level,1]$species_name)
    
  ENV <- lapply(1:nrow(only_model), function(i) { only_model$contributions[[i]][c(9:15)]})
  ENV_sd <- lapply(1:nrow(only_model), function(i) { only_model$sd_contributions[[i]][c(9:15)]})
  ENV <- do.call(rbind, ENV)
  ENV_sd <- do.call(rbind, ENV_sd)
  ENV <- tibble(value = colMedians(ENV),
                sd = colMedians(ENV_sd),
                var = c("max_sst_1year", "mean_chl_1year", "mean_DHW_1year", "mean_DHW_5year", "mean_pH_1year", "mean_sss_1year", "min_sst_1year"),
                VAR = rep("ENV",7),
                plot_level = rep(plot_level, 7))
    
  SOC <- lapply(1:nrow(only_model), function(i) { only_model$contributions[[i]][c(1,5:8,16:17)]})
  SOC_sd <- lapply(1:nrow(only_model), function(i) { only_model$sd_contributions[[i]][c(1,5:8,16:17)]})
  SOC <- do.call(rbind, SOC)
  SOC_sd <- do.call(rbind, SOC_sd)
  SOC <- tibble(value = colMedians(SOC),
                sd = colMedians(SOC_sd),
                var = c("Conflicts", "MPA Effectiveness", "GDP", "Gravity", "Human Development Index", "Natural Ressource", "Neartt"),
                VAR = rep("HUM",7),
                plot_level = rep(plot_level, 7))
    
  HAB <- lapply(1:nrow(only_model), function(i) { only_model$contributions[[i]][c(2:4,18:21)]})
  HAB_sd <- lapply(1:nrow(only_model), function(i) { only_model$sd_contributions[[i]][c(2:4,18:21)]})
  HAB <- do.call(rbind, HAB)
  HAB_sd <- do.call(rbind, HAB_sd)
  HAB <- tibble(value = colMedians(HAB),
                sd = colMedians(HAB_sd),
                var = c("Coral", "Coral_Algea(Allen)", "Depth", "Rock(Allen)", "Rubble(Allen)", "Sand(Allen)", "Reef extent(Allen)"),
                VAR = rep("HAB",7),
                plot_level = rep(plot_level,7))
    
  cont <- ENV %>% 
    full_join(HAB) %>% 
    full_join(SOC)
  cont$plot_level <- rep(plot_level, nrow(cont))

  importance_plot <- ggplot(cont) +
    geom_col(aes(x = reorder(var, value), y = value, fill = VAR)) +
    geom_errorbar(aes(x = var, y = value, ymin=value-sd, ymax=value+sd), width=.2,
                    position=position_dodge(.9)) +
    scale_fill_manual(values = c("ENV" = color[1],
                                 "HUM" = color [3],
                                 "HAB" = color [2])) +
    theme_bw() +
    coord_flip(ylim = c(0,0.22)) +
    facet_grid(~plot_level) +
    labs(y = labs_y, x = "", fill = labs_fill) +
    theme(legend.position = legend.position) +
    theme(title = element_text(size=18),
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=15),
          axis.title=element_text(size=18),
          legend.text=element_text(size=15),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20),
          strip.background = element_blank(),
          panel.background = element_rect(fill = "white", colour = "grey50",
                                          size = 1, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

}

var_max_function <- function(plot_data,
                             fitted_model
                             ){
  
  plot_level <- fitted_model
  only_model <- Contributions_biomass %>% filter(fitted_model == plot_level)
  only_model <- only_model %>% filter(species_name %in% best_models[best_models$best_model == plot_level,1]$species_name)
  only_model <- only_model[which(sapply(1:nrow(only_model), function(i) {is.null(only_model$contributions[[i]])}) == FALSE),]
  
  ENV <- lapply(1:nrow(only_model), function(i) { only_model$contributions[[i]][c(9:15)]})
  ENV_sd <- lapply(1:nrow(only_model), function(i) { only_model$sd_contributions[[i]][c(9:15)]})
  ENV <- do.call(rbind, ENV)
  ENV_sd <- do.call(rbind, ENV_sd)
  ENV <- tibble(species_name = only_model$species_name,
                value = rowMedians(ENV),
                sd = rowMedians(ENV_sd),
                var = rep("ENV",nrow(ENV)),
                plot_level = rep(plot_level, nrow(ENV)))
  
  SOC <- lapply(1:nrow(only_model), function(i) { only_model$contributions[[i]][c(1,5:8,16:17)]})
  SOC_sd <- lapply(1:nrow(only_model), function(i) { only_model$sd_contributions[[i]][c(1,5:8,16:17)]})
  SOC <- do.call(rbind, SOC)
  SOC_sd <- do.call(rbind, SOC_sd)
  SOC <- tibble(species_name = only_model$species_name,
                value = rowMedians(SOC),
                sd = rowMedians(SOC_sd),
                var = rep("HUM",nrow(SOC)),
                plot_level = rep(plot_level, nrow(SOC)))
  
  HAB <- lapply(1:nrow(only_model), function(i) { only_model$contributions[[i]][c(2:4,18:21)]})
  HAB_sd <- lapply(1:nrow(only_model), function(i) { only_model$sd_contributions[[i]][c(2:4,18:21)]})
  HAB <- do.call(rbind, HAB)
  HAB_sd <- do.call(rbind, HAB_sd)
  HAB <- tibble(species_name = only_model$species_name,
                value = rowMedians(HAB),
                sd = rowMedians(HAB_sd),
                var = rep("HAB",nrow(HAB)),
                plot_level = rep(plot_level,nrow(HAB)))
  
  cont <- ENV %>% 
    full_join(HAB) %>% 
    full_join(SOC)
  cont$plot_level <- rep(plot_level, nrow(cont))
  cont
  
  best <- tibble(species_name = unique(cont$species_name),
                 ENV = cont[cont$var == "ENV",2],
                 SOC = cont[cont$var == "HUM",2],
                 HAB = cont[cont$var == "HAB",2])
  best <- inner_join(best, best_models, by = "species_name")
  best <- as.matrix(best)
  
  best <- lapply(1:nrow(best), function(i) {
    test <- best[i,]
    testt <- which.max(test[2:4])
    testtt <- tibble(species_name = test[1],
                     varmax = testt,
                     plot_level = test[5])
  })
  
  best <- do.call(rbind, best)
  best$varmax <- as.character(best$varmax)
  best[best$varmax == "1" ,2] <- "ENV"
  best[best$varmax == "2" ,2] <- "HUM"
  best[best$varmax == "3" ,2] <- "HAB"
  
  cont <- cont %>% inner_join(best, by = c("species_name", "plot_level"))
  
  var_max <- cont %>% group_by(plot_level, varmax) %>% summarise(n = n()/3)
  var_max <- var_max %>% rename(VAR = varmax)
  
}

merged_covariates_importance_function <- function(plot_data,
                                                  fitted_model,
                                                  color,
                                                  labs_y,
                                                  labs_fill,
                                                  legend.position,
                                                  mul
                                                  ){

  # covariates relative importance by median
  
  plot_level <- fitted_model
  only_model <- plot_data %>% filter(fitted_model == plot_level)
  only_model <- only_model %>% filter(species_name %in% best_models[best_models$best_model == plot_level,1]$species_name)
  
  ENV <- lapply(1:nrow(only_model), function(i) { only_model$contributions[[i]][c(9:15)]})
  ENV_sd <- lapply(1:nrow(only_model), function(i) { only_model$sd_contributions[[i]][c(9:15)]})
  ENV <- do.call(rbind, ENV)
  ENV_sd <- do.call(rbind, ENV_sd)
  ENV <- tibble(value = colMedians(ENV),
                sd = colMedians(ENV_sd),
                var = c("max_sst_1year", "mean_chl_1year", "mean_DHW_1year", "mean_DHW_5year", "mean_pH_1year", "mean_sss_1year", "min_sst_1year"),
                VAR = rep("ENV",7),
                plot_level = rep(plot_level, 7))
  
  SOC <- lapply(1:nrow(only_model), function(i) { only_model$contributions[[i]][c(1,5:8,16:17)]})
  SOC_sd <- lapply(1:nrow(only_model), function(i) { only_model$sd_contributions[[i]][c(1,5:8,16:17)]})
  SOC <- do.call(rbind, SOC)
  SOC_sd <- do.call(rbind, SOC_sd)
  SOC <- tibble(value = colMedians(SOC),
                sd = colMedians(SOC_sd),
                var = c("Conflicts", "MPA Effectiveness", "GDP", "Gravity", "Human Development Index", "Natural Ressource", "Neartt"),
                VAR = rep("HUM",7),
                plot_level = rep(plot_level, 7))
  
  HAB <- lapply(1:nrow(only_model), function(i) { only_model$contributions[[i]][c(2:4,18:21)]})
  HAB_sd <- lapply(1:nrow(only_model), function(i) { only_model$sd_contributions[[i]][c(2:4,18:21)]})
  HAB <- do.call(rbind, HAB)
  HAB_sd <- do.call(rbind, HAB_sd)
  HAB <- tibble(value = colMedians(HAB),
                sd = colMedians(HAB_sd),
                var = c("Coral", "Coral_Algea(Allen)", "Depth", "Rock(Allen)", "Rubble(Allen)", "Sand(Allen)", "Reef extent(Allen)"),
                VAR = rep("HAB",7),
                plot_level = rep(plot_level,7))
  
  cont <- ENV %>% 
    full_join(HAB) %>% 
    full_join(SOC)
  cont$plot_level <- rep(plot_level, nrow(cont))
  
  var_max <- var_max_function(plot_data = Contributions_biomass,
                              fitted_model = plot_level)
  
  #merge contribution per var and model
  
  cont_merge <- cont %>% 
    group_by(VAR, plot_level) %>% 
    summarise(value = mean(value),
              sd = mean(sd))
  
  cont_merge <- inner_join(cont_merge, var_max, by = c("plot_level", "VAR"))
  
  merged_importance_plot <- ggplot(cont_merge) +
    geom_col(aes(x = reorder(VAR, value), y = value, fill = VAR)) +
    geom_errorbar(aes(x = VAR, y = value, ymin=value-sd, ymax=value+sd), width=.1,
                  position=position_dodge(.9)) +
    geom_text(aes(x = VAR, y = value+mul*sd, label = n), size = 5) +
    scale_fill_manual(values = c("ENV" = color[1],
                                 "HUM" = color[3],
                                 "HAB" = color[2])) +
    theme_bw() +
    coord_flip(ylim = c(0,0.13)) +
    facet_grid(~plot_level) +
    labs(y = labs_y, x = "", fill = labs_fill) +
    theme(legend.position = legend.position) +
    theme(title = element_text(size=20),
          axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title=element_text(size=20),
          legend.text=element_text(size=20),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20),
          strip.background = element_blank(),
          panel.background = element_rect(fill = "white", colour = "grey50",
                                          size = 1, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

}
