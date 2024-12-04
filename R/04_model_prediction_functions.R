
# function to rescale between 0 and 1 ----

rescale_01 = function(x){
  
  (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm=T))
  
  }

# function to unnest large data from as a data.table ----

unnest_dt2 <- function(tbl, ...) {
  
  tbl <- data.table::as.data.table(tbl)

  col <- ensyms(...)

  clnms <- syms(setdiff(colnames(tbl), as.character(col)))
  
  tbl <- data.table::as.data.table(tbl)
  
  tbl <- eval(
    expr(tbl[, lapply(.SD, unlist), by = list(!!!clnms), .SDcols = as.character(col)])
  )
  
  colnames(tbl) <- c(as.character(clnms), as.character(col))
  
  tbl
}


# input_data = read_sp_eco
# nbins = 25
# levels = c("glm", "gam", "spamm", "rf", "gbm", "sprf")

observed_predicted_plot <- function(input_data, 
                                    nbins,
                                    levels){
  
  require(ggplot2)
  require(patchwork)
  
  model_outputs <- input_data
  
  sp <- unique(model_outputs$species_name)
  
  sp_i <- pbmcapply::pbmclapply(1:length(levels), function(i) {
    
      sp_j <- model_outputs[model_outputs$model == levels[i],]
      
      sp_j$validation_predict <- as.numeric(sp_j$validation_predict)
      sp_j$validation_observed <- as.numeric(sp_j$validation_observed)
      
      sp_j <- sp_j[which(is.infinite(sp_j$validation_predict) == FALSE),]
      
      sp_j$validation_predict[sp_j$validation_predict < 0] <- 0
      
      # glm, gam and spamm can give unreal prediction. Remove prediction 20 times higher than the maximum biomass observed.
      which_to_high <- which(sp_j$validation_predict > (max(sp_j$validation_observed) * 20))
      
      if(length(which_to_high) != 0) {
        
        sp_j <- sp_j[-which_to_high,]
        
      }else{
        
        sp_j <- sp_j
        
      }
      
      sp_j$validation_predict <- rescale_01(log10(sp_j$validation_predict + 1))
      sp_j$validation_observed <- rescale_01(log10(sp_j$validation_observed + 1))
      
      sp_j

  }, mc.cores = parallel::detectCores() - 1)
  
  model_outputs <- do.call(rbind, sp_i)

  # create a transformations label
  model_outputs$transformation <- model_outputs$model
  model_outputs$predicted <- model_outputs$validation_predict
  model_outputs$observed <- model_outputs$validation_observed
  
  # truncate the data to 99th percentiles
  model_outputs$predicted[model_outputs$predicted > as.numeric(quantile(model_outputs$predicted, 0.99, na.rm =T))] <- quantile(model_outputs$predicted, 0.99, na.rm =T)
  model_outputs$predicted[model_outputs$predicted < as.numeric(quantile(model_outputs$predicted, 0.01, na.rm =T))] <- quantile(model_outputs$predicted, 0.01, na.rm =T)
  
  #  aggregate at a species level to ensure species with more data don't influence the distribution of results
  model_outputs <- model_outputs |> 
    dplyr::mutate(observed = plyr::round_any(observed, 1/nbins)) |> 
    dplyr::group_by(model, species_name, observed) |> 
    dplyr::do(predicted = mean(.$predicted, na.rm = T)) |> 
    tidyr::unnest(predicted)
  
  # plots for all plot levels ----
  plot_levels_plot <- list()
  plot_levels_plot <- lapply(1:length(levels), function(i) {
    # create basic plot
    base_plot <- model_outputs |> 
      dplyr::filter(model == levels[i]) |> 
      ggplot(aes(x = observed, y = predicted)) +
      geom_density_2d_filled(aes(x = observed, y = predicted), contour_var = 'count', 
                             contour = F, n = 100, bins = nbins, colour = 'transparent') + 
      ylim(0, 1) +
      scale_fill_viridis_d(option = 'viridis', begin = 0.1, name = 'Count') +
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            strip.background = element_rect(fill = 'grey90', colour = 'grey90'), 
            aspect.ratio = 1)
    
    # add faceting to plot level
    facet_plot <- base_plot + 
      facet_grid(~model) + #de base facet_grid, sans ncol
      geom_abline(slope = 1, intercept = 0)
    
    plot_levels_plot[[i]] <- facet_plot + 
      labs(x = "Observed", y = "Predicted") +
      theme(legend.position = "none",
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 20),
            legend.text = element_text(size = 5), 
            legend.title = element_text(size = 5),
            strip.text.x = element_text(size = 15),
            strip.text.y = element_text(size = 15),
            strip.background = element_blank(),
            panel.background = element_rect(fill = "white", colour = "grey50",
                                            size = 1, linetype = "solid"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    
  })
  
  all_plots <- (plot_levels_plot[[1]] + plot_levels_plot[[2]]) / (plot_levels_plot[[3]] + plot_levels_plot[[4]]) / (plot_levels_plot[[5]] + plot_levels_plot[[6]])

  ggsave("figures/all_predictions_pres.png", all_plots, width = 9, height = 11)
  
}



observed_predicted_best_plot <- function(input_data,
                                         nbins){

  require(ggplot2)
  require(patchwork)

  model_outputs <- input_data

  sp <- unique(model_outputs$species_name)

  model_j <- lapply(1:length(unique(model_outputs$model)), function(j) {
    
    sp_j <- model_outputs[model_outputs$model == unique(model_outputs$model)[j],]
    
    sp_j$validation_predict <- as.numeric(sp_j$validation_predict)
    sp_j$validation_observed <- as.numeric(sp_j$validation_observed)
    
    sp_j$validation_predict[sp_j$validation_predict < 0] <- 0
    
    sp_j$validation_predict <- rescale_01(log10(sp_j$validation_predict + 1))
    sp_j$validation_observed <- rescale_01(log10(sp_j$validation_observed + 1))
    
    return(sp_j)
    
  })
  
  model_outputs <- do.call(rbind, model_j)

  # create a transformations label
  model_outputs$transformation <- model_outputs$model
  model_outputs$predicted <- as.numeric(model_outputs$validation_predict)
  model_outputs$observed <- as.numeric(model_outputs$validation_observed)

  #  aggregate at a species level to ensure species with more data don't influence the distribution of results
  model_outputs <- model_outputs |>
    dplyr::mutate(observed = plyr::round_any(observed, 1/nbins)) |>
    dplyr::group_by(model, species_name, observed) |>
    dplyr::do(predicted = mean(.$predicted, na.rm = T)) |>
    tidyr::unnest(predicted)

  base_plot <- model_outputs |>
    ggplot(aes(x = observed, y = predicted)) +
    stat_density2d_filled(n = 50, bins = nbins) +
    ylim(0, 1) +
    scale_fill_viridis_d(option = 'viridis', begin = 0.1, name = 'Count') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = 'grey90', colour = 'grey90'),
          aspect.ratio = 1) +
    geom_abline(slope = 1, intercept = 0)

  plot_levels_plot <- base_plot +
    labs(x = "Observed", y = "Predicted") +
    theme(legend.position = "none",
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20),
          strip.background = element_blank(),
          panel.background = element_rect(fill = "white", colour = "grey50",
                                          size = 1, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave("figures/all_predictions_pred_best.png", plot_levels_plot, width = 7, height = 7)

}
