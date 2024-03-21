
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

# function for hex density plots of observed vs. predicted abundance ----

observed_predicted_plot <- function(input_data, 
                                    nbins = 20,
                                    levels = c('GLM', 'GAM', 'SPAMM', 'RF', 'GBM', 'SPRF')){
  
  require(ggplot2)
  require(patchwork)
  
  model_outputs <- input_data
  model_outputs <- lapply(1:nrow(model_outputs), function(i) {
    
    sp_i <- model_outputs[i,]
    sp_i_dtf <- data.frame(predict = as.numeric(unlist(sp_i$validation_predict)),
                           observe = as.numeric(unlist(sp_i$validation_observed)))
    
    sp_i_dtf <- sp_i_dtf[which(is.na(sp_i_dtf$predict) == FALSE),]
    
    sp_i_dtf <- sp_i_dtf[which(is.infinite(sp_i_dtf$predict) == FALSE),]
    
    sp_i_dtf[sp_i_dtf$predict == -1,] <- 0
    
    sp_i$validation_predict <- list(sp_i_dtf$predict)
    sp_i$validation_observed <- list(sp_i_dtf$observe)
    
    return(sp_i)

  })
  model_outputs <- do.call(rbind, model_outputs)

  # remove turn values that are <0 to NAs
  # rescale predictions
  model_outputs$predicted <- pbapply::pblapply(model_outputs$validation_predict, function(x) return(ifelse(x < 0, x, x)))
  model_outputs$predicted <- pbapply::pblapply(1:nrow(model_outputs), function(x) rescale_01(log10(unlist(model_outputs[x,]$predicted)+1)))

  # rescale observations
  model_outputs$observed <- pbapply::pblapply(model_outputs$validation_observed, function(x) rescale_01(log10(x+1)))
  
  # if fewer samples just unnest it directly
  model_outputs <- dplyr::as_data_frame(unnest_dt2(data.table::data.table(model_outputs |>  dplyr::select(-validation_observed, -validation_predict)), predicted, observed))

  # create a transformations label
  model_outputs$transformation <- model_outputs$fitted_model
  
  # truncate the data to 99th percentiles
  model_outputs$predicted[model_outputs$predicted > as.numeric(quantile(model_outputs$predicted, 0.99, na.rm =T))] <- quantile(model_outputs$predicted, 0.99, na.rm =T)
  model_outputs$predicted[model_outputs$predicted < as.numeric(quantile(model_outputs$predicted, 0.01, na.rm =T))] <- quantile(model_outputs$predicted, 0.01, na.rm =T)
  
  #  aggregate at a species level to ensure species with more data don't influence the distribution of results
  model_outputs <- model_outputs |> 
    dplyr::mutate(observed = plyr::round_any(observed, 1/nbins)) |> 
    dplyr::group_by(fitted_model, species_name, observed) |> 
    dplyr::do(predicted = mean(.$predicted, na.rm = T)) |> 
    tidyr::unnest(predicted)

  # plots for all plot levels ----
  plot_levels_plot <- list()
  plot_levels_plot <- lapply(1:length(levels), function(i) {
    # create basic plot
    base_plot <- model_outputs |> 
      dplyr::filter(fitted_model == levels[i]) |> 
      ggplot(aes(x = observed, y = predicted)) +
      geom_density_2d_filled(aes(x = observed, y = predicted), contour_var = 'count', 
                             contour = F, n = 100, bins= 10, colour = 'transparent') + 
      scale_fill_viridis_d(option = 'viridis', begin = 0.2, end = 0.9,
                           name = 'Count') +
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            strip.background = element_rect(fill = 'grey90', colour = 'grey90'), 
            aspect.ratio = 1)

    # add faceting to plot level
    facet_plot <- base_plot + 
      facet_grid(~fitted_model) + #de base facet_grid, sans ncol
      geom_abline(slope = 0.95, intercept = 0)
    
    plot_levels_plot[[i]] <- facet_plot + 
      labs(x = "Observed", y = "Predicted") +
      theme(legend.position = "none",
            axis.text=element_text(size=15),
            axis.title=element_text(size=25),
            legend.text=element_text(size=10), 
            legend.title=element_text(size=10),
            strip.text.x = element_text(size = 20),
            strip.text.y = element_text(size = 20),
            strip.background = element_blank(),
            panel.background = element_rect(fill = "white", colour = "grey50",
                                            size = 1, linetype = "solid"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    
  })
  
  all_plots <- (plot_levels_plot[[1]] + plot_levels_plot[[2]]) / (plot_levels_plot[[3]] + plot_levels_plot[[4]]) / (plot_levels_plot[[5]] + plot_levels_plot[[6]])
  ggsave("figures/all_predictions.pdf", all_plots, width = 11, height = 15)
}
