# function for building species traits figure

# data = cont
# trait = trait

kruskal_test_function <- function(data,
                                  trait
                                  ){
  
  data$newvar <- paste0(data$VAR, '_', unlist(data[trait]))
  
  res.kruskal <- with(data, agricolae::kruskal(value, newvar, p.adj = "bonferroni", group = FALSE))
  res.kruskal2 <- with(data, agricolae::kruskal(value, newvar, p.adj = "bonferroni", group = TRUE))
  
  res.kruskal2$groups <- data.frame(VAR = stringr::word(row.names(res.kruskal2$groups), 1, sep = "_"),
                                    trait = stringr::word(row.names(res.kruskal2$groups), 2, sep = "_"),
                                    value = res.kruskal2$groups$value,
                                    groups = res.kruskal2$groups$groups)
  names(res.kruskal2$groups)[2] <- trait
  
  n <- which(colnames(data) == trait)

  labs.position <- data |> 
    dplyr::group_by_at(c(n, 4)) |> 
    dplyr::summarise(mean = mean(value), quant = quantile(value, probs = 0.75))
  
  # res.kruskal2$groups <- dplyr::inner_join(res.kruskal2$groups, labs.position, by = c("VAR", trait))
  res.kruskal2$groups <- dplyr::inner_join(labs.position, res.kruskal2$groups, by = c("VAR", trait))
  res.kruskal2
  
}

# plot_data = bind_files
# data_trait = phylo
# trait = "family"
# color = pal_sp_trait
# labs_title = "family"
# x_text_angle = 65
# legend.position = c(0.85, 0.8)

species_traits_function <- function(plot_data,
                                    data_trait,
                                    trait,
                                    color,
                                    labs_title,
                                    x_text_angle,
                                    legend.position
){
  
  require(ggplot2)
  
  cont <- lapply(1:length(fitted_model), function(i) {
    
    plot_level <- fitted_model[i]
    plot_data <- plot_data |> 
      dplyr::filter(fitted_model == plot_level)
    plot_data <- plot_data |>
      dplyr::filter(species_name %in% best_models[best_models$best_model == plot_level,1]$species_name)
    
    ENV <- plot_data |> 
      dplyr::filter(variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "max_5year_nppv", "min_1year_analysed_sst", "min_5year_ph", "min_7days_o2")) |> 
      dplyr::group_by(species_name) |> 
      dplyr::summarise(value = median(Dropout_loss),
                       sd = median(sd_dropout_loss),
                       VAR = "ENV",
                       plot_level = plot_level)
    
    SOC <- plot_data |> 
      dplyr::filter(variable %in% c("effectiveness", "gdp", "gravtot2", "no_violence", "n_fishing_vessels", "natural_ressource_rent", "neartt", "marine_ecosystem_dependency")) |> 
      dplyr::group_by(species_name) |> 
      dplyr::summarise(value = median(Dropout_loss),
                       sd = median(sd_dropout_loss),
                       VAR = "HUM",
                       plot_level = plot_level)
    
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
  
  cont <- cont |> 
    dplyr::inner_join(data_trait, by = "species_name", multiple = "first")
  
  n_trait <- cont |>
    dplyr::group_by(.dots = trait) |> 
    dplyr::summarise(n = dplyr::n() / 3)
  
  kruskal_test_trait <- kruskal_test_function(cont,
                                              trait)
  
  cont[colnames(cont) %in% trait] <- sapply(1:nrow(cont[colnames(cont) %in% trait]), function(i) {
    
    row_i <- cont[colnames(cont) %in% trait][i,]
    
    which_row <- which(grepl(row_i, unlist(n_trait[colnames(n_trait) %in% trait])) == TRUE)
    
    paste0(row_i, " (n = ", n_trait$n[which_row], ")")
    
  })
  
  kruskal_test_trait$groups[colnames(kruskal_test_trait$groups) %in% trait] <- sapply(1:nrow(kruskal_test_trait$groups[colnames(kruskal_test_trait$groups) %in% trait]), function(i) {
    
    row_i <- kruskal_test_trait$groups[colnames(kruskal_test_trait$groups) %in% trait][i,]
    
    which_row <- which(grepl(row_i, unlist(n_trait[colnames(n_trait) %in% trait])) == TRUE)
    
    paste0(row_i, " (n = ", n_trait$n[which_row], ")")
    
  })
  
  if(kruskal_test_trait$statistics$p.chisq > 0.05){
    
    stop(print("No statistical differences among groups"))
    
  }else{
    
    print("p.chisq < 0.05")
    
    plot_trait <- cont |> 
      ggplot(aes_string(x = trait, y = "value", fill = "VAR")) +
      geom_boxplot(aes(fill = factor(VAR)), width = 0.6, outlier.shape = NA, position = position_dodge(width = 0.75)) +
      scale_fill_manual(values = c("ENV" = color[1],
                                   "HUM" = color[3],
                                   "HAB" = color[2])) +
      theme_bw() +
      geom_text(data = kruskal_test_trait$groups, aes_string(x = as.factor(unlist(kruskal_test_trait$groups[,1])), y = kruskal_test_trait$groups$quant, label = "groups"), vjust=-0.48, size = 5, position = position_dodge(width = 0.85)) +
      coord_cartesian(ylim = c(0,0.2)) +
      labs(y = "Relative importance (RMSE)", x = "", title = labs_title, fill = "") +
      theme(legend.position = legend.position,
            legend.direction = "horizontal",
            legend.background = element_rect(fill = "white"),
            legend.key = element_rect(fill = "white", color = NA),
            title = element_text(size = 20),
            axis.text = element_text(size = 20),
            axis.text.x = element_text(size = 18),
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

    if(!is.null(x_text_angle)) {
      
      plot_trait <- plot_trait +
        theme(axis.text.x = element_text(angle = x_text_angle, hjust = 1))
      
    }
    
  }
  
  return(plot_trait)
  
}
