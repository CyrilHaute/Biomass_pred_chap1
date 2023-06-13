# function for building species traits figure

kruskal_test_function <- function(data,
                                  trait
                                  ){
  
  data$newvar <- paste0(data$var, '_', unlist(data[trait]))
  
  res.kruskal <- with(data,kruskal(value, newvar,p.adj="bonferroni",group=FALSE))
  res.kruskal2 <- with(data,kruskal(value, newvar,p.adj="bonferroni",group=TRUE))
  
  res.kruskal2$groups <- data.frame(var = word(row.names(res.kruskal2$groups), 1, sep = "_"),
                                    trait = word(row.names(res.kruskal2$groups), 2, sep = "_"),
                                    value = res.kruskal2$groups$value,
                                    groups = res.kruskal2$groups$groups)
  names(res.kruskal2$groups)[2] <- trait
  
  n <- which(colnames(data) == trait)

  labs.position <- data %>%
    group_by_at(c(4,n)) %>%
    summarise(mean=mean(value), quant = quantile(value, probs = 0.75))
  
  res.kruskal2$groups <- inner_join(res.kruskal2$groups, labs.position, by = c("var", trait))
  res.kruskal2
  
}


species_traits_function <- function(plot_data,
                                    trait,
                                    color,
                                    labs_title
                                    ){
  
  cont <- mclapply(1:length(fitted_model), function(i) {
    
    plot_level <- fitted_model[i]
    only_model <- Contributions_biomass %>% filter(fitted_model == plot_level)
    only_model <- only_model %>% filter(species_name %in% best_models[best_models$best_model == plot_level,1]$species_name)
    only_model <- only_model[which(sapply(1:nrow(only_model), function(i) {is.null(only_model$contributions[[i]])}) == FALSE),]
    
    sapply(1:nrow(only_model), function(i) {is.null(only_model$contributions[[i]])})
    
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
    
  }, mc.cores = 1)
  
  cont <- do.call(rbind, cont)
  
  cont <- cont %>% inner_join(sp_car, by = "species_name")
  
  kruskal_test_trait <- kruskal_test_function(cont,
                                              trait)
  
  if(kruskal_test_trait$statistics$p.chisq > 0.05){
    
    stop(print("No statistical differences among groups"))
    
    }else{
      
      print("p.chisq < 0.05")
      
      plot_trait <- cont %>%
        mutate(var = fct_relevel(var, "ENV", "HUM", "HAB")) %>% 
        ggplot(aes_string(x = trait, y = "value", fill = "var")) +
        geom_boxplot(aes(fill = factor(var)), width=0.6, outlier.shape = NA) +
        scale_fill_manual(values = c("ENV" = color[1],
                                     "HUM" = color[3],
                                     "HAB" = color[2])) +
        theme_bw() +
        geom_text(data = kruskal_test_trait$groups, aes_string(x = trait, y = "quant", label = "groups"), size = 7, vjust=-0.5, hjust=-0.55, position = position_dodge(width = 0.6)) +
        coord_cartesian(ylim = c(0,0.25)) +
        labs(y = "Change in RMSE", x = "", title = labs_title) +
        theme(legend.position = 'none',
              title = element_text(size=20),
              axis.text=element_text(size=20),
              axis.text.x = element_text(size = 20),
              axis.text.y = element_text(size = 20),
              axis.title=element_text(size=20),
              legend.text=element_text(size=20),
              legend.title=element_text(size=20),
              strip.text.x = element_text(size = 20),
              strip.text.y = element_text(size = 20),
              strip.background = element_blank(),
              panel.background = element_rect(fill = "white", colour = "grey50",
                                              size = 1, linetype = "solid"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
    }
}