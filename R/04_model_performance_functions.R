# function for producing a common scale of assessment criteria ----

aggregate_metrics <- function(plot_data,
                              metrics = c("Intercept", "Slope", "Pearson", "Spearman"),
                              targets = list(Intercept  = c(0, -5, 5), 
                                             Slope      = c(1,  0, 2), 
                                             Pearson    = c(1,  0, 1), 
                                             Spearman   = c(1,  0, 1)), 
                              levels = c("GLM", "GAM", "GBM", "SPRF")){
  
  # get metrics
  metrics_data <- na.omit(plot_data[metrics])
  
  # match targets to relative scale if < 0
  metrics_data[, which(names(metrics_data) %in% c("Slope"))] <- abs(metrics_data[, which(names(metrics_data) %in% c("Slope"))] - 1)
  metrics_data[, which(names(metrics_data) %in% c("Intercept"))] <- abs(metrics_data[, which(names(metrics_data) %in% c("Intercept"))])
  
  # rank order and rescale
  rescale_01 <- function(x){(x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))}
  metrics_data <- data.frame(sapply(metrics_data, function(x) rescale_01(rank(x))))
  
  # range inversion
  invert_range <- function(x){ (max(x, na.rm = T) + min(x, na.rm = T)) - x}
  metrics_data[which(names(metrics_data) %in% c("Intercept", "Slope"))] <-
    sapply(metrics_data[which(names(metrics_data) %in% c("Intercept", "Slope"))], function(x) invert_range(x))

  # create columns for plotting

  if(length(which(names(metrics_data) %in% c("Intercept", "Slope", "Pearson", "Spearman"))) == 1){
    discrimination <- metrics_data[, which(names(metrics_data) %in% c("Intercept", "Slope", "Pearson", "Spearman"))]}else{
    discrimination <- rowSums(metrics_data[, which(names(metrics_data) %in% c("Intercept", "Slope", "Pearson", "Spearman"))])/length(which(names(metrics_data) %in% c("Intercept", "Slope", "Pearson", "Spearman")))}

  # index for ANY NAs
  NA_index <- rowSums(sapply(plot_data[metrics], is.na)) > 0
  
  plot_data$discrimination <- NA
  
  if(sum(NA_index) != 0){

  plot_data$discrimination[-which(NA_index)] <- discrimination
  
  }else{
    
    plot_data$discrimination <- discrimination

  }
  
  return(plot_data)
  
}

performance_plot <- function(plot_data,
                             metrics_sel, #options are pearson, spearman, slope and intercept
                             slope,
                             intercept,
                             color,
                             ylim,
                             legend.position,
                             plot_title
                             ){
  
  require(ggplot2)
  
  plot_perf <- plot_data |> 
    dplyr::filter(metrics == metrics_sel) |> 
    dplyr::mutate(model = forcats::fct_relevel(model, "GLM", "GAM", "SPAMM","RF", "GBM", "SPRF")) |> 
    ggplot(aes(x = model, y = value, fill = cat)) +
    geom_boxplot(outlier.shape = NA) +
    geom_abline(slope = slope, intercept = intercept, linetype = 2, size = 1.25, color = "red") +
    scale_fill_manual(values = c("Best models" = pal_perf[6],
                                 "All models" = pal_perf[2])) +
    theme_bw() +
    coord_cartesian(ylim = ylim) + 
    theme(legend.position = legend.position,
          legend.direction = "horizontal") +
    labs(y = metrics_sel, x = "", title = plot_title, fill = "") +
    theme(title = element_text(size=20),
          axis.text=element_text(size=15),
          axis.title=element_text(size=25),
          legend.text=element_text(size=20), 
          legend.title=element_text(size=25),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20),
          strip.background = element_blank(),
          panel.background = element_rect(fill = "white", colour = "grey50",
                                          size = 1, linetype = "solid"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  print(plot_perf)
  
}

# data = best_trait
# traits = "Water.column"
# metrics = "Spearman"

kruskal_test_function_perf <- function(data,
                                       traits,
                                       metrics){
  
  data$trait <- unlist(data[,traits])
  data$metric <- unlist(data[,metrics])
  
  res.kruskal <- with(data, agricolae::kruskal(metric, trait, p.adj = "bonferroni", group = FALSE))
  res.kruskal2 <- with(data, agricolae::kruskal(metric, trait, p.adj = "bonferroni", group = TRUE))
  
  res.kruskal2$groups <- data.frame(trait = stringr::word(row.names(res.kruskal2$groups), 1, sep = "_"),
                                    value = res.kruskal2$groups$metric,
                                    groups = res.kruskal2$groups$groups)
  # names(res.kruskal2$groups)[1] <- trait
  
  n <- which(colnames(data) == traits)
  
  labs.position <- data |> 
    dplyr::group_by(trait) |> 
    dplyr::summarise(mean = mean(metric),
                     quant = quantile(metric, probs = 0.75))

  res.kruskal2$groups <- dplyr::inner_join(res.kruskal2$groups, labs.position, by = "trait")
  res.kruskal2
  
}
