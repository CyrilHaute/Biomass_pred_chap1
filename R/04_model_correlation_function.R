# functions to plot correlation plot

custom_ggally_cor <- function(data, mapping, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  cor_value <- cor(x, y, use = "complete.obs")
  cor_label <- round(cor_value, 1)
  

  plot_data <- data.frame(
    xmin = 0.30, xmax = 0.9,
    ymin = 0.30, ymax = 0.9,
    cor = cor_value
  )
  
  ggplot(plot_data, aes(fill = cor_value)) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
    geom_text(aes(x = 0.60, y = 0.60, label = cor_label)) +
    scale_fill_gradient2(
      low = "#3B9AB2", mid = "#EEEEEE", high = "#F21A00",
      midpoint = 0.6, limits = c(0.35, 0.85),
      name = "Corrélation"
    ) +
    theme_void() +
    theme(legend.position = "none")
}

custom_diag_text <- function(data, mapping, ...){
  
  var_name <- as_label(mapping$x)
  
  plot_data <- data.frame(
    cor = var_name
  )
  
  ggplot(data = plot_data) +
    geom_text(aes(x = 0.5, y = 0.5, label = cor)) +
    theme_void()
}

custom_ggally_scat <- function(data, mapping, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  ggplot(data = data, mapping = mapping) +
    geom_point(aes(x = x, y = y), size = 0.3) +
    theme_minimal()
  
}