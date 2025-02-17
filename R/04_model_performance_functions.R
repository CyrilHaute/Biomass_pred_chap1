# function for producing a common scale of assessment criteria

aggregate_metrics <- function(plot_data,
                              metrics){
  
  # get metrics
  metrics_data <- na.omit(plot_data[metrics])

  # rank order and rescale
  rescale_01 <- function(x){(x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))}
  metrics_data <- data.frame(sapply(metrics_data, function(x) rescale_01(rank(x))))

  # create columns for plotting

  if(length(which(names(metrics_data) %in% c("pearson", "spearman"))) == 1){
    discrimination <- metrics_data[, which(names(metrics_data) %in% c("pearson", "spearman"))]}else{
      discrimination <- rowSums(metrics_data[, which(names(metrics_data) %in% c("pearson", "spearman"))])/length(which(names(metrics_data) %in% c("pearson", "spearman")))}

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
