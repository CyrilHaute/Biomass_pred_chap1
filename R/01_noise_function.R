# data = scenario_l
# avoid = c("survey_id", "scenario", "effectiveness", "biomass", "y", "x", "year")
# limit = 6

noise_function <- function(data,
                           avoid,
                           limit){

  n_values <- lapply(1:ncol(data[,!colnames(data) %in% avoid]), function(i) {unique(data[,!colnames(data) %in% avoid][,i])})
  
  names(n_values) <- colnames(data[,!colnames(data) %in% avoid])
  
  little_cov <- names(n_values[which(sapply(1:length(n_values), function(i) {nrow(n_values[[i]])}) <= limit)])
  
  if(sjmisc::is_empty(little_cov) == TRUE){
    
    data <- data
    
  }else{
    
    mean_value <- sapply(1:length(little_cov), function(i) { mean(unlist(data[,little_cov[i]]))})
    
    n_cov <- which(names(data) %in% little_cov)
    
    noise <- lapply(1:length(n_cov), function(i) {
      
      if(mean_value[i] != 0){
      
      abs(rnorm(nrow(data), abs(mean_value[i]) * 0.001, abs(mean_value[i]) * 0.001))
        
      }else{
        
        abs(rnorm(nrow(data), 0.001, 0.001))
        
      }
      
    })
    
    if(length(noise) == 1){
      
      data[,n_cov] <- data[,n_cov] + unlist(noise)
      
    }else{
      
      data[,n_cov] <- data[,n_cov] + noise
      
    }
    
  }
  
  print(little_cov)
  
  return(data)

}
