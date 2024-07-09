noise_function <- function(data,
                           avoid,
                           limit,
                           size){
  
  n_values <- lapply(1:ncol(data[,!colnames(data) %in% avoid]), function(i) {unique(data[,!colnames(data) %in% avoid][,i])})
  
  names(n_values) <- colnames(data[,!colnames(data) %in% avoid])
  
  little_cov <- names(n_values[which(sapply(1:length(n_values), function(i) {nrow(n_values[[i]])}) <= limit)])
  
  if(sjmisc::is_empty(little_cov) == TRUE){
    
    data <- data
    
  }else{
    
    n_cov <- which(names(data) %in% little_cov)
    
    noise <- lapply(1:length(n_cov), function(i) {
      
      abs(rnorm(nrow(data), size, size))
      
    })
    
    if(length(noise) == 1){
      
      data[,n_cov] <- data[,n_cov] + unlist(noise)
      
    }else{
      
      data[,n_cov] <- data[,n_cov] + noise
      
    }
    
  }
  
  return(data)
  
}
