
load_realm_cont_function <- function(files_path) {
  
  realm <- list.files(files_path, full.names = TRUE)
  realm_small <- list.files(files_path, full.names = FALSE)
  realm <- lapply(1:length(realm), function(i) {
    
    realm_i <- list.files(realm[i], full.names = TRUE)
    
    realm_j <- lapply(1:length(realm_i), function(j) {
      
      load(realm_i[j])
      assign(paste0(j), extracted_contributions)
      
    })
    realm_j <- do.call(rbind, realm_j)
    realm_j$realm <- realm_small[i]
    
    return(realm_j)
    
  })

}

load_realm_partial_function <- function(files_path) {
  
  realm <- list.files(files_path, full.names = TRUE)
  realm_small <- list.files(files_path, full.names = FALSE)
  realm <- lapply(1:length(realm), function(i) {
    
    realm_i <- list.files(realm[i], full.names = TRUE)
    
    realm_j <- lapply(1:length(realm_i), function(j) {
      
      load(realm_i[j])
      assign(paste0(j), partial_plot)
      
    })

    return(realm_j)
    
  })
  
}
