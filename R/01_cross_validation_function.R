
#' Title : scv_function
#' 
#' This function create a k fold spatial cross validation.
#'
#' @param dats a dataframe that is to be split into k folds
#' @param n.folds a numeric corresponding to the number of k fold of the cross validation
#'
#' @return a list with k element. Each element of the list is a fold of the cross validation. In each element, you have two dataframes, one is for model training (named "fitting")
#' the second is for model testing (named "validation")
#' @export
#'
#' @examples

dats = rls_biomass
n.folds = 10

scv_function <- function(dats, 
                         n.folds){
  
  # flexible object for storing folds
  folds <- list()
  
  # Store levels of categorical variables
  cat_levels <- levels(dats$effectiveness)

  test <- sapply(1:length(cat_levels), function(i) { 
    
    length(which(dats$effectiveness %in% cat_levels[i] == TRUE))/n.folds
    
    })
  names(test) <- cat_levels
  
  positive_indices <- which.min(unlist(lapply(dats[, -c(1:4, 597:620)], function(col) length(which(col > 0))))) + 4
  
  fold.size <- nrow(dats)/n.folds
  fold.size.pos <- nrow(dats[which(dats[,positive_indices] > 0),])/n.folds
  fold.size.zero <- nrow(dats[which(dats[,positive_indices] == 0),])/n.folds
  
  # all obs are in
  remain.pos <- which(dats[,positive_indices] > 0)
  remain.zero <- which(dats[,positive_indices] == 0)

  for(i in 1:n.folds) {
    
    # randomly sample “fold_size” from the ‘remaining observations’
    select.pos <- sample(remain.pos, fold.size.pos, replace = FALSE)
    select.zero <- sample(remain.zero, fold.size.zero, replace = FALSE)
    
    # store indices
    folds[[i]] <- c(select.pos, select.zero)
    
    if (i == n.folds){
      
      folds[[i]] <- c(remain.pos, remain.zero)
      
    }
    
    # update remaining indices to reflect what was taken out
    remain.pos <- setdiff(remain.pos, select.pos)
    remain.zero <- setdiff(remain.zero, select.zero)
    
  }
  
  train_test <- list()
  
  for(i in 1:n.folds) {
    
    # fold i
    # unpack into a vector
    indis <- folds[[i]]
    
    # split into train and test sets
    train <- dats[-indis,]
    test <- dats[indis,]
    
    # Restore levels of categorical variables
    for (j in 1:length(cat_levels)) {
      
      if (!is.null(cat_levels[[j]])) {
        
        train[, j] <- factor(train[, j], levels = cat_levels[[j]])
        test[, j] <- factor(test[, j], levels = cat_levels[[j]])
        
      }
      
    }
    
    train_test[[i]] <- list(train, test)
    names(train_test[[i]]) <- c("fitting", "validation")
    
    # delete from the train set, transects that appears in the same site than the point in the validation set
    train_test[[i]]$fitting <- train_test[[i]]$fitting |> 
      dplyr::filter(!site_code %in% train_test[[i]]$validation$site_code)
    
    train_test[[i]]$fitting <- train_test[[i]]$fitting |> 
      dplyr::select(-site_code)
    
    train_test[[i]]$validation <- train_test[[i]]$validation |> 
      dplyr::select(-site_code)
    
  }
  
  return(train_test)
  
}