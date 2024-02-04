
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

scv_function <- function(dats, 
                         n.folds){
  
  # flexible object for storing folds
  folds <- list()
  
  fold.size <- nrow(dats)/n.folds
  
  # all obs are in
  remain <- 1:nrow(dats)
  
  for(i in 1:n.folds) {
    
    # randomly sample “fold_size” from the ‘remaining observations’
    select <- sample(remain, fold.size, replace = FALSE)
    
    # store indices
    folds[[i]] <- select
    
    if (i == n.folds){
      
      folds[[i]] <- remain
      
    }
    
    # update remaining indices to reflect what was taken out
    remain <- setdiff(remain, select)
    
    remain
  }
  
  train_test <- list()
  
  for(i in 1:n.folds) {
    
    # fold i
    # unpack into a vector
    indis <- folds[[i]]
    
    # split into train and test sets
    train <- dats[-indis,]
    test <- dats[indis,]
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