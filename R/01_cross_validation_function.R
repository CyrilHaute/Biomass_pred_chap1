
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

# dats = sp
# n.folds = 10

scv_function <- function(dats, 
                         n.folds){
  
  # flexible object for storing folds
  folds <- list()
  
  fold.size <- nrow(dats)/n.folds
  
  # all obs are in
  remain <- 1:nrow(dats)
  
  for(i in 1:n.folds) {
    
    # randomly sample “fold_size” from the ‘remaining observations’
    select.val <- sample(remain, fold.size, replace = FALSE)
    
    # store indices
    folds[[i]] <- select.val
    
    if (i == n.folds){
      
      folds[[i]] <-  remain
      
    }
    
    # update remaining indices to reflect what was taken out
    remain <- setdiff(remain, select.val)
    
  }
  
  train_test_val <- list()
  
  for(i in 1:n.folds) {
    
    # fold i
    # unpack into a vector
    indis <- folds[[i]]
    
    train <- dats[-indis,]
    
    test <- dats[indis,]
    
    train_test_val[[i]] <- list(train, test)
    names(train_test_val[[i]]) <- c("train", "test")
    
    # delete from the train set, transects that appears in the same site than the point in the validation set
    train_test_val[[i]]$train <- train_test_val[[i]]$train |>
      dplyr::filter(!site_code %in% unique(train_test_val[[i]]$test$site_code))
    
    train_test_val[[i]]$train <- train_test_val[[i]]$train |> 
      dplyr::select(-site_code)
    
    train_test_val[[i]]$test <- train_test_val[[i]]$test |> 
      dplyr::select(-site_code)
    
    
  }
  
  return(train_test_val)
  
}
