# Function to create the spatial cross validation dataset

#' Title
#'
#' @param dats 
#' @param n.folds 
#'
#' @return
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
    # train$set <- rep("fitting", nrow(train))
    # test$set <- rep("validation", nrow(test))
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