# This function create a k fold spatial cross validation

scv_function <- function(dats, 
                         n.folds){
  
  folds <- list()
  
  fold.size <- nrow(dats)/n.folds
  
  # all obs are in
  remain <- 1:nrow(dats)
  
  for(i in 1:n.folds) {
    
    select.val <- sample(remain, fold.size, replace = FALSE)
    
    folds[[i]] <- select.val
    
    if (i == n.folds){
      
      folds[[i]] <-  remain
      
    }
    
    remain <- setdiff(remain, select.val)
    
  }
  
  train_test_val <- list()
  
  for(i in 1:n.folds) {
    
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
