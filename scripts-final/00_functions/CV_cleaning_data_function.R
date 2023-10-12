# Function to create the spatial cross validation dataset

CV <- function(dats, n.folds){
  folds <- list() # flexible object for storing folds
  fold.size <- nrow(dats)/n.folds
  remain <- 1:nrow(dats) # all obs are in
  
  for (i in 1:n.folds){
    select <- sample(remain, fold.size, replace = FALSE)
    #randomly sample “fold_size” from the ‘remaining observations’
    
    folds[[i]] <- select # store indices
    
    if (i == n.folds){
      folds[[i]] <- remain
    }
    
    #update remaining indices to reflect what was taken out
    remain <- setdiff(remain, select)
    remain
  }
  
  train_test <- list()
  
  for (i in 1:n.folds){
    # fold i
    indis <- folds[[i]] #unpack into a vector
    train <- dats[-indis,] #split into train and test sets
    test <- dats[indis,]
    train$set <- rep("fitting", nrow(train))
    test$set <- rep("validation", nrow(test))
    train_test[[i]] <- list(train, test)
    names(train_test[[i]]) <- c("fitting", "validation")
    train_test[[i]]$fitting <- train_test[[i]]$fitting %>% #delete from the train set transects that appears in the same site than the point in the validation set
      filter(!SiteCode %in% train_test[[i]]$validation$SiteCode)
    train_test[[i]]$fitting <- train_test[[i]]$fitting %>% dplyr::select(-c(SiteCode, set))
    train_test[[i]]$validation <- train_test[[i]]$validation %>% dplyr::select(-c(SiteCode, set))
  }
  
  return(train_test)
  
}