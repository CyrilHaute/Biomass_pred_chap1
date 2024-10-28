# function to fit boosted regression tree

# biomass = biomass_scv
# covariates = rls_covariates
# species_name = colnames(biomass_scv[[1]]$fitting)[!colnames(biomass_scv[[1]]$fitting) %in% c("survey_id", "latitude", "longitude")]
# n.cores = 1
# base_dir = base_dir

#' Title brt_function
#' 
#' This function fit a gradient boosting machine model using the R package `gbm` with a k fold spatial cross validation procedure
#'
#' @param biomass a list in which each elements is a fold of the spatial cross validation procededure. Each fold is split into two subset, the first one named "fitting" to
#' train the model and the second one named "validation" to test the model
#' @param covariates a datagrame containg all covariates to fit the model
#' @param species_name a vector containg the name of all species contain in @param biomass
#' @param base_dir the path to save the data
#'
#' @return a dataframe with as many row as the length of @param species_name . Each row is a species with its biomass observation and prediction from each cross validation fold
#' @export
#'
#' @examples

# biomass = rls_biomass_i
# covariates = rls_covariates
# species_name = species_name[i]
# base_dir = base_dir

brt_function <- function(biomass, 
                         covariates, 
                         species_name,
                         n.cores,
                         base_dir){
  
  brt_formula <- as.formula(paste0("biomass ~ ", paste0(colnames(covariates)[!colnames(covariates) %in% "survey_id"], collapse = " + ")))
  
  # First, select only zero within the living area of the species considered
  
  sp <- biomass
  colnames(sp)[colnames(sp) %in% species_name] <- "biomass" # Rename the species name for consistency
  sp$biomass <- as.numeric(sp$biomass)
  sp$latitude_arr <- round(sp$latitude, 1) # Do so by considering only 0 where presences are recorded within degree of latitude
  
  pres <- sp |> 
    dplyr::filter(biomass != 0)
  abs <- sp |> 
    dplyr::filter(biomass == 0) |> 
    dplyr::filter(latitude_arr %in% unique(pres$latitude_arr))
  
  sp <- rbind(pres, abs) |> 
    dplyr::select(! latitude_arr) # merge presences and selected absences
  
  # Even with the absences selection above, it can still result in to much absences in the dataset resulting in zero inflation.
  # Thus, select randomly as much absences as presences
  
  biomass_only <- sp[which(sp[, "biomass"] > 0),]
  
  n_subsample <- nrow(sp[which(sp[, "biomass"] > 0),])
  
  absence <- sp[which(sp[, "biomass"] == 0),]
  
  if(nrow(absence) > 0) {
    
    replacement <- ifelse(length(which(sp[, "biomass"] == 0)) < n_subsample, T, F)
    
    absence <- absence[sample(which(absence[, "biomass"] == 0), n_subsample, replace = replacement),]
    
  }
  
  # Combine absences and presences
  sp <- rbind(biomass_only, absence)
  
  covariates_sp <- covariates |> 
    dplyr::filter(survey_id %in% sp$survey_id) |> 
    noise_function(avoid = c("survey_id", "effectiveness"),
                   limit = 6)
  
  covariates_sp[!colnames(covariates_sp) %in% c("survey_id", "effectiveness")] <- scale(covariates_sp[!colnames(covariates_sp) %in% c("survey_id", "effectiveness")], center = TRUE, scale = TRUE)
  
  # Create spatial k-fold cross-validation dataset, here with 5 fold, each fold being splited in 80% for training and 20% for testing. The spatial compenent can resulting in less than 80% of data in the training set
  biomass_scv <- scv_function(sp,
                              5)
  
  cv_j <- pbmcapply::pbmclapply(1:length(biomass_scv), function(j) {
    
    # select the jth species from the training set
    training <- biomass_scv[[j]]$train[,c("survey_id", "biomass")]
    
    # select the jth species from the validation set
    testing <- biomass_scv[[j]]$test[,c("survey_id", "biomass")]
    
    training <- training |> 
      dplyr::inner_join(rls_surveys[,colnames(rls_surveys) %in% c("survey_id", "latitude", "longitude")]) |> 
      dplyr::rename(X = longitude, Y = latitude)
    
    testing <- testing |> 
      dplyr::inner_join(rls_surveys[,colnames(rls_surveys) %in% c("survey_id", "latitude", "longitude")]) |> 
      dplyr::rename(X = longitude, Y = latitude)
    
    # As some covariates are at the country level, it means you can have very few or even only one value for these covariates
    # Check for the number of values in each covariates and add noise if < 6 values
    
    train_covariates <- covariates_sp |> 
      dplyr::filter(survey_id %in% training$survey_id) |> 
      noise_function(avoid = c("survey_id", "effectiveness"),
                     limit = 6)
    test_covariates <- covariates_sp |> 
      dplyr::filter(survey_id %in% testing$survey_id) |> 
      noise_function(avoid = c("survey_id", "effectiveness"),
                     limit = 6)
    
    # add covariates
    training <- dplyr::inner_join(training, train_covariates, by = "survey_id")
    
    # add covariates
    testing <- dplyr::inner_join(testing, test_covariates, by = "survey_id")
    
    ### FITTING BOOSTED REGRESSION TREES
    
    model_fit <- gbm::gbm(formula = brt_formula,
                          data = training, 
                          distribution = "gaussian", 
                          n.trees = 10000,
                          interaction.depth = 3, 
                          shrinkage = 0.001,
                          bag.fraction = 0.5,
                          n.cores = 1)
    
    validation_predict  <- predict(object = model_fit, testing)
    
    validation_predict <- data.frame(survey_id = testing$survey_id,
                                     validation_predict = validation_predict)
    
    validation_observed <- testing[,c("survey_id", "biomass")]
    
    validation_obs_prd <- validation_predict |>
      dplyr::inner_join(validation_observed, multiple = "first") |> 
      dplyr::mutate(species_name = species_name,
                    cv = j,
                    model = "gbm")
    
    colnames(validation_obs_prd)[colnames(validation_obs_prd) == "biomass"] <- "validation_observed"
    
    validation_obs_prd$validation_predict <- 10^validation_obs_prd$validation_predict - 1
    validation_obs_prd$validation_observed <- 10^validation_obs_prd$validation_observed - 1
    
    return(validation_obs_prd)
    
  }, mc.cores = 1)
  
  cv_j_bind <- do.call(rbind, cv_j)
  
  # save prediciton output in same file structure
  
  model_dir <- "gbm"
  
  dir.create(base_dir)
  
  save(cv_j_bind, file = paste0(base_dir, stringr::str_replace_all(species_name, " ", "_"), ".Rdata"))
  
}
