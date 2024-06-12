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

# biomass = rls_biomass
# covariates = rls_covariates
# species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")]
# base_dir = base_dir

brt_function <- function(biomass, 
                         covariates, 
                         species_name,
                         n.cores,
                         base_dir){
  
  brt_formula <- as.formula(paste0("biomass ~ ", paste0(colnames(covariates)[!colnames(covariates) %in% "survey_id"], collapse = " + ")))
  
  sp_i <- pbmcapply::pbmclapply(1:length(species_name), function(i) {
    
    sp <- biomass[colnames(biomass) %in% c("survey_id", "latitude", "longitude", "site_code", species_name[i])]
    colnames(sp)[colnames(sp) %in% species_name[i]] <- "biomass"
    sp$latitude_arr <- round(sp$latitude, 1)
    
    pres <- sp |> 
      dplyr::filter(biomass != 0)
    abs <- sp |> 
      dplyr::filter(biomass == 0) |> 
      dplyr::filter(latitude_arr %in% unique(pres$latitude_arr))
    
    sp <- rbind(pres, abs) |> 
      dplyr::select(! latitude_arr)
    
    biomass_only <- sp[which(sp[, "biomass"] > 0),]
    
    n_subsample <- nrow(sp[which(sp[, "biomass"] > 0),])
    
    absence <- sp[which(sp[, "biomass"] == 0),]
    
    if(nrow(absence) > 0) {
      
      replacement <- ifelse(length(which(sp[, "biomass"] == 0)) < n_subsample, T, F)
      
      absence <- absence[sample(which(absence[, "biomass"] == 0), n_subsample, replace = replacement),]
      
    }
    
    # combine absence and presence
    sp <- rbind(biomass_only, absence)
    
    sp$biomass <- as.numeric(sp$biomass)
    
    biomass_scv <- scv_function(sp,
                                10)
    
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
      
      train_covariates <- covariates |> 
        dplyr::filter(survey_id %in% training$survey_id) |> 
        noise_function(avoid = c("survey_id", "effectiveness"),
                       limit = 6,
                       size = 0.01)
      test_covariates <- covariates |> 
        dplyr::filter(survey_id %in% testing$survey_id) |> 
        noise_function(avoid = c("survey_id", "effectiveness"),
                       limit = 6,
                       size = 0.01)
      
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
                            bag.fraction = 0.8,
                            n.cores = 1)
      
      validation_predict  <- predict(object = model_fit, testing)
      
      validation_predict <- data.frame(survey_id = testing$survey_id,
                                       validation_predict = validation_predict)
      
      validation_observed <- testing[,c("survey_id", "biomass")]
      
      validation_obs_prd <- validation_predict |>
        dplyr::inner_join(validation_observed, multiple = "first") |> 
        dplyr::mutate(species_name = species_name[i],
                      cv = j,
                      model = "gbm")
      
      colnames(validation_obs_prd)[colnames(validation_obs_prd) == "biomass"] <- "validation_observed"
      
      validation_obs_prd$validation_predict <- 10^validation_obs_prd$validation_predict - 1
      validation_obs_prd$validation_observed <- 10^validation_obs_prd$validation_observed - 1
      
      return(validation_obs_prd)
      
    }, mc.cores = 1)
    
    cv_j_bind <- do.call(rbind, cv_j)
    
  }, mc.cores = parallel::detectCores() - 1)
  
  # save prediciton output in same file structure
  
  model_dir <- "gbm"
  
  dir.create(base_dir, recursive = T)
  
  save(sp_i, file = paste0(base_dir, model_dir, "_extracted_predictions.RData"))
  
  rm(list=ls())
  gc()
  
}
