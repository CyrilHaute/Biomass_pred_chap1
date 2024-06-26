#' Title spatialrf_function
#' 
#' This function fit a spatial random forest using the R package `SpatialML` with a k fold spatial cross validation procedure
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
# base_dir = eco_base_dir

spatialrf_function <- function(biomass, 
                               covariates,
                               base_dir){
  
  fmla <<- as.formula(paste0("biomass ~ ", paste0(colnames(covariates)[!colnames(covariates) %in% "survey_id"], collapse = " + ")))
  
  sp <- biomass
  sp_name <- colnames(sp)[5]
  colnames(sp)[colnames(sp) %in% sp_name] <- "biomass"
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
    training <<- biomass_scv[[j]]$train[,c("survey_id", "biomass")]
    
    # select the jth species from the validation set
    testing <- biomass_scv[[j]]$test[,c("survey_id", "biomass")]
    
    training <<- training |> 
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
    training <<- dplyr::inner_join(training, train_covariates, by = "survey_id")
    
    # add covariates
    testing <- dplyr::inner_join(testing, test_covariates, by = "survey_id")
    
    coords <- training |> 
      dplyr::select(X, Y) |> 
      as.data.frame()
    
    ### FITTING MODELS 
    # fit the spatial random forests
    
    model_fit <- SpatialML::grf(formula = fmla,
                                dframe = training,
                                bw = 100,
                                kernel = "fixed",
                                coords = coords,
                                ntree = 100,
                                geo.weighted = FALSE)
    
    validation_predict  <- SpatialML::predict.grf(object = model_fit,
                                                  new.data = testing,
                                                  x.var.name = "X",
                                                  y.var.name = "Y",
                                                  local.w = 0.5,
                                                  global.w = 0.5)
    
    validation_predict <- data.frame(survey_id = testing$survey_id,
                                     validation_predict = validation_predict)
    
    validation_observed <- testing[,c("survey_id", "biomass")]
    
    validation_obs_prd <- validation_predict |>
      dplyr::inner_join(validation_observed, multiple = "first") |> 
      dplyr::mutate(species_name = sp_name,
                    cv = j,
                    model = "sprf")
    
    colnames(validation_obs_prd)[colnames(validation_obs_prd) == "biomass"] <- "validation_observed"
    
    validation_obs_prd$validation_predict <- 10^validation_obs_prd$validation_predict - 1
    validation_obs_prd$validation_observed <- 10^validation_obs_prd$validation_observed - 1
    
    return(validation_obs_prd)
    
  }, mc.cores = parallel::detectCores() - 1)
  
  cv_j_bind <- do.call(rbind, cv_j)
  
  model_dir <- "sprf/"
  
  dir.create(base_dir)
  dir.create(paste0(base_dir, model_dir))
  
  save(cv_j_bind, file = paste0(base_dir, model_dir, stringr::str_replace_all(sp_name, " ", "_"), "_extracted_predictions.RData"))
  
}
