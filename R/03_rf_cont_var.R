# function to fit Random Forest and assess covariates relative importance

rf_function_cont <- function(biomass, 
                             covariates, 
                             species_name, 
                             base_dir_cont){
  
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
    noise_function(avoid = c("survey_id", "protection_status2"),
                   limit = 6)
  
  covariates_sp[!colnames(covariates_sp) %in% c("survey_id", "protection_status2")] <- scale(covariates_sp[!colnames(covariates_sp) %in% c("survey_id", "protection_status2")], center = TRUE, scale = TRUE)
  
  # add covariates
  sp <- dplyr::inner_join(sp, covariates_sp, by = "survey_id")
  
  # Fit the random forest model
  
  model_fit <- ranger::ranger(x = sp[!colnames(sp) %in% c("survey_id", "biomass", "longitude", "latitude", "site_code")],
                              y = unlist(sp[colnames(sp) %in% "biomass"]),
                              num.trees = ncol(sp[!colnames(sp) %in% c("survey_id", "biomass", "longitude", "latitude", "site_code")]) * 10)
  
  # Use the package DALEX to assess covariates relative importance
  # First create an explain object (a representation of your model, depend on the structure of the algorithm used)
  explainer_rf <- DALEX::explain(model = model_fit, 
                                 data = sp[!colnames(sp) %in% c("survey_id", "latitude", "longitude", "site_code", "biomass")],
                                 y = unlist(sp[colnames(sp) %in% "biomass"]),
                                 label = "ranger")
  
  # Compute a 10-permutation-based value of the RMSE for all explanatory variables
  vip.10_rf <- DALEX::model_parts(explainer = explainer_rf,
                                  loss_function = DALEX::loss_root_mean_square, # Here we used the RMSE as our loss function
                                  B = 10, # Number of permutation
                                  type = "difference")
  
  # From the model_parts function you get 10 RMSE values for each covariates. 
  # Take the mean and assess the standard-deviation of the RMSE for each covariates to assess the error of the permutation method
  vip.10_rf <- vip.10_rf |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(Dropout_loss = mean(dropout_loss),
                     sd_dropout_loss = sd(dropout_loss))
  
  vip.10_rf <- vip.10_rf |> 
    dplyr::filter(!variable %in% c("_baseline_", "_full_model_"))
  
  
  extracted_contributions <- dplyr::tibble(species_name = species_name, 
                                           fitted_model = "rf", 
                                           # estimate contribution
                                           vip.10_rf)
  
  # save contribution output in same file structure
  
  dir.create(base_dir_cont)
  
  save(extracted_contributions, file = paste0(base_dir_cont, stringr::str_replace_all(species_name, " ", "_"), ".Rdata"))
  
}
