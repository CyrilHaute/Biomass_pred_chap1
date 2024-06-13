# function to fit Random Forest and assess covariates relative importance

biomass = rls_biomass
covariates = rls_covariates
species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")]
base_dir_cont = base_dir

rf_function_cont <- function(biomass, 
                             covariates, 
                             species_name, 
                             base_dir_cont){

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
    
    covariates <- covariates |> 
      dplyr::filter(survey_id %in% sp$survey_id) |> 
      noise_function(avoid = c("survey_id", "effectiveness"),
                     limit = 6,
                     size = 0.01)
    
    # add covariates
    sp <- dplyr::inner_join(sp, covariates, by = "survey_id")
    
    # Fit the model
    
    model_fit <- ranger::ranger(x = sp[!colnames(sp) %in% c("survey_id", "latitude", "longitude", "site_code", "biomass")],
                                y = unlist(sp[colnames(sp) %in% "biomass"]),
                                num.trees = 100,
                                mtry = round(length(colnames(sp)[!colnames(sp) %in% c("survey_id", "biomass")]) / 3))
  
    # Use the package DALEX to assess covariates relative importance
    # First create an explain object (a representation of your model, depend on the structure of the algorithm used)
    explainer_rf <- DALEX::explain(model = model_fit, 
                                   data = sp[!colnames(sp) %in% c("survey_id", "latitude", "longitude", "site_code", "biomass")],
                                   y = unlist(sp[colnames(sp) %in% "biomass"]),
                                   label = "ranger")
      
    # Compute a 25-permutation-based value of the RMSE for all explanatory variables
    vip.25_rf <- DALEX::model_parts(explainer = explainer_rf,
                                    loss_function = DALEX::loss_root_mean_square, # Here we used the RMSE as our loss function
                                    B = 25, # Number of permutation
                                    type = "difference")
      
    # From the model_parts function you get 25 RMSE values for each covariates. 
    # Take the mean and assess the standard-deviation of the RMSE for each covariates to assess the error of the permutation method
    vip.25_rf <- vip.25_rf |> 
      dplyr::group_by(variable) |> 
      dplyr::summarise(Dropout_loss = mean(dropout_loss),
                       sd_dropout_loss = sd(dropout_loss))
      
    vip.25_rf <- vip.25_rf |> 
      dplyr::filter(!variable %in% c("_baseline_", "_full_model_"))

    }, mc.cores = parallel::detectCores() - 1)
  
  extracted_contributions <- dplyr::tibble(species_name = species_name, 
                                           fitted_model = "rf", 
                                           # estimate contribution
                                           contributions_and_sd = sp_i)
  
  # save contribution output in same file structure
  
  model_dir <- "rf"
  
  dir.create(base_dir_cont)
  
  save(extracted_contributions, file = paste0(base_dir_cont, model_dir, "_extracted_contributions.RData"))
  
  rm(list=ls())
  gc()
  
}
