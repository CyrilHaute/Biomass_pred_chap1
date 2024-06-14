# function to fit spatial Random Forest and assess covariates relative importance

biomass = rls_biomass
covariates = rls_covariates
species_name = colnames(rls_biomass)[!colnames(rls_biomass) %in% c("survey_id", "latitude", "longitude", "site_code")]
base_dir_cont = base_dir

spatialrf_function_cont <- function(biomass,
                                    covariates,
                                    species_name,
                                    base_dir_cont){
  
  fmla <<- as.formula(paste0("biomass ~ ", paste0(colnames(covariates)[!colnames(covariates) %in% "survey_id"], collapse = " + ")))
  
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
    
    coords <<- sp |>
      dplyr::rename(Y = latitude,
                    X = longitude) |> 
      dplyr::select(Y, X) |>
      as.data.frame()
    
    sp <<- sp |>
      dplyr::rename(Y = latitude,
                    X = longitude)
    
    ### FITTING MODELS 
    # fit the spatial random forests
    
    model_fit <- SpatialML::grf(formula = fmla,
                                dframe = sp,
                                bw = 100,
                                kernel = "fixed",
                                coords = coords,
                                ntree = 100,
                                geo.weighted = FALSE)
    
    covnames_new_new <- colnames(covariates)[which(colnames(covariates) %in% "survey_id" == FALSE)]
    
    # Use the package DALEX to assess covariates relative importance
    # First create an explain object (a representation of your model, depend on the structure of the algorithm used)
    global_explainer_sprf <- DALEX::explain(model = model_fit[[1]],
                                            data = sp[covnames_new_new],
                                            y = sp[,"biomass"],
                                            label = "ranger")
    
    # Compute a 10-permutation-based value of the RMSE for all explanatory variables  
    global_vip.10_sprf <- DALEX::model_parts(explainer = global_explainer_sprf,
                                             loss_function = DALEX::loss_root_mean_square, # Here we used the RMSE as our loss function
                                             B = 10, # Number of permutation
                                             type = "difference")
    
    # From the model_parts function you get 10 RMSE values for each covariates. 
    # Take the mean and assess the standard-deviation of the RMSE for each covariates to assess the error of the permutation method
    global_vip.10_sprf <- global_vip.10_sprf |> 
      dplyr::group_by(variable) |> 
      dplyr::summarise(global_dropout_loss = mean(dropout_loss),
                       global_sd_dropout_loss = sd(dropout_loss))
    
    global_vip.10_sprf <- global_vip.10_sprf |> 
      dplyr::filter(!variable %in% c("_baseline_", "_full_model_"))
    
    local_importance <- pbmcapply::pbmclapply(1:length(model_fit$Forests), function(j) {
      
      explainer_sprf <- DALEX::explain(model = model_fit$Forests[[j]],
                                       data = sp[covnames_new_new],
                                       y = sp[,"biomass"],
                                       label = "ranger")
      
      vip.10_sprf <- DALEX::model_parts(explainer = explainer_sprf,
                                        loss_function = DALEX::loss_root_mean_square,
                                        B = 10,
                                        type = "difference")
      
      vip.10_sprf <- vip.10_sprf |> 
        dplyr::group_by(variable) |> 
        dplyr::summarise(local_dropout_loss = mean(dropout_loss),
                         local_sd_dropout_loss = sd(dropout_loss)) |> 
        dplyr::inner_join(global_vip.10_sprf)
      
      vip.10_sprf <- vip.10_sprf |> 
        dplyr::mutate(dropout_loss = global_dropout_loss * 0.5 + local_dropout_loss * 0.5,
                      sd_dropout_loss = global_sd_dropout_loss * 0.5 + local_sd_dropout_loss * 0.5) |> 
        dplyr::select(variable, dropout_loss, sd_dropout_loss)
      
      
    }, mc.cores = parallel::detectCores() - 1)
    
    local_importance_mean <- average_dfs(local_importance)
    
  }, mc.cores = 1)
  
  extracted_contributions <- dplyr::tibble(species_name = species_name, 
                                           fitted_model = "sprf", 
                                           # estimate contribution
                                           contributions_and_sd = sp_i)
  
  # save contribution output in same file structure
  
  model_dir <- "sprf"
  
  dir.create(base_dir_cont)
  
  save(extracted_contributions, file = paste0(base_dir_cont, model_dir, "_extracted_contributions.RData"))
  
  rm(list=ls())
  gc()
  
}
