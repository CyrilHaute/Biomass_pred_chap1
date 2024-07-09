# function to fit gam and assess covariates relative importance

# biomass = realm
# covariates = rls_covariates
# species_name = new_species_name
# base_dir_cont = eco_base_dir

gam_function_cont <- function(biomass, 
                              covariates,
                              species_name,
                              base_dir_cont){
  
  # response variable name
  response <- "biomass ~ "
  
  # rename covariates
  covnames_new <- names(covariates)
  covnames_new <- covnames_new[-which(covnames_new %in% c("survey_id", "effectiveness"))]
  covnames_new <- c(covnames_new, "factor(effectiveness)")
  covnames_new_bis <- covnames_new[-which(covnames_new %in% c("factor(effectiveness)"))]
  
  # create formula with new covariate names
  covnames_splines <- paste0(rep("s(", length(covnames_new)),  covnames_new, rep(", k = 3)", length(covnames_new)))
  covnames_splines <- covnames_splines[-which(covnames_splines %in% c("s(factor(effectiveness), k = 3)"))]
  covnames_splines <- c(covnames_splines, "factor(effectiveness)")
  covnames_splines2 <- paste0(rep("s(", length(covnames_new_bis)),  covnames_new_bis, rep(", k = 3)", length(covnames_new_bis)))
  
  covnames_combined <- paste0(covnames_splines, collapse = " + ")
  covnames_combined2 <- paste0(covnames_splines2, collapse = " + ")
  
  # find the full model formula
  model_formula <- as.formula(paste0(response, covnames_combined))
  model_formula2 <- as.formula(paste0(response, covnames_combined2))

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
    
    if(length(unique(sp$effectiveness)) == 1){
      
      model_fit <- mgcv::gam(model_formula2, data = sp, family = gaussian, select = FALSE, method = 'ML')
      
    }else{
      
      model_fit <- mgcv::gam(model_formula, data = sp, family = gaussian, select = FALSE, method = 'ML')
      
    }
    
    covnames_new_new <- colnames(covariates)[which(colnames(covariates) %in% "survey_id" == FALSE)]

    # Use the package DALEX to assess covariates relative importance
    # First create an explain object (a representation of your model, depend on the structure of the algorithm used)
    explainer_gam <- DALEX::explain(model = model_fit, 
                                    data = sp[covnames_new_new], 
                                    y = sp[,"biomass"],
                                    label = "gam")
      
    # Compute a 25-permutation-based value of the RMSE for all explanatory variables 
    vip.25_gam <- DALEX::model_parts(explainer = explainer_gam,
                                     loss_function = DALEX::loss_root_mean_square, # Here we used the RMSE as our loss function
                                     B = 25, # Number of permutation
                                     type = "difference")
      
    # From the model_parts function you get 25 RMSE values for each covariates. 
    # Take the mean and assess the standard-deviation of the RMSE for each covariates to assess the error of the permutation method
    vip.25_gam <- vip.25_gam |> 
      dplyr::group_by(variable) |> 
      dplyr::summarise(Dropout_loss = mean(dropout_loss),
                       sd_dropout_loss = sd(dropout_loss))
      
    vip.25_gam <- vip.25_gam |> 
      dplyr::filter(!variable %in% c("_baseline_", "_full_model_"))

    }, mc.cores = parallel::detectCores() - 1)
  
  extracted_contributions <- dplyr::tibble(species_name = species_name, 
                                           fitted_model = "gam", 
                                           # estimate contribution
                                           contributions_and_sd = sp_i)
  
  # save contribution output in same file structure
  
  model_dir <- "gam"
  
  dir.create(base_dir_cont)
  
  save(extracted_contributions, file = paste0(base_dir_cont, model_dir, "_extracted_contributions.RData"))
  
  rm(list=ls())
  gc()
  
}
