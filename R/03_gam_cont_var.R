# function to fit gam and assess covariates relative importance

# biomass = realm_j
# covariates = rls_covariates
# species_name = new_species_name[j]
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
  
  # add covariates
  sp <- dplyr::inner_join(sp, covariates_sp, by = "survey_id")
  
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
  
  # Compute a 10-permutation-based value of the RMSE for all explanatory variables 
  vip.10_gam <- DALEX::model_parts(explainer = explainer_gam,
                                   loss_function = DALEX::loss_root_mean_square, # Here we used the RMSE as our loss function
                                   B = 10, # Number of permutation
                                   type = "difference")
  
  # From the model_parts function you get 10 RMSE values for each covariates. 
  # Take the mean and assess the standard-deviation of the RMSE for each covariates to assess the error of the permutation method
  vip.10_gam <- vip.10_gam |> 
    dplyr::group_by(variable) |> 
    dplyr::summarise(Dropout_loss = mean(dropout_loss),
                     sd_dropout_loss = sd(dropout_loss))
  
  vip.10_gam <- vip.10_gam |> 
    dplyr::filter(!variable %in% c("_baseline_", "_full_model_"))
  
  extracted_contributions <- dplyr::tibble(species_name = species_name, 
                                           fitted_model = "gam", 
                                           # estimate contribution
                                           vip.10_gam)
  
  # save contribution output in same file structure
  
  dir.create(base_dir_cont)
  
  save(extracted_contributions, file = paste0(base_dir_cont, stringr::str_replace_all(species_name, " ", "_"), ".Rdata"))
  
}
