#' Title rf_function
#' 
#' This function fit a random forest using the R package `randomForest` with a k fold spatial cross validation procedure
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

# biomass = realm_j
# covariates = rls_covariates
# species_name = new_species_name[j]
# base_dir = eco_base_dir

partial_function_spamm <- function(biomass, 
                                   covariates, 
                                   species_name,
                                   base_dir){
  
  # response variable name
  response <- "biomass ~ "
  
  # rename covariates
  covnames_new <- names(covariates)
  covnames_new <- covnames_new[-which(covnames_new %in% c("survey_id", "protection_status2"))]
  covnames_new <- c(covnames_new, "factor(protection_status2)", "Matern(1 | longitude + latitude)")
  covnames_new_bis <- covnames_new[-which(covnames_new %in% c("factor(protection_status2)"))]

  covnames_combined <- paste0(covnames_new, collapse = " + ")
  covnames_combined2 <- paste0(covnames_new_bis, collapse = " + ")
  
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
    noise_function(avoid = c("survey_id", "protection_status2"),
                   limit = 6)
  
  # add covariates
  sp <- dplyr::inner_join(sp, covariates_sp, by = "survey_id")
  
  sp$protection_status2 <- as.factor(as.character(sp$protection_status2))
  
  # Fit the model
  
  if(length(unique(sp$protection_status2)) == 1){
    
    model_fit <- spaMM::fitme(model_formula2, data = sp, method = "ML")
    
  }else{
    

    model_fit <- spaMM::fitme(model_formula, data = sp, method = "ML")
    
  }
  
  cov <- sp[!colnames(sp) %in% c("survey_id", "biomass", "longitude", "latitude", "site_code")]
  
  partial_plot <- lapply(1:ncol(cov), function(j) {
    
    partial_plot_j <- pdp::partial(model_fit,
                                   train = sp, 
                                   pred.var = colnames(cov)[j],
                                   prob = TRUE, 
                                   rug = TRUE,
                                   type = "regression") |> 
      dplyr::rename(biomass = yhat) |> 
      dplyr::mutate(species_name = species_name,
                    realm = unique(realm$realm))
    
  })
  
  dir.create(base_dir)
  
  save(partial_plot, file = paste0(base_dir, stringr::str_replace_all(species_name, " ", "_"), ".Rdata"))
  
}