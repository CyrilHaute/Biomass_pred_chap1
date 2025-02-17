# function to fit rf and assess covariates partial dependence

partial_function_rf <- function(biomass, 
                                covariates, 
                                species_name,
                                base_dir){
  
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
  
  # Fit the model
  
  model_fit <- ranger::ranger(x = sp[!colnames(sp) %in% c("survey_id", "biomass", "longitude", "latitude", "site_code")],
                              y = unlist(sp[colnames(sp) %in% "biomass"]),
                              num.trees = ncol(sp[!colnames(sp) %in% c("survey_id", "biomass", "longitude", "latitude", "site_code")]) * 10)
  
  cov <- sp[!colnames(sp) %in% c("survey_id", "biomass", "longitude", "latitude", "site_code")]
  
  partial_plot <- lapply(1:ncol(cov), function(j) {
    
    partial_plot_j <- pdp::partial(model_fit,
                                   train = sp, 
                                   pred.var = colnames(cov)[j]) |> 
      dplyr::rename(biomass = yhat) |> 
      dplyr::mutate(species_name = species_name,
                    realm = unique(realm$realm))
    
  })

  dir.create(base_dir)
  
  save(partial_plot, file = paste0(base_dir, stringr::str_replace_all(species_name, " ", "_"), ".Rdata"))
  
}