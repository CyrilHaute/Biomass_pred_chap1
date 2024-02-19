# function to fit spatial Random Forest and assess covariates relative importance

# biomass = biomass_contribution
# covariates = rls_covariates
# species_name = colnames(biomass_contribution)[!colnames(biomass_contribution) %in% c("survey_id", "latitude", "longitude")]
# base_dir_cont = base_dir

spatialrf_function_cont <- function(biomass,
                                    covariates,
                                    species_name,
                                    base_dir_cont){

  # create raw biomass object and select cross validation set i
  raw_biomass <- biomass
  
  fmla <<- as.formula(paste0("Biomass ~ ", paste0(colnames(covariates)[!colnames(covariates) %in% "survey_id"], collapse = " + ")))

  contribution <- pbmcapply::pbmclapply(1:length(species_name), function(j){
    
    # select the jth species from the fitting set
    fitting <- raw_biomass[,c("survey_id", species_name[j])]
    
    # add covariates
    fitting <- dplyr::inner_join(fitting, covariates, by = "survey_id")
    
    # log10(x+1) transform biomass
    fitting[,species_name[j]] <- log10(fitting[,species_name[j]] + 1)
    
    # keep only absences from species life area 
    # load rls surveys info, we need ecoregion 
    load("data/new_raw_data/00_rls_surveys.Rdata")
    rls_surveys$survey_id <- as.character(rls_surveys$survey_id)
    
    fitting <- dplyr::inner_join(fitting, rls_surveys[!colnames(rls_surveys) %in% "depth"])

    zone_geo_fit <- fitting[which(fitting[,species_name[j]] > 0),]
    zone_geo_fit <- unique(zone_geo_fit$ecoregion)

    fitting <- fitting |>  
      dplyr::filter(ecoregion %in% zone_geo_fit) |> 
      dplyr::select("survey_id", "longitude", "latitude", species_name[j], colnames(rls_covariates)[!colnames(rls_covariates) %in% c("survey_id")]) |> 
      dplyr::rename(X = longitude,
                    Y = latitude)

    # get biomass data
    biomass_only <- fitting[which(fitting[,species_name[j]] > 0),]

    # keep only two times more absences than observation  
    # get absence
    
    n_subsample_fit <- nrow(fitting[which(fitting[, species_name[j]] > 0),]) * 2

    absence_fit <- fitting[which(fitting[, species_name[j]] == 0),]

    if(nrow(absence_fit) > 0) {
      
      replacement_fit <- ifelse(length(which(fitting[, species_name[j]] == 0)) < n_subsample_fit, T, F)
      
      absence_fit <- absence_fit[sample(which(absence_fit[, species_name[j]] == 0), n_subsample_fit, replace = replacement_fit),]
      
    }

    # combine absence and presence
    biomass_final <<- rbind(biomass_only, absence_fit) |> 
      as.data.frame()

    names(biomass_final)[names(biomass_final) == species_name[j]] <<- "Biomass"

    # As some covariates are at the country level, it means you can have very few or even only one value for these covariates
    # Check for the number of values in each covariates and add noise if < 6 values
    
    n_values <- lapply(1:ncol(biomass_final[,!colnames(biomass_final) %in% c("survey_id", "Biomass", "effectiveness")]), function(i) {unique(biomass_final[,!colnames(biomass_final) %in% c("survey_id", "Biomass", "effectiveness")][,i])})

    names(n_values) <- colnames(biomass_final[,!colnames(biomass_final) %in% c("survey_id", "Biomass", "effectiveness")])

    little_cov <- names(n_values[which(sapply(1:length(n_values), function(i) {length(n_values[[i]])}) <= 6)])

    if(sjmisc::is_empty(little_cov) == TRUE){
      
      biomass_final <- biomass_final
      
    }else{
      
      n_cov <- which(names(biomass_final) %in% little_cov)
      
      noise <- lapply(1:length(n_cov), function(i) {
        
        abs(rnorm(nrow(biomass_final), 0.01, 0.01))
        
      })
      
      if(length(noise) == 1){
        
        biomass_final[,n_cov] <- biomass_final[,n_cov] + unlist(noise)
        
      }else{
        
        biomass_final[,n_cov] <- biomass_final[,n_cov] + noise
        
      }
      
    }
    
    coords <<- biomass_final |> 
      dplyr::select(X, Y) |> 
      as.data.frame()
    
    ### FITTING MODELS 
    # fit the spatial random forests
    
    model_fit <- tryCatch(SpatialML::grf(formula = fmla,
                                         dframe = biomass_final,
                                         bw = 15,
                                         kernel = "adaptive",
                                         coords = coords,
                                         ntree = 1000,
                                         geo.weighted = FALSE), error = function(e) NA)
    
    covnames_new_new <- colnames(covariates)[which(colnames(covariates) %in% "survey_id" == FALSE)]
   
    # Use the package DALEX to assess covariates relative importance
    # First create an explain object (a representation of your model, depend on the structure of the algorithm used)
    explainer_sprf <- DALEX::explain(model = model_fit[[1]],
                                     data = biomass_final[covnames_new_new],
                                     y = biomass_final[,"Biomass"],
                                     label = "ranger")
    
    # Compute a 25-permutation-based value of the RMSE for all explanatory variables  
    vip.25_sprf <- DALEX::model_parts(explainer = explainer_sprf,
                                      loss_function = DALEX::loss_root_mean_square, # Here we used the RMSE as our loss function
                                      B = 25, # Number of permutation
                                      type = "difference")
    
    # From the model_parts function you get 25 RMSE values for each covariates. 
    # Take the mean and assess the standard-deviation of the RMSE for each covariates to assess the error of the permutation method
    vip.25_sprf <- vip.25_sprf |> 
      dplyr::group_by(variable) |> 
      dplyr::summarise(Dropout_loss = mean(dropout_loss),
                       sd_dropout_loss = sd(dropout_loss))
    
    vip.25_sprf <- vip.25_sprf |> 
      dplyr::filter(!variable %in% c("_baseline_", "_full_model_"))

    }, mc.cores = 15)
  
  extracted_contributions <- dplyr::tibble(species_name = species_name, 
                                           fitted_model = "SPRF", 
                                           # estimate contribution
                                           contributions_and_sd = contribution)
  
  # save contribution output in same file structure
  
  model_dir <- "sprf"
  
  dir.create(base_dir_cont)
  
  save(extracted_contributions, file = paste0(base_dir_cont, model_dir, "_extracted_contributions.RData"))
  
  rm(list=ls())
  gc()
  
}
