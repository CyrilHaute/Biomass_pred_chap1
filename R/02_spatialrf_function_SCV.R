# function to fit spatial Random Forest

biomass = biomass_scv
covariates = rls_covariates
species_name = colnames(biomass_scv[[1]]$fitting)[!colnames(biomass_scv[[1]]$fitting) %in% c("survey_id", "latitude", "longitude")]
base_dir = base_dir

spatialrf_function <- function(biomass, 
                               covariates,
                               species_name,
                               base_dir){
  
  species_j <- list()
  
  for(i in 1:length(biomass)) {
    
    print(paste0("cv ", i))
    
    # create raw biomass object and select cross validation set i
    raw_biomass <- biomass[[i]]
    
    fmla <<- as.formula(paste0("Biomass ~ ", paste0(colnames(covariates)[!colnames(covariates) %in% "survey_id"], collapse = " + ")))

    species_j[[i]] <- pbmcapply::pbmclapply(1:length(species_name), function(j){
      
      # select the jth species from the fitting set
      fitting <- raw_biomass$fitting[,c("survey_id", species_name[j])]
      
      # add covariates
      fitting <- dplyr::inner_join(fitting, covariates, by = "survey_id")
      
      # log10(x+1) transform biomass
      fitting[,species_name[j]] <- log10(fitting[,species_name[j]] + 1)
      
      # select the jth species from the validation set
      validation <- raw_biomass$validation[,c("survey_id", species_name[j])]
      
      # add covariates
      validation <- dplyr::inner_join(validation, covariates, by = "survey_id")
      
      # log10(x+1) transform biomass
      validation[,species_name[j]] <- log10(validation[,species_name[j]]  + 1) 
      
      # keep only absences from species life area 
      # load rls surveys info, we need ecoregion 
      load("new_data/new_raw_data/00_rls_surveys.Rdata")
      rls_surveys$survey_id <- as.character(rls_surveys$survey_id)

      fitting <- dplyr::inner_join(fitting, rls_surveys)
      validation <- dplyr::inner_join(validation, rls_surveys)
      
      zone_geo_fit <- fitting[which(fitting[,species_name[j]] > 0),]
      zone_geo_fit <- unique(zone_geo_fit$ecoregion)
      
      zone_geo_val <- validation[which(validation[,species_name[j]] > 0),]
      zone_geo_val <- unique(zone_geo_val$ecoregion)
      
      if(sjmisc::is_empty(zone_geo_val)){
        
        zone_geo_val <- zone_geo_fit
        
      }
      
      fitting <- fitting |>  
        dplyr::filter(ecoregion %in% zone_geo_fit) |> 
        dplyr::select("survey_id", "longitude", "latitude", species_name[j], colnames(rls_covariates)[!colnames(rls_covariates) %in% c("survey_id")]) |> 
        dplyr::rename(X = longitude,
                      Y = latitude)
      
      validation <- validation |>  
        dplyr::filter(ecoregion %in% zone_geo_val) |> 
        dplyr::select("survey_id", "longitude", "latitude", species_name[j], colnames(rls_covariates)[!colnames(rls_covariates) %in% c("survey_id")]) |> 
        dplyr::rename(X = longitude,
                      Y = latitude)
      
      # get biomass data
      biomass_only <- fitting[which(fitting[,species_name[j]] > 0),]
      biomass_only_val <- validation[which(validation[,species_name[j]] > 0),]
      
      if(nrow(biomass_only_val) == 0){
        
        biomass_only_val <- validation
        
      }
      
      # keep only two times more absences than observation  
      # get absence
      
      n_subsample_fit <- nrow(fitting[which(fitting[, species_name[j]] > 0),]) * 2
      n_subsample_val <- nrow(validation[which(validation[, species_name[j]] > 0),]) * 2
      
      absence_fit <- fitting[which(fitting[, species_name[j]] == 0),]
      absence_val <- validation[which(validation[, species_name[j]] == 0),]
      
      if(nrow(absence_fit) > 0) {
        
        replacement_fit <- ifelse(length(which(fitting[, species_name[j]] == 0)) < n_subsample_fit, T, F)
        
        absence_fit <- absence_fit[sample(which(absence_fit[, species_name[j]] == 0), n_subsample_fit, replace = replacement_fit),]
        
      }
      
      if(nrow(absence_val) > 0) {
        
        replacement_val <- ifelse(length(which(validation[, species_name[j]] == 0)) < n_subsample_val, T, F)
        
        absence_val <- absence_val[sample(which(absence_val[, species_name[j]] == 0), n_subsample_val, replace = replacement_val),]
        
      }
      
      # combine absence and presence
      biomass_final <<- rbind(biomass_only, absence_fit) |> 
        as.data.frame()
      biomass_validation <- rbind(biomass_only_val, absence_val) |> 
        as.data.frame()
      
      names(biomass_final)[names(biomass_final) == species_name[j]] <<- "Biomass"
      names(biomass_validation)[names(biomass_validation) == species_name[j]] <- "Biomass"
      
      # As some covariates are at the country level, it means you can have very few or even only one value for these covariates
      # Check for the number of values in each covariates and add noise if < 6 values
      
      n_values <- lapply(1:ncol(biomass_final[,!colnames(biomass_final) %in% c("survey_id", "Biomass", "effectiveness")]), function(i) {unique(biomass_final[,!colnames(biomass_final) %in% c("survey_id", "Biomass", "effectiveness")][,i])})
      n_values_val <- lapply(1:ncol(biomass_validation[,!colnames(biomass_validation) %in% c("survey_id", "Biomass", "effectiveness")]), function(i) {unique(biomass_validation[,!colnames(biomass_validation) %in% c("survey_id", "Biomass", "effectiveness")][,i])})
      
      names(n_values) <- colnames(biomass_final[,!colnames(biomass_final) %in% c("survey_id", "Biomass", "effectiveness")])
      names(n_values_val) <- colnames(biomass_validation[,!colnames(biomass_validation) %in% c("survey_id", "Biomass", "effectiveness")])
      
      little_cov <- names(n_values[which(sapply(1:length(n_values), function(i) {length(n_values[[i]])}) <= 6)])
      little_cov_val <- names(n_values_val[which(sapply(1:length(n_values_val), function(i) {length(n_values_val[[i]])}) <= 6)])
      
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
      
      if(sjmisc::is_empty(little_cov_val) == TRUE){
        
        biomass_validation <- biomass_validation
        
      }else{
        
        n_cov_val <- which(names(biomass_validation) %in% little_cov_val)
        
        noise_val <- lapply(1:length(n_cov_val), function(i) {
          
          abs(rnorm(nrow(biomass_validation), 0.01, 0.01))
          
        })
        
        if(length(noise_val) == 1){
          
          biomass_validation[,n_cov_val] <- biomass_validation[,n_cov_val] + unlist(noise_val)
          
        }else{
          
          biomass_validation[,n_cov_val] <- biomass_validation[,n_cov_val] + noise_val
          
        }
        
      }
      
      coords <<- biomass_final |> 
        dplyr::select(X, Y) |> 
        as.data.frame()

      ### FITTING MODELS 
      # fit the spatial random forests

      model_fit <- SpatialML::grf(formula = fmla,
                                  dframe = biomass_final,
                                  bw = 15,
                                  kernel = "adaptive",
                                  coords = coords,
                                  ntree = 1000,
                                  geo.weighted = FALSE)
      
      if(!any(is.na(model_fit) == TRUE)){
        
        validation_predict  <- tryCatch(SpatialML::predict.grf(object = model_fit, new.data = biomass_validation, x.var.name = "X" , y.var.name = "Y", local.w = 1, global.w = 0), error = function(e) NA)
        
        # back transform predictions
        validation_predict <- 10^(validation_predict) - 1
        
        validation_predict <- data.frame(survey_id = biomass_validation$survey_id,
                                         validation_predict = validation_predict)
        
        validation_observed <- biomass_validation[,c("survey_id", "Biomass")]
        
        validation_observed <- validation_observed |>
          dplyr::rename(validation_observed = Biomass)
        
        validation_observed$validation_observed <- 10^(validation_observed$validation_observed) - 1
        
        validation_obs_prd <- validation_predict |>
          dplyr::inner_join(validation_observed)
        
        validation_obs_prd
        
      }else{
        
        validation_obs_prd  <- NA
        
      }
      
      rm(model_fit)
      gc()

    }, mc.cores = parallel::detectCores() - 1)
    
  }

  validation_prediction <- parallel::mclapply(1:length(species_j[[1]]), function(i){
    
    species_i <- lapply(species_j, `[[`, i)
    
    species_i_bind <- do.call(rbind, species_i)
    
  }, mc.cores = 10)
  
  
  extracted_predictions <- dplyr::tibble(species_name = species_name, 
                                         fitted_model = "SPRF", 
                                         validation_observed = lapply(validation_prediction, '[[', "validation_observed"),
                                         validation_predict = lapply(validation_prediction, '[[', "validation_predict"))
  
  # save prediciton output in same file structure
  
  model_dir <- "sprf"
  
  dir.create(base_dir, recursive = T)
  
  save(extracted_predictions, file = paste0(base_dir, model_dir, "_extracted_predictions.RData"))
  
  rm(list=ls())
  gc()

}