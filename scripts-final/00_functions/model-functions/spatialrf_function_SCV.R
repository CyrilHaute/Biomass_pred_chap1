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
      
      load("new_data/new_raw_data/00_rls_surveys.Rdata")
      
      rls_surveys$survey_id <- as.character(rls_surveys$survey_id)
      
      # select the jth species from the fitting set
      fitting <- raw_biomass$fitting[,c("survey_id", species_name[j])]
      
      # add covariates
      fitting <- fitting |> 
        dplyr::inner_join(covariates, by = "survey_id") |> 
        dplyr::inner_join(rls_surveys[,c("survey_id", "latitude", "longitude")], by = "survey_id") |> 
        dplyr::select("survey_id", "latitude", "longitude", species_name[j], colnames(rls_covariates)[!colnames(rls_covariates) %in% c("survey_id")])
      
      # log10(x+1) transform biomass
      fitting[,species_name[j]] <- log10(fitting[,species_name[j]] + 1)
      
      # select the jth species from the validation set
      validation <- raw_biomass$validation[,c("survey_id", species_name[j])]
      
      # add covariates
      validation <- validation |> 
        dplyr::inner_join(covariates, by = "survey_id") |> 
        dplyr::inner_join(rls_surveys[,c("survey_id", "latitude", "longitude")], by = "survey_id") |> 
        dplyr::select("survey_id", "latitude", "longitude", species_name[j], colnames(rls_covariates)[!colnames(rls_covariates) %in% c("survey_id")])
      
      # log10(x+1) transform biomass
      validation[,species_name[j]] <- log10(validation[,species_name[j]]  + 1) 
      
      # get biomass data
      biomass_only <- fitting[which(fitting[,species_name[j]] > 0),]
      biomass_only_val <- validation[which(validation[,species_name[j]] > 0),]
      
      # keep only absences from species life area 
      # load rls surveys info, we need ecoregion 

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
        dplyr::select("survey_id", "latitude", "longitude", species_name[j], colnames(rls_covariates)[!colnames(rls_covariates) %in% c("survey_id")])
      
      validation <- validation |>  
        dplyr::filter(ecoregion %in% zone_geo_val) |> 
        dplyr::select("survey_id", "latitude", "longitude", species_name[j], colnames(rls_covariates)[!colnames(rls_covariates) %in% c("survey_id")])
      
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
      biomass_final <- rbind(biomass_only, absence_fit)
      biomass_validation <- rbind(biomass_only_val, absence_val)
      
      names(biomass_final)[names(biomass_final) == species_name[j]] <- "Biomass"

      coords <<- biomass_final[,c(3,4)]

      ### FITTING MODELS 
      # fit the spatial random forests

      model_fit <- grf(formula = fmla,
                       dframe = biomass_final,
                       bw = 15,
                       kernel = "adaptive",
                       coords = coords,
                       ntree = 1000,
                       weighted = FALSE)

      verification_predict  <- tryCatch(predict.grf(model_fit, biomass_only, x.var.name = "X" , y.var.name = "Y", local.w=1, global.w=0), error = function(e) NA)
      validation_predict  <- tryCatch(predict.grf(model_fit, biomass_only_val, x.var.name = "X" , y.var.name = "Y", local.w=1, global.w=0), error = function(e) NA)
        
      # back transform predictions
      if(!is.na(validation_predict)){
      verification_predict <- 10^(verification_predict)-1
      validation_predict <- 10^(validation_predict)-1
      validation_predict <- data.frame(SurveyID = biomass_only_val$SurveyID,
                                       validation_predict = validation_predict)}
 
      predictions <- list(verification_predict, validation_predict)
      names(predictions) <- c("verification_predict", "validation_predict")
      predictions
        
    }, mc.cores = detectCores() - 1)
  }, mc.cores = 1)

  validation_prediction <- mclapply(1:length(predictions[[1]]), function(i){ # for each species, make mean, median and sd of fitting prediction across cross validation
    
    species_i <- lapply(predictions, `[[`, i)
    
    validation_prediction <- lapply(species_i, `[[`, 2)
    test <- do.call(rbind, validation_prediction)
    test <- test[which(!is.na(test$validation_predict)),]
    
    }, mc.cores = 10)
  
  validation <- lapply(biomass, '[[', 2)
  
  validation_observed <- mclapply(2:ncol(validation[[1]]), function(i){
    
    names_col <- c("SurveyID", colnames(validation[[1]])[i])
    test <- lapply(validation, '[', names_col)
    test <- do.call(rbind, test)
    names(test) <- c("SurveyID", "Biomass")
    test <- test[test$Biomass > 0,]
    test <- test %>% filter(SurveyID %in% validation_prediction[[i-1]]$SurveyID)
    
  }, mc.cores = 1)
  
  extracted_predictions <- tibble(species_name = species_name, 
                                  fitted_model = 'SPRF', 
                                  # estimate mean predictions
                                  validation_observed = lapply(validation_observed, '[[', 2),
                                  validation_predict = lapply(validation_prediction, '[[', 2),
                                  MPA = NA)
  
  # create prediction object to save

  extracted_predictions <- setNames(split(extracted_predictions, seq(nrow(extracted_predictions))), extracted_predictions$species_name)

  model_dir <- "sprf"
  
  dir.create(base_dir, recursive = T)
  names.list <- species_name
  names(extracted_predictions) <- names.list
  lapply(names(extracted_predictions), function(df)
    saveRDS(extracted_predictions[[df]], file = paste0(base_dir, '/', model_dir, '_', df, '.rds')))

  rm(list=ls())
  gc()

}