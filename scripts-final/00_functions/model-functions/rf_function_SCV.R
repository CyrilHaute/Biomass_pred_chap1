# Function to fit a Random Forest

rf_function <- function(biomass = biomass, 
                        covariates = covariates, 
                        species_name = species_name,
                        base_dir   = 'results/rls'){
  
  require(randomForest)
  require(pbmcapply)

  predictions <- pbmclapply(1:length(biomass), function(i){
    
    # create raw biomass object
    raw_biomass <- biomass[[i]]

    # rename covariates
  
    covNames_org <- names(covariates)
    covNames_org <- covNames_org[-which(covNames_org %in% c('SurveyID'))]
    for(i in 1:length(covNames_org)){names(covariates)[1+i] <- paste0('cov', i)}
    covNames_new <- names(covariates) # randomForests take matrix which can be subset with this object 
    covNames_new <- covNames_new[-which(covNames_new %in% c('SurveyID'))]

    species_j <- mclapply(2:length(raw_biomass$fitting), function(j){
      
      biomass <- raw_biomass$fitting[,c(1,j)] # select the jth species from the fitting set
      biomass <- inner_join(biomass, covariates, by = "SurveyID") # add covariates
      biomass[,2] <- log10(biomass[,2]+1) # log10(x+1) transorm biomass
      validation <- raw_biomass$validation[,c(1,j)] # select the jth species from the validation set
      validation <- inner_join(validation, covariates, by = "SurveyID")
      validation[,2] <- log10(validation[,2]+1) 
      
      # get biomass data
      biomass_only <- biomass[which(biomass[,2] > 0),]
      biomass_only_val <- validation[which(validation[,2] > 0),]
      
      # keep only absences from species life area 
      rls_sitesInfos <- readRDS("data/Cyril_data/RLS_sitesInfos.rds")
      biomass <- inner_join(biomass, rls_sitesInfos, by = "SurveyID")
      biomass <- biomass[,-c(24:31,33:35)]
      zone_geo <- biomass[which(biomass[,2] > 0),]
      zone_geo <- unique(zone_geo$Ecoregion)
      biomass <- biomass %>% filter(Ecoregion %in% zone_geo)
      biomass <- biomass[,-24]
      
      # keep only two times more absences than observation 
      # get absence

      n_subsample <- length(biomass[which(biomass[,2] > 0),2])*2
      
      absence <- biomass[which(biomass[,2] == 0),]
      
      replacement <- ifelse(length(which(biomass[,2] == 0)) < n_subsample, T, F)
      
      absence <- absence[sample(which(absence[,2] == 0), n_subsample, replace = replacement),]
      
      # combine absence and presence
      biomass_final <- rbind(biomass_only, absence)
      namesp <- colnames(biomass_final[2])
      names(biomass_final)[names(biomass_final) == namesp] <- "Biomass"

      model_fit <- randomForest(x = biomass_final[covNames_new],
                                y = biomass_final$Biomass,
                                ntree = 1000,
                                importance=FALSE)

      verification_predict  <- as.numeric(predict(model_fit, biomass_only, type = 'response'))
      validation_predict  <- as.numeric(predict(model_fit, biomass_only_val, type = 'response'))
            
      # back transform predictions
      verification_predict <- 10^(verification_predict)-1
      validation_predict <- 10^(validation_predict)-1
            
      predictions <- list(verification_predict, validation_predict)
      names(predictions) <- c("verification_predict", "validation_predict")
      predictions
      
      }, mc.cores = detectCores() - 5)
    }, mc.cores = 1)

  validation_prediction <- mclapply(1:length(predictions[[1]]), function(i){ # for each species, make mean, median and sd of fitting prediction across cross validation
    
    species_i <- lapply(predictions, `[[`, i)
    
    validation_prediction <- lapply(species_i, `[[`, 2)
    
    max_prediction <- max(sapply(1:length(validation_prediction), function(j) {length(validation_prediction[[j]])}))
    
    validation_prediction <- lapply(1:length(validation_prediction), function(k) {
      cv_k <- as.vector(validation_prediction[[k]])
      length(cv_k) <- max_prediction
      cv_k})
    
    validation_prediction <- as.matrix(do.call(cbind, validation_prediction))
    validation_prediction[which(is.finite(validation_prediction) == FALSE)] <- NA
    means_prediction  <- validation_prediction %>% rowMeans(na.rm = TRUE)
    means_prediction <- as.vector(means_prediction)
    medians_prediction <- validation_prediction %>% rowMedians(na.rm = TRUE)
    medians_prediction <- as.vector(medians_prediction)
    sd_prediction <- validation_prediction %>% rowSds(na.rm = TRUE)
    sd_prediction <- mean(sd_prediction, na.rm = TRUE)
    
    final_object <- list(means_prediction, medians_prediction, sd_prediction)
    names(final_object) <- c("means_prediction", "medians_prediction", "sd_prediction")
    final_object
    
  }, mc.cores = 1)
  
  validation_observed <- mclapply(1:length(biomass), function(i){
    
    cv_i <- biomass[[i]]$validation
    
    validation_observed <- mclapply(2:ncol(cv_i), function(j){
      
      species_j <- cv_i[,c(1,j)]
      
      species_j <- species_j[which(species_j[,2] > 0),]
      
    }, mc.cores = 1)
  }, mc.cores = 1)
  
  validation_observed <- mclapply(1:length(validation_observed[[1]]), function(i){
    
    species_i <- lapply(validation_observed, `[[`, i)
    
    validation_observed <- lapply(species_i, `[[`, 2)
    
    max_observation <- max(sapply(1:length(validation_observed), function(j) {length(validation_observed[[j]])}))
    
    validation_observed <- lapply(1:length(validation_observed), function(k) {
      cv_k <- as.vector(validation_observed[[k]])
      length(cv_k) <- max_observation
      cv_k
    })
    
    validation_observed  <- as.matrix(do.call(cbind, validation_observed ))
    means_observed <- validation_observed  %>% rowMeans(na.rm = TRUE)
    medians_observed <- validation_observed  %>% rowMedians(na.rm = TRUE)
    
    final_object <- list(means_observed, medians_observed)
    names(final_object) <- c("means_observed", "medians_observed")
    final_object
    
  }, mc.cores = 1)

  extracted_predictions <- tibble(species_name = species_name, 
                                  fitted_model = 'RF', 
                                  # estimate mean predictions
                                  validation_observed_mean = lapply(validation_observed, '[[', 1),
                                  validation_predict_mean = lapply(validation_prediction, '[[', 1),
                                  # estimate median predictions
                                  validation_observed_median = lapply(validation_observed, '[[', 2),
                                  validation_predict_median = lapply(validation_prediction, '[[', 2),
                                  # the amount of variation caused by cross validation
                                  sd_validation = lapply(validation_prediction, '[[', 3),
                                  MPA = NA)

  # create prediction object to save
  
  extracted_predictions <- setNames(split(extracted_predictions, seq(nrow(extracted_predictions))), extracted_predictions$species_name)
  
  model_dir <- "rf"
  
  dir.create(base_dir, recursive = T)
  names.list <- species_name
  names(extracted_predictions) <- names.list
  lapply(names(extracted_predictions), function(df)
    saveRDS(extracted_predictions[[df]], file = paste0(base_dir, '/', model_dir, '_', df, '.rds')))
  
  rm(list=ls())
  gc()
  
  }
  