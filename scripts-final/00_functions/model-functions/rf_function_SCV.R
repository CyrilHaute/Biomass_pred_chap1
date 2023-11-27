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

    species_j <- pbmclapply(2:length(raw_biomass$fitting), function(j){
      
      biomass <- raw_biomass$fitting[,c(1,j)] # select the jth species from the fitting set
      biomass <- inner_join(biomass, covariates, by = "SurveyID") # add covariates
      biomass[,2] <- log10(biomass[,2]+1) # log10(x+1) transform biomass
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
      # biomass <- biomass[,-c(15:22,24:26)]
      zone_geo <- biomass[which(biomass[,2] > 0),]
      zone_geo <- unique(zone_geo$Ecoregion)
      biomass <- biomass %>% filter(Ecoregion %in% zone_geo)
      biomass <- biomass[,-24]
      # biomass <- biomass[,-15]

      # keep only two times more absences than observation
      # get absence

      n_subsample <- length(biomass[which(biomass[,2] > 0),2])*2

      absence <- biomass[which(biomass[,2] == 0),]

      if(nrow(absence) > 0) {

        replacement <- ifelse(length(which(biomass[,2] == 0)) < n_subsample, T, F)

        absence <- absence[sample(which(absence[,2] == 0), n_subsample, replace = replacement),]

      }

      # combine absence and presence
      biomass_final <- rbind(biomass_only, absence)
      namesp <- colnames(biomass_final[2])
      names(biomass_final)[names(biomass_final) == namesp] <- "Biomass"

      model_fit <- tryCatch(randomForest(x = biomass_final[covNames_new],
                                         y = biomass_final$Biomass,
                                         ntree = 1000,
                                         importance=FALSE), error = function(e) NA)

      verification_predict  <- tryCatch(predict(model_fit, biomass_only, type = 'response'), error = function(e) NA)
      validation_predict <- tryCatch(predict(model_fit, biomass_only_val, type = 'response'), error = function(e) NA)

      # back transform predictions
      if(!any(is.na(validation_predict))){
      verification_predict <- 10^(verification_predict)-1
      validation_predict <- 10^(validation_predict)-1
      validation_predict <- data.frame(SurveyID = biomass_only_val$SurveyID,
                                       validation_predict = validation_predict)}

      predictions <- list(verification_predict, validation_predict)
      names(predictions) <- c("verification_predict", "validation_predict")
      predictions
      
      validation_predict <- inner_join(validation_predict, validation, by = "SurveyID")
      validation_predict$`Bodianus diplotaenia` <- 10^(validation_predict$`Bodianus diplotaenia`)-1
      ggplot(validation_predict, aes(y = validation_predict, x = `Bodianus diplotaenia`)) + geom_point() + geom_smooth(method = 'lm')
      ggplot(validation_predict, aes(y = validation_predict, x = cov21)) + geom_point() + geom_smooth(method = 'lm')
      
      
      }, mc.cores = detectCores() - 1)
    }, mc.cores = 1)

  validation_prediction <- mclapply(1:length(predictions[[1]]), function(i){ # for each species, make mean, median and sd of fitting prediction across cross validation
    
    species_i <- lapply(predictions, `[[`, i)
    
    validation_prediction <- lapply(species_i, `[[`, 2)
    test <- do.call(rbind, validation_prediction)
    if(ncol(test) == 2){
    test <- test[which(!is.na(test$validation_predict)),]}else{test = tibble(SurveyID = NA,
                                                                             validation_predict = NA)}
  }, mc.cores = 10)
  
  validation <- lapply(biomass, '[[', 2)
  
  validation_observed <- mclapply(2:ncol(validation[[1]]), function(i){
    
    names_col <- c("SurveyID", colnames(validation[[1]])[i])
    test <- lapply(validation, '[', names_col)
    test <- do.call(rbind, test)
    names(test) <- c("SurveyID", "Biomass")
    test <- test[test$Biomass > 0,]
    if(!is.na(validation_prediction[[i-1]])){
    test <- test %>% filter(SurveyID %in% validation_prediction[[i-1]]$SurveyID)}else{test <- tibble(SurveyID = NA,
                                                                                                     Biomass = NA)}
    
  }, mc.cores = 1)
  
    extracted_predictions <- tibble(species_name = species_name, 
                                  fitted_model = 'RF', 
                                  # estimate mean predictions
                                  validation_observed = lapply(validation_observed, '[[', 2),
                                  validation_predict = lapply(validation_prediction, '[[', 2),
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
  