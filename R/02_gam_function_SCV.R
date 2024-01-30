# function to fit gams 

gam_function <- function(biomass = biomass, 
                         covariates = covariates,
                         species_name = species_name,
                         base_dir = 'results/rls'){
  
  require(mgcv)
  require(pbmcapply)

  predictions <- pbmclapply(1:length(biomass), function(i){
    
    # create raw biomass object
    raw_biomass <- biomass[[i]]
    
    # response variable name
    response <- 'Biomass~'
  
    # rename covariates
    covNames_org <- names(covariates)
    covNames_org <- covNames_org[-which(covNames_org %in% c('SurveyID'))]
    for(i in 1:length(covNames_org)){names(covariates)[1+i] <- paste0('cov', i)}
    covNames_new <- names(covariates)
    covNames_new <- covNames_new[-which(covNames_new %in% c('SurveyID', 'cov21'))]
    covNames_new <- c(covNames_new, "factor(cov21)")
    covNames_new_bis <- covNames_new[-which(covNames_new %in% c('factor(cov21)'))]
    
    # create formula with new covariate names
    covNames_splines <- paste0(rep('s(', length(covNames_new)),  covNames_new, rep(', k = 3)', length(covNames_new)))
    covNames_splines <- covNames_splines[-21]
    covNames_splines <- c(covNames_splines, "factor(cov21)")
    covNames_splines2 <- paste0(rep('s(', length(covNames_new_bis)),  covNames_new_bis, rep(', k = 3)', length(covNames_new_bis)))
    covNames_combined <- paste0(covNames_splines, collapse = '+')
    covNames_combined2 <- paste0(covNames_splines2, collapse = '+')
  
    # find the full model formula
    model_formula <- as.formula(paste0(response, covNames_combined))
    model_formula2 <- as.formula(paste0(response, covNames_combined2))

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
      
      if(nrow(absence) > 0) {
        
      replacement <- ifelse(length(which(biomass[,2] == 0)) < n_subsample, T, F)
        
      absence <- absence[sample(which(absence[,2] == 0), n_subsample, replace = replacement),]
      
      }
        
      # combine absence and presence
      biomass_final <- rbind(biomass_only, absence)
      namesp <- colnames(biomass_final[2])
      names(biomass_final)[names(biomass_final) == namesp] <- "Biomass"
        
      # As some covariates are at the country level, it means you can have very few or even only one value for these covariates
      # Check for the number of values in each covariates and add noise if < 6 values
      
      n_values <- lapply(3:ncol(biomass_final[,-23]), function(i) {unique(biomass_final[,i])})
      
      names(n_values) <- names(covariates)[-c(1,22)]
      
      little_cov <- names(n_values[which(sapply(1:length(n_values), function(i) {length(n_values[[i]])}) <= 6)])
      
      if(is_empty(little_cov) == TRUE){biomass_final <- biomass_final
      
      }else{
        
        n_cov <- gsub("cov", "", little_cov)
        n_cov <- as.numeric(n_cov)
        n_cov <- n_cov + 2
        noise <- lapply(1:length(n_cov), function(i) {
          abs(rnorm(nrow(biomass_final), 0.01, 0.01))
        })
        
        if(length(noise) == 1){
          
          biomass_final[,n_cov] <- biomass_final[,n_cov] + unlist(noise)
          
        }else{
          
          biomass_final[,n_cov] <- biomass_final[,n_cov] + noise
          
        }
        
      }
      
      # Fit the model

      if(length(unique(biomass_final$cov21)) == 1){
        
        model_fit <- gam(model_formula2, data=biomass_final, family=gaussian, select=FALSE, method = 'ML')
          
      }else{
        
        test <- unique(biomass_only_val$cov21) %in% unique(biomass_final$cov21)
        
        if(any(test == FALSE)){biomass_only_val <- biomass_only_val %>% filter(cov21 %in% unique(biomass_final$cov21))}
        
        model_fit <- gam(model_formula, data=biomass_final, family=gaussian, select=FALSE, method = 'ML')

        }
      
      verification_predict <- predict(model_fit, newdata = biomass_only, type = 'response')
      validation_predict <- predict(model_fit, newdata = biomass_only_val, type = 'response')
  
      # back transform predictions
      verification_predict <- 10^(verification_predict)-1
      validation_predict <- 10^(validation_predict)-1
      validation_predict <- data.frame(SurveyID = biomass_only_val$SurveyID,
                                       validation_predict = validation_predict)
        
      MPA <- ifelse(length(unique(biomass_final$cov21)) == 1, "no", "yes")
  
      predictions <- list(verification_predict, validation_predict, MPA)
      names(predictions) <- c("verification_predict", "validation_predict", "MPA")
      predictions
      
      }, mc.cores = detectCores() - 1)
    }, mc.cores = 1)

  validation_prediction <- mclapply(1:length(predictions[[1]]), function(i){ # for each species, make mean, median and sd of fitting prediction across cross validation
    
    species_i <- lapply(predictions, `[[`, i)
    
    validation_prediction <- lapply(species_i, `[[`, 2)
    test <- do.call(rbind, validation_prediction)
    
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
  
  MPA_test <- mclapply(1:length(predictions), function(i){
    
    cv_i <- predictions[[i]]
    
    MPA_test = lapply(cv_i, `[[`, 3)
    
  },mc.cores = 1)
  
  extracted_predictions <- tibble(species_name = species_name, 
                                  fitted_model = 'GAM', 
                                  # estimate mean predictions
                                  validation_observed = lapply(validation_observed, '[[', 2),
                                  validation_predict = lapply(validation_prediction, '[[', 2),
                                  MPA = MPA_test[[1]])

  # save prediction output in the same file structure
  
  extracted_predictions <- setNames(split(extracted_predictions, seq(nrow(extracted_predictions))), extracted_predictions$species_name)
  
  model_dir <- "gam"
  
  dir.create(base_dir, recursive = T)
  names.list <- species_name
  names(extracted_predictions) <- names.list
  lapply(names(extracted_predictions), function(df)
         saveRDS(extracted_predictions[[df]], file = paste0(base_dir, '/', model_dir, '_', df, '.rds')))
  
  rm(list=ls())
  gc()
  
}
