# function to fit boosted regression tree

brt_function <- function(biomass = biomass, 
                         covariates = covariates, 
                         species_name = species_name,
                         n.cores=1,
                         base_dir = 'results/rls'){
  
  require(gbm)
  require(pbmcapply)

  tryCatch({

    predictions <- pbmclapply(1:length(rls_biomass_SCV), function(i){
      
      # create raw biomass object
      raw_biomass <- biomass[[i]]
      
      # response variable name
      response <- 'Biomass~'

      # rename covariates
      
      covNames_org <- names(covariates)
      covNames_org <- covNames_org[-which(covNames_org %in% c('SurveyID'))]
      for(i in 1:length(covNames_org)){names(covariates)[1+i] <- paste0('cov', i)}
      covNames_new <- names(covariates)
      covNames_new <- covNames_new[-which(covNames_new %in% c('SurveyID'))]

      # create formula
      covNames_combined <- paste(covNames_new, collapse = '+')
      brt_formula <- as.formula(paste(response, covNames_combined))

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

        ### FITTING BOOSTED REGRESSION TREES
        
        # boosted regression tree modelling in 'gbm'

        # fit models with gbm packages
        # fit boosted regression tree with stochastic gradient boosting (i.e., bag.fraction != 1)
        model_fit <- tryCatch(gbm(formula = brt_formula,
                                  data = biomass_final, 
                                  distribution = 'gaussian', 
                                  n.trees = 10000,
                                  interaction.depth = 3, 
                                  shrinkage = 0.001,
                                  bag.fraction = 0.8, 
                                  cv.folds = 10, 
                                  n.cores = n.cores), error = function(e) NA)
          
        if(!is.na(model_fit[1])){
            
        # selecting the best number of trees from cross validations
          gbm.mod.perf <- gbm.perf(model_fit, method = "cv", plot.it = F) 
            
          # fit model to all data
          model_fit <- gbm(formula = brt_formula,
                           data = biomass_final, 
                           distribution = 'gaussian', 
                           n.trees = gbm.mod.perf,
                           bag.fraction = 0.8,
                           interaction.depth = 3, 
                           shrinkage = 0.001)

        }else{ # if the optimal number of trees cannot be identified in gbm

          # find the best model using gbm
          times = 0
          gbm.mod.perf <- NULL
          while(is.null(gbm.mod.perf)){
            gbm.mod.perf <- tryCatch(dismo::gbm.step(data = data.frame(biomass_final),
                                                     gbm.x = 3:ncol(data.frame(biomass_final)),
                                                     gbm.y = 2,
                                                     family = 'gaussian', # cannot have multinomial in this package, but also cannot fit models with single covariate and cross validation in gbm
                                                     tree.complexity = 3,
                                                     learning.rate = 0.001-(0.0001*times),
                                                     n.folds = 10, 
                                                     bag.fraction = 0.8, 
                                                     plot.main = F)$n.trees, error = function(e) NULL)
            times <- times + 1
            if(times == 500){stop()}
          }

          model_fit <- tryCatch(gbm(formula = brt_formula,
                                    data = biomass_final, 
                                    distribution = 'gaussian', 
                                    n.trees = gbm.mod.perf,
                                    bag.fraction = 0.8,
                                    interaction.depth = 3, 
                                    shrinkage = 0.001), error = identity)
            
          if(class(model_fit)[1] == 'simpleError'){
            verification_predict  <- NA
            validation_predict <- NA
            next
          }
            
        } # end of if the optimal number of trees cannot be identified in gbm

        predictions_verification <- predict(model_fit, newdata = biomass_only, n.trees = gbm.mod.perf, type = 'response')
        verification_predict  <- as.numeric(as.character(predictions_verification))
        predictions_validation <- predict(model_fit, newdata = biomass_only_val, n.trees = gbm.mod.perf, type = 'response')
        validation_predict  <- as.numeric(as.character(predictions_validation))
        
        # back transform predictions
        verification_predict <- 10^(verification_predict)-1
        validation_predict <- 10^(validation_predict)-1
        validation_predict <- data.frame(SurveyID = biomass_only_val$SurveyID,
                                         validation_predict = validation_predict)
        
        predictions <- list(verification_predict, validation_predict)
        names(predictions) <- c("verification_predict", "validation_predict")
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
    
    extracted_predictions <- tibble(species_name = species_name, 
                                    fitted_model = 'GBM', 
                                    # estimate mean predictions
                                    validation_observed = lapply(validation_observed, '[[', 2),
                                    validation_predict = lapply(validation_prediction, '[[', 2),
                                    MPA = NA)

  ### SAVE

  # create prediction object to save

  extracted_predictions <- setNames(split(extracted_predictions, seq(nrow(extracted_predictions))), extracted_predictions$species_name)
  
  model_dir <- "brt"
  
  dir.create(base_dir, recursive = T)
  names.list <- species_name
  names(extracted_predictions) <- names.list
  lapply(names(extracted_predictions), function(df)
    saveRDS(extracted_predictions[[df]], file = paste0(base_dir, '/', model_dir, '_', df, '.rds')))
  
  }, 
  error = function(e) NA)
  
  rm(list=ls())
  gc()
  
} # end of function
