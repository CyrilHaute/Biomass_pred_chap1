# Function for fitting boosted regression tree abundance models

# biomass = rls_biomass_cont
# covariates = covariates
# species_name = names(rls_biomass_cont)[-1]
# base_dir_cont   = base_dir_cont
# contribution_path = 'contributions_biomass'
# n.cores=1

brt_function_cont <- function(biomass = biomass, 
                              covariates = covariates, 
                              species_name = species_name, 
                              n.cores=1,
                              base_dir_cont = base_dir_cont){
  
  require(gbm)
  require(DALEX)
  require(pbmcapply)

  tryCatch({
    
    # create raw biomass object
    raw_biomass <- biomass
    
    # response variable name
    response <- 'Biomass~'
    
    # rename and get general names for covariates for generalism
    
    covNames_org <- names(covariates)
    covNames_org <- covNames_org[-which(covNames_org %in% c('SurveyID'))]

    # create formula
    covs     <- paste(covNames_org, collapse = '+')
    brt_formula <- as.formula(paste(response, covs))

    contribution <- pbmclapply(2:length(raw_biomass), function(i){
      
      biomass <- raw_biomass[,c(1,i)] # select the ith species
      biomass <- inner_join(biomass, covariates, by = "SurveyID") # add covariates
      biomass[,2] <- log10(biomass[,2]+1) # log10(x+1) transorm biomass
      
      # keep only absences from species life area 
      rls_sitesInfos <- readRDS("data/Cyril_data/RLS_sitesInfos.rds")
      biomass <- inner_join(biomass, rls_sitesInfos, by = "SurveyID")
      biomass <- biomass[,-c(24:31,33:35)]
      zone_geo <- biomass[which(biomass[,2] > 0),]
      zone_geo <- unique(zone_geo$Ecoregion)
      biomass <- biomass %>% filter(Ecoregion %in% zone_geo)
      biomass <- biomass[,-24]
      
      # get biomass data
      biomass_only <- biomass[which(biomass[,2] > 0),]
      
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
      
      # Fit the model
        
      model_fit <- tryCatch(gbm(formula = brt_formula,
                                data = biomass_final, 
                                distribution = "gaussian", 
                                n.trees = 10000,
                                interaction.depth = 3, 
                                shrinkage = 0.001,
                                bag.fraction = 0.8, 
                                cv.folds = 10, 
                                n.cores = n.cores), error = function(e) NA)
        
      if(!is.na(model_fit[[1]])){
        
        # selecting the best number of trees from cross validations
        gbm.mod.perf <- gbm.perf(model_fit, method = "cv", plot.it = F) 
          
        # fit model to all data
        model_fit <- gbm(formula = brt_formula,
                         data = biomass_final, 
                         distribution = "gaussian", 
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
                                                    family = "gaussian", # cannot have multinomial in this package, but also cannot fit models with single covariate and cross validation in gbm
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
                                  distribution = "gaussian", 
                                  n.trees = gbm.mod.perf,
                                  bag.fraction = 0.8,
                                  interaction.depth = 3, 
                                  shrinkage = 0.001), error = identity)
        
        if(class(model_fit)[1] == 'simpleError'){
          verification_predict <- NA
          validation_predict <- NA
          next
          
        }
        
        }

      # Use the package DALEX to assess covariates relative importance
      # First create an explain object (a representation of your model, depend on the structure of the algorithm used)
      explainer_gbm <- DALEX::explain(model = model_fit, 
                                      data = biomass_final[covNames_org], 
                                      y = biomass_final[,2],
                                      label = "gbm")
        
      # Compute a 25-permutation-based value of the RMSE for all explanatory variables
      vip.25_gbm <- DALEX::model_parts(explainer = explainer_gbm, 
                                       loss_function = loss_root_mean_square, # Here we used the RMSE as our loss function
                                       B = 25, # Number of permutation
                                       type = "difference")
        
      # From the model_parts function you get 25 RMSE values for each covariates. 
      # Take the mean and assess the standard-deviation of the RMSE for each covariates to assess the error of the permutation method
      vip.25_gbm <- vip.25_gbm %>% 
        dplyr::group_by(variable) %>% 
        dplyr::summarise(Dropout_loss = mean(dropout_loss),
                         sd_dropout_loss = sd(dropout_loss))
        
      vip.25_gbm <- vip.25_gbm %>% 
        dplyr::filter(!variable %in% c("SurveyID", "_baseline_", "_full_model_"))

      }, mc.cores = detectCores() - 1)

    extracted_contributions <- tibble(species_name = species_name, 
                                      fitted_model = 'GBM', 
                                      # estimate contribution
                                      contributions = lapply(contribution, '[[', 2),
                                      # standard-deviation contribution between permutations
                                      sd_contributions = lapply(contribution, '[[', 3))
    
    # save contribution output in same file structure
    
    extracted_contributions <- setNames(split(extracted_contributions, seq(nrow(extracted_contributions))), extracted_contributions$species_name)
    
    model_dir <- "brt"
    dir.create(base_dir_cont, recursive = T)
    names.list <- species_name
    names(extracted_contributions) <- names.list
    lapply(names(extracted_contributions), function(df)
      saveRDS(extracted_contributions[[df]], file = paste0(base_dir_cont, '/', model_dir, '_', df, '.rds')))
    
  }, 
  error = function(e) NA)
  
  rm(list=ls())
  gc()
  
} # end of function
