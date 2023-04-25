# function to fit random forest

biomass = rls_biomass_cont
covariates = covariates
species_name = names(rls_biomass_cont)[-1]
base_dir_cont   = base_dir_cont

rf_function_cont <- function(biomass = biomass, 
                             covariates = covariates, 
                             species_name = species_name, 
                             base_dir_cont = 'results/rls'){
  
  require(randomForest)
  require(pbmcapply)
  require(DALEX)

  # create raw biomass object
  raw_biomass <- biomass
  
  # response variable name
  response <- 'Biomass~'
  
  # rename and get general names for covariates for generalism
  
  covNames_new <- names(covariates) # randomForests take matrix which can be subset with this object 
  covNames_new <- covNames_new[-which(covNames_new %in% c('SurveyID'))]
i=3
  contribution <- pbmclapply(2:length(raw_biomass), function(i){
    
    biomass <- raw_biomass[,c(1,i)] # select the ith species
    biomass <- inner_join(biomass, covariates, by = "SurveyID") # add covariates
    biomass[,2] <- log10(biomass[,2]+1) # log10(x+1) transorm biomass
    
    # get biomass data
    biomass_only <- biomass[which(biomass[,2] > 0),]
      
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
    
    model_fit <- randomForest(x = biomass_final[covNames_new],
                              y = biomass_final[,2],
                              ntree = 1000,
                              importance=FALSE)
  
    # Use the package DALEX to assess covariates relative importance
    # First create an explain object (a representation of your model, depend on the structure of the algorithm used)
    explainer_rf <- DALEX::explain(model = model_fit, 
                                   data = biomass_final[covNames_new],
                                   y = biomass_final[,2],
                                   label = "randomForest")
      
    # Compute a 25-permutation-based value of the RMSE for all explanatory variables
    vip.25_rf <- DALEX::model_parts(explainer = explainer_rf,
                                    loss_function = loss_root_mean_square, # Here we used the RMSE as our loss function
                                    B = 25, # Number of permutation
                                    type = "difference")
      
    # From the model_parts function you get 25 RMSE values for each covariates. 
    # Take the mean and assess the standard-deviation of the RMSE for each covariates to assess the error of the permutation method
    vip.25_rf <- vip.25_rf %>% 
      dplyr::group_by(variable) %>% 
      dplyr::summarise(Dropout_loss = mean(dropout_loss),
                       sd_dropout_loss = sd(dropout_loss))
      
    vip.25_rf <- vip.25_rf %>% 
      dplyr::filter(!variable %in% c("SurveyID", "_baseline_", "_full_model_"))

    }, mc.cores = detectCores() - 1)
  
  extracted_contributions <- tibble(species_name = species_name, 
                                    fitted_model = 'RF', 
                                    # estimate contribution
                                    contributions = lapply(contribution, '[[', 2),
                                    # standard-deviation contribution between permutations
                                    sd_contributions = lapply(contribution, '[[', 3))
  
  # save contribution output in same file structure
  
  extracted_contributions <- setNames(split(extracted_contributions, seq(nrow(extracted_contributions))), extracted_contributions$species_name)

  model_dir <- 'rf'
  dir.create(base_dir_cont, recursive = T)
  names.list <- species_name
  names(extracted_contributions) <- names.list
  lapply(names(extracted_contributions), function(df)
    saveRDS(extracted_contributions[[df]], file = paste0(base_dir_cont, '/', model_dir, '_', df, '.rds')))
  
  rm(list=ls())
  gc()
  
}
