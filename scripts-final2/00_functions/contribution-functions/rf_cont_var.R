# biomass = rls_biomass_cont
# covariates = covariates
# species_name = names(rls_biomass_cont)[-1]
# base_dir_cont   = base_dir_cont
# model_path      = 'model_abunocc'
# contribution_path = 'contributions_abunocc'

rf_function_cont <- function(biomass = biomass, 
                             covariates = covariates, 
                             species_name = species_name, 
                             base_dir_cont        = 'results/rls',
                             model_path      = 'model', 
                             contribution_path = 'contributions'){
  
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

  contribution <- pbmclapply(2:length(raw_biomass), function(i){
    
    biomass <- raw_biomass[,c(1,i)] # select the ith species
    biomass <- inner_join(biomass, covariates, by = "SurveyID") # add covariates
    biomass[,2] <- log10(biomass[,2]+1) # log10(x+1) transorm biomass
    
    # get biomass data
    biomass_only <- biomass[which(biomass[,2] > 0),]
      
    # keep only absences from species life area 
    rls_sitesInfos <- readRDS("../Biomass_prediction/data/Cyril_data/RLS_sitesInfos.rds")
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

    model_fit <- randomForest(x = biomass_final[covNames_new],
                              y = biomass_final[,2],
                              ntree = 1000,
                              importance=FALSE)
      
    explainer_rf <- DALEX::explain(model = model_fit, 
                                   data = biomass_final[covNames_new], 
                                   y = biomass_final[,2])
      
    vip.25_rf <- model_parts(explainer = explainer_rf,
                             loss_function = loss_root_mean_square,
                             B = 25,
                             type = "difference")
      
    vip.25_rf <- vip.25_rf %>% 
      group_by(variable) %>% 
      dplyr::summarise(Dropout_loss = mean(dropout_loss),
                       sd_dropout_loss = sd(dropout_loss))
      
    vip.25_rf <- vip.25_rf %>% 
      filter(!variable %in% c("SurveyID", "_baseline_", "_full_model_"))

    }, mc.cores = detectCores() - 1)
  
  extracted_contributions <- tibble(species_name = species_name, 
                                    fitted_model = 'RF', 
                                    # estimate contribution
                                    contributions = lapply(contribution, '[[', 2),
                                    # estimate median contribution
                                    sd_contributions = lapply(contribution, '[[', 3))
  
  # save prediciton output in same file structure
  
  path = (here::here("results", "rls", "predictions"))
  
  extracted_contributions <- setNames(split(extracted_contributions, seq(nrow(extracted_contributions))), extracted_contributions$species_name)

  model_dir <- 'rf'
  contribution_final_path <- paste0(base_dir_cont, '/', contribution_path)
  dir.create(contribution_final_path, recursive = T)
  names.list <- species_name
  names(extracted_contributions) <- names.list
  lapply(names(extracted_contributions), function(df)
    saveRDS(extracted_contributions[[df]], file = paste0(contribution_final_path, '/', model_dir, '_', df, '.rds')))
  
  rm(list=ls())
  gc()
  
}
