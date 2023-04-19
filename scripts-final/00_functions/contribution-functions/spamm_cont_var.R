# biomass = rls_biomass_cont
# covariates = spatial_covariates_cont
# species_name = names(rls_biomass_cont)[-1]
# base_dir_cont   = base_dir_cont
# contribution_path = 'contributions_biomass'

spamm_function_cont<- function(biomass = biomass, 
                               covariates = covariates,
                               species_name = species_name,
                               base_dir_cont   = 'results/rls',
                               contribution_path = 'predictions'){
  
  require(spaMM)
  require(pbmcapply)
  require(DALEX)

  # create raw biomass object
  raw_biomass <- biomass
  
  # response variable name
  response <- 'Biomass~'
    
  # rename covariates
  covNames_org <- names(covariates)
  covNames_org <- covNames_org[-which(covNames_org %in% c('SurveyID', 'Y', 'X'))]
  covNames_new <- names(covariates)
  covNames_new <- covNames_new[-which(covNames_new %in% c('SurveyID', 'Y', 'X'))]
  covNames_new <- c(covNames_new, "Matern(1 | X + Y)")
  covNames_new_bis <- covNames_new[-which(covNames_new %in% c('Effectiveness'))]
  
  # create formula with new covariate names
  fmla <- as.formula(paste("Biomass ~ ", paste(covNames_new, collapse= "+")))
  fmla2 <- as.formula(paste("Biomass ~ ", paste(covNames_new_bis, collapse= "+")))

  contribution <- pbmclapply(2:length(raw_biomass), function(i){
    
    biomass <- raw_biomass[,c(1,i)] # select the ith species
    biomass <- inner_join(biomass, covariates, by = "SurveyID") # add covariates
    biomass[,2] <- log10(biomass[,2]+1) # log10(x+1) transorm biomass
    
    # get biomass data
    biomass_only <- biomass[which(biomass[,2] > 0),]
    
    # keep only absences from species life area 
    rls_sitesInfos <- readRDS("data/Cyril_data/RLS_sitesInfos.rds")
    biomass <- inner_join(biomass, rls_sitesInfos, by = "SurveyID")
    biomass <- biomass[,-c(26:33,35:37)]
    zone_geo <- biomass[which(biomass[,2] > 0),]
    zone_geo <- unique(zone_geo$Ecoregion)
    biomass <- biomass %>% filter(Ecoregion %in% zone_geo)
    biomass <- biomass %>% dplyr::select(-Ecoregion)

    # get absence

    n_subsample <- length(biomass[which(biomass[,2] > 0),2])*2
        
    absence <- biomass[which(biomass[,2] == 0),]
        
    replacement <- ifelse(length(which(biomass[,2] == 0)) < n_subsample, T, F)
        
    absence <- absence[sample(which(absence[,2] == 0), n_subsample, replace = replacement),]

    # combine absence and presence
    biomass_final <- rbind(biomass_only, absence)
    namesp <- colnames(biomass_final[2])
    names(biomass_final)[names(biomass_final) == namesp] <- "Biomass"
    
    # As some covariates are at the country level, it means you can have very few or even only one value for these covariates
    # Check for the number of values in each covariates and add noise if < 6 values

    n_values <- lapply(5:ncol(biomass_final[,-25]), function(i) {unique(biomass_final[,i])})
    
    names(n_values) <- names(covariates)[-c(1:3,24)]
    
    little_cov <- names(n_values[which(sapply(1:length(n_values), function(i) {length(n_values[[i]])}) <= 6)])
    
    if(is_empty(little_cov) == TRUE){biomass_final <- biomass_final
    
    }else{
      
      n_cov <- which(names(covariates) %in% little_cov)
      n_cov <- n_cov + 1
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
    
    if(length(unique(biomass_final$Effectiveness)) == 1){
      
      model_fit <- fitme(fmla2, data = biomass_final, method = "ML")
      
    }else{
      
      model_fit <- fitme(fmla, data = biomass_final, method = "ML")
      
    }
    
    explainer_spamm <- DALEX::explain(model = model_fit, 
                                      data = biomass_final[,3:25], 
                                      y = biomass_final[,2])
         
    vip.25_spamm <- model_parts(explainer = explainer_spamm, 
                                loss_function = loss_root_mean_square,
                                B = 25,
                                type = "difference")
          
    vip.25_spamm <- vip.25_spamm %>% 
      group_by(variable) %>% 
      dplyr::summarise(Dropout_loss = mean(dropout_loss),
                       sd_dropout_loss = sd(dropout_loss))
          
    vip.25_spamm <- vip.25_spamm %>% 
      filter(!variable %in% c("SurveyID", "_baseline_", "_full_model_", "X", "Y"))
    
    }, mc.cores = detectCores() - 1)

  extracted_contributions <- tibble(species_name = species_name, 
                                    fitted_model = 'SPAMM', 
                                    # estimate contribution
                                    contributions = lapply(contribution, '[[', 2),
                                    # estimate median contribution
                                    sd_contributions = lapply(contribution, '[[', 3))
  
  # create prediction object to save
  
  path = (here::here("results", "rls", "predictions"))
  
  extracted_contributions <- setNames(split(extracted_contributions, seq(nrow(extracted_contributions))), extracted_contributions$species_name)
  
  model_dir <- 'spaMM'
  contribution_final_path <- paste0(base_dir_cont, '/', contribution_path)
  dir.create(contribution_final_path, recursive = T)
  names.list <- species_name
  names(extracted_contributions) <- names.list
  lapply(names(extracted_contributions), function(df)
    saveRDS(extracted_contributions[[df]], file = paste0(contribution_final_path, '/', model_dir, '_', df, '.rds')))

  rm(list=ls())
  gc()
  
}