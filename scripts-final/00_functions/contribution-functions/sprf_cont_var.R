# biomass = rls_biomass_cont
# covariates = spatial_covariates_cont
# species_name = names(rls_biomass_cont)[-1]
# base_dir_cont   = base_dir_cont
# model_path      = 'model_abunocc'
# contribution_path = 'contributions_biomass'

spatialrf_function_cont <- function(biomass = biomass,
                                    covariates = covariates,
                                    species_name = species_name,
                                    base_dir_cont   = 'results/rls',
                                    model_path      = 'model', 
                                    contribution_path = 'contributions'){

  require(SpatialML)
  require(DALEX)
  require(pbmcapply)
    
  # create raw biomass object
  raw_biomass <- biomass
  
  # response variable name
  response <- 'Biomass~'

  covNames_org <- names(covariates)
  covNames_org <- covNames_org[-which(covNames_org %in% c('SurveyID', "Y",  "X"))]
  covNames_new <- names(covariates)
  covNames_new <- covNames_new[-which(covNames_new %in% c('SurveyID', 'Y', 'X'))]
  fmla <<- as.formula(paste("Biomass ~ ", paste(covNames_new, collapse= "+")))
    
  contribution <- pbmclapply(2:length(raw_biomass), function(i){
    
    biomass <- raw_biomass[,c(1,i)] # select the ith species
    biomass <- inner_join(biomass, covariates, by = "SurveyID") # add covariates
    biomass[,2] <- log10(biomass[,2]+1) # log10(x+1) transorm biomass
    
    # get biomass data
    biomass_only <- biomass[which(biomass[,2] > 0),]
    
    # keep only absences from species life area
    RLS_sitesInfos <- readRDS("../Biomass_prediction/data/Cyril_data/RLS_sitesInfos.rds")
    biomass <- inner_join(biomass, RLS_sitesInfos, by = "SurveyID")
    biomass <- biomass[,-c(26:33,35:37)]
    zone_geo <- biomass[which(biomass[,2] > 0),]
    zone_geo <- unique(zone_geo$Ecoregion)
    biomass <- biomass %>% filter(Ecoregion %in% zone_geo)
    biomass <- biomass[,-26]

    # keep only two times more absences than observation 
    # get absence

    n_subsample <- length(biomass[which(biomass[,2] > 0),2])*2
        
    absence <- biomass[which(biomass[,2] == 0),]
        
    replacement <- ifelse(length(which(biomass[,2] == 0)) < n_subsample, T, F)
        
    boot_absence <- absence[sample(which(absence[,2] == 0), n_subsample, replace = replacement),]
        
    # combine absence and presence
    biomass_final <<- rbind(biomass_only, absence)
    namesp <<- colnames(biomass_final[2])
    names(biomass_final)[names(biomass_final) == namesp] <<- "Biomass"
    
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
                      
    explainer_sprf <- DALEX::explain(model = model_fit[[1]],
                                     data = biomass_final[covNames_new],
                                     y = biomass_final[,2])
    
    vip.25_sprf <- model_parts(explainer = explainer_sprf,
                               loss_function = loss_root_mean_square,
                               B = 25,
                               type = "difference")
    
    vip.25_sprf <- vip.25_sprf %>%
      group_by(variable) %>% 
      dplyr::summarise(Dropout_loss = mean(dropout_loss),
                       sd_dropout_loss = sd(dropout_loss))
    
    vip.25_sprf <- vip.25_sprf %>% 
      filter(!variable %in% c("SurveyID", "_baseline_", "_full_model_"))

    },mc.cores = detectCores() - 1)
  
  extracted_contributions <- tibble(species_name = species_name, 
                                    fitted_model = 'SPRF', 
                                    # estimate contribution
                                    contributions = lapply(contribution, '[[', 2),
                                    # estimate median contribution
                                    sd_contributions = lapply(contribution, '[[', 3))
  
  # create prediction object to save
      
  path = (here::here("results", "rls", "predictions"))
      
  extracted_contributions <- setNames(split(extracted_contributions, seq(nrow(extracted_contributions))), extracted_contributions$species_name)
      
  model_dir <- 'Sprf'
  contribution_final_path <- paste0(base_dir_cont, '/', contribution_path)
  dir.create(contribution_final_path, recursive = T)
  names.list <- species_name
  names(extracted_contributions) <- names.list
  lapply(names(extracted_contributions), function(df)
    saveRDS(extracted_contributions[[df]], file = paste0(contribution_final_path, '/', model_dir, '_', df, '.rds')))
}
