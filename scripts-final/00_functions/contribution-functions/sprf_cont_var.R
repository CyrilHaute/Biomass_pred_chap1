# function to fit spatial Random Forest and assess covariates relative importance

spatialrf_function_cont <- function(biomass = biomass,
                                    covariates = covariates,
                                    species_name = species_name,
                                    base_dir_cont = base_dir_cont){

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
    RLS_sitesInfos <- readRDS("data/Cyril_data/RLS_sitesInfos.rds")
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
        
    absence <- absence[sample(which(absence[,2] == 0), n_subsample, replace = replacement),]
        
    # combine absence and presence
    biomass_final <<- rbind(biomass_only, absence)
    namesp <<- colnames(biomass_final[2])
    names(biomass_final)[names(biomass_final) == namesp] <- "Biomass"
    
    coords <<- biomass_final[,c(3,4)]

    # Fit the model

    model_fit <- grf(formula = fmla,
                     dframe = biomass_final,
                     bw = 15, # number of nearest neighbours
                     kernel = "adaptive", # adaptive kernel allow to select a number of nearest neighbours. Other option : "fixed" => buffer
                     coords = coords, # set X and Y coordinates (warning, order matters !)
                     ntree = 1000,
                     weighted = FALSE)
   
    # Use the package DALEX to assess covariates relative importance
    # First create an explain object (a representation of your model, depend on the structure of the algorithm used)
    explainer_sprf <- DALEX::explain(model = model_fit[[1]],
                                     data = biomass_final[covNames_new],
                                     y = biomass_final[,2],
                                     label = "ranger")
    
    # Compute a 25-permutation-based value of the RMSE for all explanatory variables  
    vip.25_sprf <- DALEX::model_parts(explainer = explainer_sprf,
                                      loss_function = loss_root_mean_square, # Here we used the RMSE as our loss function
                                      B = 25, # Number of permutation
                                      type = "difference")
    
    # From the model_parts function you get 25 RMSE values for each covariates. 
    # Take the mean and assess the standard-deviation of the RMSE for each covariates to assess the error of the permutation method
    vip.25_sprf <- vip.25_sprf %>%
      dplyr::group_by(variable) %>% 
      dplyr::summarise(Dropout_loss = mean(dropout_loss),
                       sd_dropout_loss = sd(dropout_loss))
    
    vip.25_sprf <- vip.25_sprf %>% 
      dplyr::filter(!variable %in% c("SurveyID", "_baseline_", "_full_model_"))

    },mc.cores = detectCores() - 1)
  
  extracted_contributions <- tibble(species_name = species_name, 
                                    fitted_model = 'SPRF', 
                                    # estimate contribution
                                    contributions = lapply(contribution, '[[', 2),
                                    # standard-deviation contribution between permutations
                                    sd_contributions = lapply(contribution, '[[', 3))
  
  # save contribution output in same file structure
      
  extracted_contributions <- setNames(split(extracted_contributions, seq(nrow(extracted_contributions))), extracted_contributions$species_name)
      
  model_dir <- 'Sprf'
  dir.create(base_dir_cont, recursive = T)
  names.list <- species_name
  names(extracted_contributions) <- names.list
  lapply(names(extracted_contributions), function(df)
    saveRDS(extracted_contributions[[df]], file = paste0(base_dir_cont, '/', model_dir, '_', df, '.rds')))
}
