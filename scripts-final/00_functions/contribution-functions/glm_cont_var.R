# function to fit glms 

# buildmer for model fitting and stepwise model selection of glmmTMB
# remotes::install_github("cvoeten/buildmer"); https://github.com/cvoeten/buildmer

# biomass = rls_biomass_cont
# covariates = covariates_cont
# species_name = names(rls_biomass_cont)[-1]
# base_dir_cont   = base_dir_cont

glm_function_cont <- function(biomass = biomass, 
                              covariates = covariates, 
                              species_name = species_name, 
                              base_dir_cont = base_dir_cont){
  
  require(DALEX)
  require(pbmcapply)
  
  # create raw biomass object
  raw_biomass <- biomass

  # response variable name
  response <- 'Biomass~'
    
  # rename covariates
  covNames_org <- names(covariates)
  covNames_org <- covNames_org[-which(covNames_org %in% c('SurveyID'))]
  covNames_org <- names(covariates)
  covNames_org <- covNames_org[-which(covNames_org %in% c('SurveyID', 'Effectiveness'))]
  covNames_new <- c(covNames_org, "factor(Effectiveness)")
  covNames_new_bis <- names(covariates)[-c(1,22)]
    
  # create formula with new covariate names
  covNames_new_2 <- paste0('I(', covNames_new, '^2',')')
  covNames_new_2 <- covNames_new_2[!covNames_new_2 == 'I(factor(Effectiveness)^2)']
  covNames_combined <- paste0(c(covNames_new, covNames_new_2), collapse = '+')
  covNames_combined2 <- paste0(c(covNames_new_bis, covNames_new_2), collapse = '+')
  model_formula <- as.formula(paste0(response, covNames_combined))
  model_formula2 <- as.formula(paste0(response, covNames_combined2))
  terms <- tabulate.formula(model_formula)
  terms2 <- tabulate.formula(model_formula2)
  form <- build.formula(dep='Biomass',terms)
  form2 <- build.formula(dep='Biomass',terms2)
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
        
    # As some covariates are at the country level, it means you can have very few or even only one value for these covariates
    # Check for the number of values in each covariates and add noise if < 6 values
      
    n_values <- lapply(3:ncol(biomass_final[,-23]), function(i) {unique(biomass_final[,i])})
      
    names(n_values) <- names(covariates)[-c(1,22)]
      
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
        
      model_fit <- glm(formula = form2, family = gaussian, data = biomass_final)
        
    }else{
        
      model_fit <- glm(formula = form, family = gaussian, data = biomass_final)
        
    }
      
    covNames_new <- covNames_new[-21]
    covNames_new <- c(covNames_new, "Effectiveness")

    # Use the package DALEX to assess covariates relative importance
    # First create an explain object (a representation of your model, depend on the structure of the algorithm used)
    explainer_glm <- DALEX::explain(model = model_fit, 
                                    data = biomass_final[covNames_new], 
                                    y = biomass_final[,2],
                                    label = "glm")
    
    # Compute a 25-permutation-based value of the RMSE for all explanatory variables  
    vip.25_glm <- DALEX::model_parts(explainer = explainer_glm, 
                                     loss_function = loss_root_mean_square, # Here we used the RMSE as our loss function
                                     B = 25, # Number of permutation
                                     type = "difference")

    # From the model_parts function you get 25 RMSE values for each covariates. 
    # Take the mean and assess the standard-deviation of the RMSE for each covariates to assess the error of the permutation method
    vip.25_glm <- vip.25_glm %>% 
      dplyr::group_by(variable) %>% 
      dplyr::summarise(Dropout_loss = mean(dropout_loss),
                       sd_dropout_loss = sd(dropout_loss))
        
    vip.25_glm <- vip.25_glm %>% 
      dplyr::filter(!variable %in% c("SurveyID", "_baseline_", "_full_model_"))

    }, mc.cores = detectCores() -1)
  
  extracted_contributions <- tibble(species_name = species_name, 
                                    fitted_model = 'GLM', 
                                    # estimate contribution
                                    contributions = lapply(contribution, '[[', 2),
                                    # standard-deviation contribution between permutations
                                    sd_contributions = lapply(contribution, '[[', 3))

  # save contribution output in same file structure
    
  extracted_contributions <- setNames(split(extracted_contributions, seq(nrow(extracted_contributions))), extracted_contributions$species_name)
    
  model_dir <- 'glm'
  dir.create(base_dir_cont)
  names.list <- species_name
  names(extracted_contributions) <- names.list
  lapply(names(extracted_contributions), function(df)
    saveRDS(extracted_contributions[[df]], file = paste0(base_dir_cont, '/', model_dir, '_', df, '.rds')))

  rm(list=ls())
  gc()
  
}
