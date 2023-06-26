# function to fit spatial Random Forest

spatialrf_function <- function(biomass = biomass, 
                               covariates = covariates,
                               species_name = species_name,
                               base_dir   = 'results/rls'){
  
  require(SpatialML)
  require(pbmcapply)

  predictions <- pbmclapply(1:length(biomass), function(i){
    
    # create raw biomass object
    raw_biomass <- biomass[[i]]
    
    # rename covariates
    
    covNames_org <- names(covariates)
    covNames_org <- covNames_org[-which(covNames_org %in% c('SurveyID', "Y",  "X"))]
    for(i in 1:length(covNames_org)){names(covariates)[3+i] <- paste0('cov', i)}
    covNames_new <- names(covariates) # randomForests take matrix which can be subset with this object 
    covNames_new <- covNames_new[-which(covNames_new %in% c('SurveyID', 'Y', 'X'))]
    fmla <<- as.formula(paste("Biomass ~ ", paste(covNames_new, collapse= "+")))

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
      biomass <- biomass[,-c(26:33,35:37)]
      zone_geo <- biomass[which(biomass[,2] > 0),]
      zone_geo <- unique(zone_geo$Ecoregion)
      biomass <- biomass %>% filter(Ecoregion %in% zone_geo)
      biomass <- biomass %>% dplyr::select(-Ecoregion)

      # keep only two times more absences than observation 
      # get absence

      n_subsample <- length(biomass[which(biomass[,2] > 0),2])*2
      
      absence <- biomass[which(biomass[,2] == 0),]
      
      replacement <- ifelse(length(which(biomass[,2] == 0)) < n_subsample, T, F)
      
      absence <- absence[sample(which(absence[,2] == 0), n_subsample, replace = replacement),]
      
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

      verification_predict  <- predict.grf(model_fit, biomass_only, x.var.name = "X" , y.var.name = "Y", local.w=1, global.w=0)
      validation_predict  <- predict.grf(model_fit, biomass_only_val, x.var.name = "X" , y.var.name = "Y", local.w=1, global.w=0)
        
      # back transform predictions
      verification_predict <- 10^(verification_predict)-1
      validation_predict <- 10^(validation_predict)-1
      validation_predict <- data.frame(SurveyID = biomass_only_val$SurveyID,
                                       validation_predict = validation_predict)
 
      predictions <- list(verification_predict, validation_predict)
      names(predictions) <- c("verification_predict", "validation_predict")
      predictions
        
    }, mc.cores = detectCores() - 5)
  }, mc.cores = 1)

  validation_prediction <- mclapply(1:length(predictions[[1]]), function(i){ # for each species, make mean, median and sd of fitting prediction across cross validation
    
    species_i <- lapply(predictions, `[[`, i)
    
    validation_prediction <- lapply(species_i, `[[`, 2)
    SurveyID <- lapply(validation_prediction, `[[`, 1)
    SurveyID <- unlist(SurveyID)
    SurveyID <- as.data.frame(sort(unique(SurveyID))) %>% rename(SurveyID = "sort(unique(SurveyID))")
    
    CV <- lapply(1:length(validation_prediction), function(i) {full_join(validation_prediction[[i]], SurveyID, by = "SurveyID")})
    CV <- lapply(1:length(CV), function(i) {
      
      cv_i <- CV[[i]]
      colnames(cv_i)[2] <- paste0("validation_predict_cv",i)
      cv_i
      
    })
    
    CV <- CV[[1]] %>% 
      inner_join(CV[[2]], by = "SurveyID") %>%
      inner_join(CV[[3]], by = "SurveyID") %>% 
      inner_join(CV[[4]], by = "SurveyID") %>% 
      inner_join(CV[[5]], by = "SurveyID") %>% 
      inner_join(CV[[6]], by = "SurveyID") %>% 
      inner_join(CV[[7]], by = "SurveyID") %>% 
      inner_join(CV[[8]], by = "SurveyID") %>% 
      inner_join(CV[[9]], by = "SurveyID") %>% 
      inner_join(CV[[10]], by = "SurveyID")
    CV <- CV[,-1]
    means_prediction <- as.matrix(CV) %>% rowMeans(na.rm = TRUE)
    medians_prediction <- as.matrix(CV) %>% rowMedians(na.rm = TRUE)
    sd_prediction<- mean((as.matrix(CV) %>% rowSds(na.rm = TRUE)), na.rm = TRUE)
    
    final_object <- list(means_prediction, medians_prediction, sd_prediction)
    names(final_object) <- c("means_prediction", "medians_prediction", "sd_prediction")
    final_object
    
    
  }, mc.cores = 10)
  
  validation_observed <- mclapply(1:length(biomass), function(i){
    
    cv_i <- biomass[[i]]$validation
    
    validation_observed <- mclapply(2:ncol(cv_i), function(j){
      
      species_j <- cv_i[,c(1,j)]
      
      species_j <- species_j[which(species_j[,2] > 0),]
      
    }, mc.cores = 1)
  }, mc.cores = 1)
  
  
  validation_observed <- mclapply(1:length(validation_observed[[1]]), function(i){
    
    species_i <- lapply(validation_observed, `[[`, i)
    
    SurveyID <- lapply(species_i, `[[`, 1)
    SurveyID <- unlist(SurveyID)
    SurveyID <- as.data.frame(sort(unique(SurveyID))) %>% rename(SurveyID = "sort(unique(SurveyID))")
    
    CV <- lapply(1:length(species_i), function(i) {full_join(species_i[[i]], SurveyID, by = "SurveyID")})
    CV <- lapply(1:length(CV), function(i) {
      
      cv_i <- CV[[i]]
      colnames(cv_i)[2] <- paste0("validation_predict_cv",i)
      cv_i
      
    })
    
    CV <- CV[[1]] %>% 
      inner_join(CV[[2]], by = "SurveyID") %>%
      inner_join(CV[[3]], by = "SurveyID") %>% 
      inner_join(CV[[4]], by = "SurveyID") %>% 
      inner_join(CV[[5]], by = "SurveyID") %>% 
      inner_join(CV[[6]], by = "SurveyID") %>% 
      inner_join(CV[[7]], by = "SurveyID") %>% 
      inner_join(CV[[8]], by = "SurveyID") %>% 
      inner_join(CV[[9]], by = "SurveyID") %>% 
      inner_join(CV[[10]], by = "SurveyID")
    CV <- CV[,-1]
    means_observed <- as.matrix(CV) %>% rowMeans(na.rm = TRUE)
    medians_observed <- as.matrix(CV) %>% rowMedians(na.rm = TRUE)
    
    final_object <- list(means_observed, medians_observed)
    names(final_object) <- c("means_observed", "medians_observed")
    final_object
    
  }, mc.cores = 10)
  
  extracted_predictions <- tibble(species_name = species_name, 
                                  fitted_model = 'SPRF', 
                                  # estimate mean predictions
                                  validation_observed_mean = lapply(validation_observed, '[[', 1),
                                  validation_predict_mean = lapply(validation_prediction, '[[', 1),
                                  # estimate median predictions
                                  validation_observed_median = lapply(validation_observed, '[[', 2),#list(validation_observed),
                                  validation_predict_median = lapply(validation_prediction, '[[', 2),
                                  # the amount of variation caused by cross validation
                                  sd_validation = lapply(validation_prediction, '[[', 3),
                                  MPA = NA)

  # create prediction object to save

  extracted_predictions <- setNames(split(extracted_predictions, seq(nrow(extracted_predictions))), extracted_predictions$species_name)

  model_dir <- "sprf"
  
  dir.create(base_dir, recursive = T)
  names.list <- species_name
  names(extracted_predictions) <- names.list
  lapply(names(extracted_predictions), function(df)
    saveRDS(extracted_predictions[[df]], file = paste0(base_dir, '/', model_dir, '_', df, '.rds')))

  rm(list=ls())
  gc()

}