# function to fit glms 

biomass = biomass_scv
covariates = rls_covariates
species_name = colnames(biomass_scv[[1]]$fitting)[!colnames(biomass_scv[[1]]$fitting) %in% c("survey_id", "latitude", "longitude")]
base_dir = base_dir

glm_function <- function(biomass, 
                         covariates,
                         species_name,
                         base_dir){
  
  # require(pbmcapply)

  predictions <- pbmcapply::pbmclapply(1:length(biomass), function(i){ #####################" mettre une boucle for ici
    
    print(paste0("cv ", i))
    
    # create raw biomass object and select cross validation set i
    raw_biomass <- biomass[[i]]
    
    # response variable name
    response <- 'Biomass~'

    # rename covariates
    covnames_new <- names(covariates)
    covnames_new <- covnames_new[-which(covnames_new %in% c("survey_id", "effectiveness"))]
    covnames_new <- c(covnames_new, "factor(effectiveness)")
    covnames_new_bis <- covnames_new[-which(covnames_new %in% c("factor(effectiveness)"))]
      
    # create formula with new covariate names
    covnames_new_poly <- paste0('I(', covnames_new, '^2',')')
    covnames_new_poly <- covnames_new_poly[-which(covnames_new_poly %in% c('I(factor(effectiveness)^2)'))]
    
    covnames_combined <- paste0(c(covnames_new, covnames_new_poly), collapse = " + ")
    covnames_combined2 <- paste0(c(covnames_new_bis, covnames_new_poly), collapse = " + ")
    
    model_formula <- as.formula(paste0(response, covnames_combined))
    model_formula2 <- as.formula(paste0(response, covnames_combined2))
    
    terms <- buildmer::tabulate.formula(model_formula)
    terms2 <- buildmer::tabulate.formula(model_formula2)
    
    form <- buildmer::build.formula(dep = "Biomass", terms)
    form2 <- buildmer::build.formula(dep = "Biomass", terms2)

    species_j <- pbmcapply::pbmclapply(1:length(species_name), function(j){
      
      # select the jth species from the fitting set
      fitting <- raw_biomass$fitting[,c("survey_id", species_name[j])]
      
      # add covariates
      fitting <- dplyr::inner_join(fitting, covariates, by = "survey_id")
      
      # log10(x+1) transform biomass
      fitting[,species_name[j]] <- log10(fitting[,species_name[j]] + 1)
      
      # select the jth species from the validation set
      validation <- raw_biomass$validation[,c("survey_id", species_name[j])]
      
      # add covariates
      validation <- dplyr::inner_join(validation, covariates, by = "survey_id")
      
      # log10(x+1) transform biomass
      validation[,species_name[j]] <- log10(validation[,species_name[j]]  + 1) 
  
      # get biomass data
      biomass_only <- fitting[which(fitting[,species_name[j]] > 0),]
      biomass_only_val <- validation[which(validation[,species_name[j]] > 0),]
        
      # keep only absences from species life area 
      # load rls surveys info, we need ecoregion 
      load("new_data/new_raw_data/00_rls_surveys.Rdata")
      rls_surveys$survey_id <- as.character(rls_surveys$survey_id)
      
      fitting <- dplyr::inner_join(fitting, rls_surveys)

      zone_geo <- fitting[which(fitting[,species_name[j]] > 0),]
      zone_geo <- unique(zone_geo$ecoregion)
      
      fitting <- fitting |>  
        dplyr::filter(ecoregion %in% zone_geo) |> 
        dplyr::select("survey_id", species_name[j], colnames(rls_covariates)[!colnames(rls_covariates) %in% c("survey_id")])

      # keep only two times more absences than observation  
      # get absence

      n_subsample <- nrow(fitting[which(fitting[, species_name[j]] > 0),]) * 2
          
      absence <- fitting[which(fitting[, species_name[j]] == 0),]
      
      if(nrow(absence) > 0) {
          
      replacement <- ifelse(length(which(fitting[, species_name[j]] == 0)) < n_subsample, T, F)
  
      absence <- absence[sample(which(absence[, species_name[j]] == 0), n_subsample, replace = replacement),]
      
      }
        
      # combine absence and presence
      biomass_final <- rbind(biomass_only, absence)

      names(biomass_final)[names(biomass_final) == species_name[j]] <- "Biomass"
  
      # As some covariates are at the country level, it means you can have very few or even only one value for these covariates
      # Check for the number of values in each covariates and add noise if < 6 values

      # n_values <- lapply(1:ncol(biomass_final[,colnames(covariates)[!colnames(covariates) %in% c("survey_id", "effectiveness")]]), function(i) {unique(biomass_final[,colnames(covariates)[!colnames(covariates) %in% c("survey_id", "effectiveness")]][,i])})

      n_values <- lapply(1:ncol(biomass_final[,!colnames(biomass_final) %in% c("survey_id", "Biomass", "effectiveness")]), function(i) {unique(biomass_final[,!colnames(biomass_final) %in% c("survey_id", "Biomass", "effectiveness")][,i])})
      
      names(n_values) <- colnames(biomass_final[,!colnames(biomass_final) %in% c("survey_id", "Biomass", "effectiveness")])
      
      little_cov <- names(n_values[which(sapply(1:length(n_values), function(i) {nrow(n_values[[i]])}) <= 6)])
      
      if(sjmisc::is_empty(little_cov) == TRUE){
        
        biomass_final <- biomass_final
      
      }else{
        
        n_cov <- which(names(biomass_final) %in% little_cov)
   
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
      
      if(length(unique(biomass_final$effectiveness)) == 1){
        
        model_fit <- tryCatch(glm(formula = form2, family = gaussian, data = biomass_final), error = function(e) NA)
        
        }else{
          
          test <- unique(biomass_only_val$effectiveness) %in% unique(biomass_final$effectiveness) # test if you have the same MPA effectivness factors in the fitting and validation set
          
          if(any(test == FALSE)){
            
            biomass_only_val <- biomass_only_val |>
              
              dplyr::filter(effectiveness %in% unique(biomass_final$effectiveness))
            
            }
          
          model_fit <- tryCatch(glm(formula = form, family = gaussian, data = biomass_final), error = function(e) NA)
          
        }
      
      if(!any(is.na(model_fit) == TRUE)){
        
        verification_predict  <- predict(model_fit, biomass_only, type = 'response')
        
        validation_predict <- predict(model_fit, biomass_only_val, type = 'response')
        
        # back transform predictions
        verification_predict <- 10^(verification_predict)-1
        
        validation_predict <- 10^(validation_predict)-1
        
        validation_predict <- data.frame(survey_id = biomass_only_val$survey_id,
                                         validation_predict = validation_predict)
          
        MPA <- ifelse(length(unique(biomass_final$effectiveness)) == 1, "no", "yes")
    
        predictions <- list(verification_predict, validation_predict, MPA)
        names(predictions) <- c("verification_predict", "validation_predict", "MPA")
        predictions
        
        }else{
          
          predictions <- NA
          
          }
        
        }, mc.cores = parallel::detectCores() - 1)
    
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
                                  fitted_model = 'GLM', 
                                  # estimate mean predictions
                                  validation_observed = lapply(validation_observed, '[[', 2),
                                  validation_predict = lapply(validation_prediction, '[[', 2),
                                  MPA = MPA_test[[1]])
  
  # save prediciton output in same file structure
  
  extracted_predictions <- setNames(split(extracted_predictions, seq(nrow(extracted_predictions))), extracted_predictions$species_name)
  
  model_dir <- "glm"
  
  dir.create(base_dir, recursive = T)
  names.list <- species_name
  names(extracted_predictions) <- names.list
  lapply(names(extracted_predictions), function(df)
    saveRDS(extracted_predictions[[df]], file = paste0(base_dir, '/', model_dir, '_', df, '.rds'))) 

  rm(list=ls())
  gc()
  
  }
