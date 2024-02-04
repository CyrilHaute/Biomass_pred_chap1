# function for evaluating metrics

biomass_assessment_metrics <- function(predictions, observations, scale = NULL){
  
  # check lengths are the same
  if(length(observations) != length(predictions)){return(data.frame(Armse = NA, 
                                                                    Amae  = NA, 
                                                                    Intercept = NA, 
                                                                    Slope = NA, 
                                                                    Pearson = NA, 
                                                                    Spearman = NA, 
                                                                    Psd = NA, 
                                                                    Pdispersion = NA, 
                                                                    Pr2 = NA, 
                                                                    Evaluation_number = NA, 
                                                                    Scale = if(is.null(scale)){0}else{scale}, 
                                                                    Evaluation_message = 'warning: length of observations and predictions does not match locations'))}
  
  # keep only observations that are biomasses in the first input
  to_keep <- observations>0
  observations <- observations[to_keep]
  predictions  <- predictions[to_keep]

  # create dataframe of locations predictions and observations
  if(!is.null(scale)){
    
    # create dataframe
    ddata <- data.frame(observations, predictions)
    
    # aggregate by scale
    rescaled_observations <- ddata %>% 
      mutate(SiteLatitude = plyr::round_any(.$SurveyID, scale)) %>% 
      group_by(SurveyID) %>% 
      do(observations = mean(.$observations), 
         predictions = mean(.$predictions)) %>% 
      unnest(c('observations', 'predictions'))
    
    observations = rescaled_observations$observations
    predictions  = rescaled_observations$predictions
    
    }
  
  # Linear model between values
  lm_test <- tryCatch(lm(predictions ~ observations), error = function(e) NA)

  # if the lm test is NA then some metrics cannot be calculated
  if(!is.na(lm_test[[1]][1])){
    sum_lm            <- summary(lm_test)
    coef_lm           <- coef(lm_test)
    # summaries from linear model
    residual_standard_error <- lm_test$sigma
    Intercept <- coef_lm[1]
    Slope     <- coef_lm[2]
    Pr2         <- sum_lm$r.squared
  }else{
    Intercept <- NA
    Slope     <- NA
    Pr2        <- NA
  }
  
  tryCatch(cor.test(observations, predictions, method = 'pearson'), error = function(e) NA)
  cor.test_pearson  <-   tryCatch(cor.test(observations, predictions, method = 'pearson'), error = function(e) NA)
  cor.test_spearman <-   tryCatch(cor.test(observations, predictions, method = 'spearman'), error = function(e) NA)
  sd_predictions    <- tryCatch(sd(predictions, na.rm = T), error = function(e) NA)
  

  # Accuracy (A-overall)
  # root mean squared error between average predicted biomass
  # across all sites and observed biomass at each site
  Armse <- sqrt(mean((predictions - observations)^2, na.rm = T)) # root mean squared error gives more weight to big deviations
  Amae  <- mean(abs((predictions  - observations)), na.rm = T)    # mean absolute error weights all errors the same - if positive the value of observations is > the values of predictions
  
  # express as relative to the mean and median
  Armse_rel_mean <- Armse / mean(observations, na.rm = T)
  Armse_rel_median <- Armse / median(observations, na.rm = T)
  Amae_rel_mean <- Amae / mean(observations, na.rm = T)
  Amae_rel_median <- Amae / median(observations, na.rm = T)
  
  # Discrimination
  if(is.na(cor.test_pearson[[1]])){Pearson <- NA}else{Pearson <- cor.test_pearson$estimate}
  if(is.na(cor.test_spearman[[1]])){Spearman <- NA}else{Spearman  <- cor.test_spearman$estimate}
  
  # Precision
  if(is.na(sd_predictions)){Psd <- NA}else{Psd <- sd_predictions}
  if(is.na(sd_predictions)){Pdispersion <- NA}else{Pdispersion <- sd_predictions / sd(observations)}
  
  metric_summary <- data.frame(Armse = Armse, 
                               Armse_rel_mean = Armse_rel_mean, 
                               Armse_rel_median = Armse_rel_median,
                               Amae_rel_mean = Amae_rel_mean, 
                               Amae_rel_median = Amae_rel_median,
                               Amae  = Amae, 
                               Intercept = Intercept, 
                               Slope = Slope, 
                               Pearson = Pearson, 
                               Spearman = Spearman, 
                               Psd = Psd, 
                               Pdispersion = Pdispersion, 
                               Pr2 = Pr2, 
                               Evaluation_number = length(observations), 
                               Scale = if(is.null(scale)){0}else{scale}, 
                               Evaluation_message = 'none')
  
  return(metric_summary)
  
}

