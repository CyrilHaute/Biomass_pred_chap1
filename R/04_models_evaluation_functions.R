# function for evaluating metrics

# predictions <- model_i$validation_predict
# predictions <- predictions[[1]]
# observations <- model_i$validation_observed
# observations <- observations[[1]]
# scale <- NULL

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

  # Linear model between values
  lm_test <- tryCatch(lm(predictions ~ observations), error = function(e) NA)

  # if the lm test is NA then some metrics cannot be calculated
  if(!is.na(lm_test[[1]][1])){
    
    sum_lm <- summary(lm_test)
    coef_lm <- coef(lm_test)
    
    # summaries from linear model
    residual_standard_error <- lm_test$sigma
    Intercept <- coef_lm[1]
    Slope <- coef_lm[2]
    Pr2 <- sum_lm$r.squared
    
    }else{
      
      Intercept <- NA
      Slope <- NA
      Pr2 <- NA
      
      }

  cor.test_pearson  <- tryCatch(cor.test(observations, predictions, method = 'pearson'), error = function(e) NA)
  cor.test_spearman <- tryCatch(cor.test(observations, predictions, method = 'spearman'), error = function(e) NA)
  sd_predictions <- tryCatch(sd(predictions, na.rm = T), error = function(e) NA)
  
  # Discrimination
  if(is.na(cor.test_pearson[[1]])){Pearson <- NA}else{Pearson <- cor.test_pearson$estimate}
  if(is.na(cor.test_spearman[[1]])){Spearman <- NA}else{Spearman  <- cor.test_spearman$estimate}
  
  metric_summary <- data.frame(Intercept = Intercept, 
                               Slope = Slope, 
                               Pearson = Pearson, 
                               Spearman = Spearman, 
                               Evaluation_number = length(observations), 
                               Scale = if(is.null(scale)){0}else{scale}, 
                               Evaluation_message = 'none')
  
  return(metric_summary)
  
}

