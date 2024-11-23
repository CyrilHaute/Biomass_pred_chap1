
################## extract evaluation model global ##################

glm <- list.files("outputs/glm_prediction3", full.names = TRUE)
rf <- list.files("outputs/rf_prediction3", full.names = TRUE)
brt <- list.files("outputs/brt_prediction3", full.names = TRUE)
gam <- list.files("outputs/gam_prediction3", full.names = TRUE)
spamm <- list.files("outputs/spamm_prediction3", full.names = TRUE)
sprf <- list.files("outputs/sprf_prediction3", full.names = TRUE)

file_model <- c(glm, rf, brt, gam, spamm, sprf)

read_sp_eco <- lapply(1:length(file_model), function(i) {
  
  load(file_model[i])
  assign(as.character(i), cv_j_bind)
  
})

performance <- pbmcapply::pbmclapply(1:length(read_sp_eco), function(i) {
  
  sp <- read_sp_eco[[i]]
  
  if(is.list(sp) == FALSE) {
    
    metric_summary <- data.frame(species_name = NA,
                                 r2 = NA,
                                 r2_sd = NA,
                                 sd_pearson = NA,
                                 pearson = NA,
                                 sd_spearman = NA,
                                 spearman = NA,
                                 model = NA)
    
  }else{
    
    
    
    sp$survey_id <- as.numeric(sp$survey_id)
    sp$validation_observed <- as.numeric(sp$validation_observed)
    
    remove_na <- which(sapply(sp$survey_id, is.na))
    
    if(length(remove_na) == 0){
      
      sp <- sp
      
    }else{
      
      sp <- sp[-remove_na,]
      
    }
    
    sp$validation_predict <- as.numeric(sp$validation_predict)
    
    remove_infinite <- which(sapply(sp$validation_predict, is.infinite))
    
    if(length(remove_infinite) == 0){
      
      sp <- sp
      
    }else{
      
      sp <- sp[-remove_infinite,]
      
    }
    
    # Linear model between values
    lm_test <- lapply(1:length(unique(sp$cv)), function(k) {
      
      tryCatch(lm(validation_predict ~ validation_observed, data = sp[sp$cv %in% unique(sp$cv)[k],]), error = function(e) NA)
      
    })
    
    
    Pr2 <- lapply(1:length(lm_test), function(k) {
      
      sum_lm <- summary(lm_test[[k]])
      Pr2 <- sum_lm$r.squared
      
    })
    
    cor.test_pearson  <- lapply(1:length(unique(sp$cv)), function(l) {
      
      cor.test_pearson  <- tryCatch(cor.test(as.numeric(sp[sp$cv %in% unique(sp$cv)[l],]$validation_predict), 
                                             as.numeric(sp[sp$cv %in% unique(sp$cv)[l],]$validation_observed), method = 'pearson'), error = function(e) NA)
      
      Pearson <- cor.test_pearson$estimate
      
    })
    
    cor.test_spearman  <- lapply(1:length(unique(sp$cv)), function(l) {
      
      cor.test_spearman <- tryCatch(cor.test(as.numeric(sp[sp$cv %in% unique(sp$cv)[l],]$validation_predict), 
                                             as.numeric(sp[sp$cv %in% unique(sp$cv)[l],]$validation_observed), method = 'spearman'), error = function(e) NA)
      
      Spearman <- cor.test_spearman$estimate
      
    })
    
    metric_summary <- lapply(1:length(unique(sp$cv)), function(j) {
      
      data.frame(species_name = unique(sp$species_name),
                 R2 = unlist(Pr2[[j]]),
                 pearson = unlist(cor.test_pearson[[j]]),
                 spearman = unlist(cor.test_spearman[[j]]),
                 cv = unique(sp$cv)[j])
      
    })
    metric_summary <- do.call(rbind, metric_summary)
    
    metric_summary <- metric_summary |> 
      dplyr::group_by(species_name) |> 
      dplyr::summarise(r2 = mean(R2, na.rm = TRUE),
                       r2_sd = round(sd(R2, na.rm = TRUE), digits = 2),
                       sd_pearson = round(sd(pearson, na.rm = TRUE), digits = 2),
                       pearson = mean(pearson, na.rm = TRUE),
                       sd_spearman = round(sd(spearman, na.rm = TRUE), digits = 2),
                       spearman = mean(spearman, na.rm = TRUE)) |> 
      dplyr::mutate(model = unique(sp$model))
    
  }
  
  return(metric_summary)
  
}, mc.cores = parallel::detectCores() - 1)

performance_bind <- do.call(rbind, performance)

save(performance_bind, file = "outputs/performance_model.Rdata")

perf_model_gouped <- performance_bind |> 
  tidyr::drop_na() |> 
  dplyr::group_by(model) |>
  dplyr::summarise(medianpearson = median(pearson, na.rm = TRUE),
                   q95_pearson = quantile(pearson, probs = 0.95, na.rm = TRUE),
                   q05_pearson = quantile(pearson, probs = 0.05, na.rm = TRUE),
                   maxpearson = max(pearson, na.rm = TRUE),
                   minpearson = min(pearson, na.rm = TRUE),
                   sdpearson = sd(pearson, na.rm = TRUE),
                   medianspearman = median(spearman, na.rm = TRUE),
                   q95_spearman = quantile(spearman, probs = 0.95, na.rm = TRUE),
                   q05_spearman = quantile(spearman, probs = 0.05, na.rm = TRUE),
                   maxspearman = max(spearman, na.rm = TRUE),
                   minspearman = min(spearman, na.rm = TRUE),
                   sdspearman = sd(spearman, na.rm = TRUE))
perf_model_all <- performance_bind |> 
  tidyr::drop_na() |>
  dplyr::summarise(medianpearson = median(pearson, na.rm = TRUE),
                   q95_pearson = quantile(pearson, probs = 0.95, na.rm = TRUE),
                   q05_pearson = quantile(pearson, probs = 0.05, na.rm = TRUE),
                   maxpearson = max(pearson, na.rm = TRUE),
                   minpearson = min(pearson, na.rm = TRUE),
                   sdpearson = sd(pearson, na.rm = TRUE),
                   medianspearman = median(spearman, na.rm = TRUE),
                   q95_spearman = quantile(spearman, probs = 0.95, na.rm = TRUE),
                   q05_spearman = quantile(spearman, probs = 0.05, na.rm = TRUE),
                   maxspearman = max(spearman, na.rm = TRUE),
                   minspearman = min(spearman, na.rm = TRUE),
                   sdspearman = sd(spearman, na.rm = TRUE))
