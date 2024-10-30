
################## extract evaluation model global ##################

glm <- list.files("outputs/glm_prediction", full.names = TRUE)
rf <- list.files("outputs/rf_prediction", full.names = TRUE)
brt <- list.files("outputs/brt_prediction", full.names = TRUE)
gam <- list.files("outputs/gam_prediction", full.names = TRUE)
spamm <- list.files("outputs/spamm_prediction", full.names = TRUE)
sprf <- list.files("outputs/sprf_prediction", full.names = TRUE)

file_model <- c(glm, rf, brt, gam, spamm, sprf)

read_sp_eco <- lapply(1:length(file_model), function(i) {
  
  load(file_model[i])
  assign(as.character(i), cv_j_bind)
  
})

performance <- pbmcapply::pbmclapply(1:length(read_sp_eco), function(i) {
  
  sp <- read_sp_eco[[i]]
  sp$survey_id <- as.numeric(sp$survey_id)
  
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
  
  return(metric_summary)
  
}, mc.cores = parallel::detectCores() - 1)

performance_bind <- do.call(rbind, performance)

perf_model <- performance_bind |> 
  dplyr::group_by(model) |> 
  dplyr::summarise(medianr2 = median(r2, na.rm = TRUE),
                   meanr2 = mean(r2, na.rm = TRUE),
                   maxr2 = max(r2, na.rm = TRUE),
                   minr2 = min(r2, na.rm = TRUE),
                   sdr2 = sd(r2, na.rm = TRUE),
                   medianpearson = median(pearson, na.rm = TRUE),
                   meanpearson = mean(pearson, na.rm = TRUE),
                   maxpearson = max(pearson, na.rm = TRUE),
                   minpearson = min(pearson, na.rm = TRUE),
                   sdpearson = sd(pearson, na.rm = TRUE),
                   medianspearman = median(spearman, na.rm = TRUE),
                   meanspearman = mean(spearman, na.rm = TRUE),
                   maxspearman = max(spearman, na.rm = TRUE),
                   minspearman = min(spearman, na.rm = TRUE),
                   sdspearman = sd(spearman, na.rm = TRUE))


################## extract evaluation model per realm ##################


realm_full <- list.files("outputs/biomass_prediction_ecoregion", full.names = T)
realm_small <- list.files("outputs/biomass_prediction_ecoregion", full.names = F)

load_outputs <- lapply(1:length(realm_full), function(i) {
  
  realms_all <- list.files(realm_full[i], full.names = T)
  realms <- realms_all[which(grepl(pattern = "sprf", realms_all) == FALSE)]
  
  realm_j <- lapply(1:length(realms), function(j) {
    
    load(realms[j])
    assign(paste0("model_", j), sp_i)
    
  })
  
  sprf <- list.files(realms_all[which(grepl(pattern = "sprf", realms_all) == TRUE)], full.names = T)
  
  load_outputs_sprf <- list(lapply(1:length(sprf), function(i) {
    
    load(sprf[i])
    assign(paste0("model_", i), cv_j_bind)
    
  }))
  
  load_outputs <- c(realm_j, load_outputs_sprf)

})

names(load_outputs) <- realm_small

# perform validation (out the bag)

model_assessment <- lapply(1:length(load_outputs), function(i){
  
  realm_i <- load_outputs[[i]]
  
  models <- lapply(1:length(realm_i), function(j) {
    
    model_j <- realm_i[[j]]
    
    sp_k <- pbmcapply::pbmclapply(1:length(model_j), function(k) {
      
      sp <- model_j[[k]]
      
      if(!is.data.frame(sp)){
        
        cv_k_bind_mean <- data.frame(species_name = NA,
                                     model = NA,
                                     Intercept = NA, 
                                     Slope = NA,
                                     Pearson = NA,
                                     Spearman = NA)
        
      }else{
        
        cv_l <- lapply(1:length(unique(sp$cv)), function(l) {
          
          sp[sp$cv == unique(sp$cv)[l],] |>
            dplyr::do(metrics = biomass_assessment_metrics(predictions   = .$validation_predict, 
                                                           observations  = .$validation_observed)) |> 
            tidyr::unnest(metrics) |> 
            dplyr::mutate(species_name = unique(sp$species_name),
                          cv = l,
                          model = unique(sp$model))
          
        })
        
        cv_l_bind <- do.call(rbind, cv_l)
        
        cv_l_bind_mean <- cv_l_bind |> 
          dplyr::group_by(species_name, model) |> 
          dplyr::summarise(dplyr::across("Intercept":"Spearman", mean))
        
      }
      
    }, mc.cores = parallel::detectCores() - 1)
    
    sp_k_bind <- do.call(rbind, sp_k)
    sp_k_bind$realm <- names(load_outputs)[i]
    
    return(sp_k_bind)
    
  })
  
  models_bind <- do.call(rbind, models)
  
})

model_assessment_realm <- do.call(rbind, model_assessment)

model_assessment_realm <- model_assessment_realm |> 
  dplyr::group_split(species_name)

dir.create("outputs/model_assessment_validation", recursive = T)

save(model_assessment_realm, file = "outputs/model_assessment_validation/validation_realm.Rdata")

