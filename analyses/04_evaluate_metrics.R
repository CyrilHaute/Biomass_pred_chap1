
################## extract evaluation model global ##################

# all functions for evaluating outputs
source("R/04_models_evaluation_functions.R")

sprf <- list.files("outputs/biomass_prediction/sprf", full.names = T)

load_outputs_sprf <- list(lapply(1:length(sprf), function(i) {
  
  load(sprf[i])
  assign(paste0("model_", i), cv_j_bind)
  
}))

output_files <- list.files("outputs/biomass_prediction", full.names = T)
output_files <- output_files[which(grepl(pattern = "sprf", output_files) == FALSE)]

load_outputs <- lapply(1:length(output_files), function(i) {
  
  load(output_files[i])
  assign(paste0("model_", i), sp_i)
  
})
load_outputs <- c(load_outputs, load_outputs_sprf)

# perform validation (out the bag)

model_assessment <- pbmcapply::pbmclapply(1:length(load_outputs), function(i){
  
  model_i <- load_outputs[[i]]
  
  sp_j <- pbmcapply::pbmclapply(1:length(model_i), function(j) {
    
    sp <- model_i[[j]]
    
    if(is.null(sp)){
      
      cv_k_bind_mean <- data.frame(species_name = NA,
                                   model = NA,
                                   Intercept = NA, 
                                   Slope = NA,
                                   Pearson = NA,
                                   Spearman = NA)
      
    }else{
        
      cv_k <- lapply(1:length(unique(sp$cv)), function(k) {
        
        sp[sp$cv == unique(sp$cv)[k],] |>
          dplyr::do(metrics = biomass_assessment_metrics(predictions   = .$validation_predict, 
                                                         observations  = .$validation_observed)) |> 
          tidyr::unnest(metrics) |> 
          dplyr::mutate(species_name = unique(sp$species_name),
                        cv = k,
                        model = unique(sp$model))
        
      })
      
      cv_k_bind <- do.call(rbind, cv_k)
      
      cv_k_bind_mean <- cv_k_bind |> 
        dplyr::group_by(species_name, model) |> 
        dplyr::summarise(dplyr::across("Intercept":"Spearman", mean))
      
      }

  }, mc.cores = parallel::detectCores() - 1)
  
  sp_j_bind <- do.call(rbind, sp_j)
  
}, mc.cores = 1)

model_assessment <- do.call(rbind, model_assessment)

model_assessment <- model_assessment |> 
  dplyr::group_split(species_name)

dir.create("outputs/model_assessment_validation", recursive = T)

save(model_assessment, file = "outputs/model_assessment_validation/validation.Rdata")


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

