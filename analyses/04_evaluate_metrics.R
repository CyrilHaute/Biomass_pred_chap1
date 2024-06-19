# extract evaluation 

# all functions for evaluating outputs
source("R/04_models_evaluation_functions.R")

output_files <- list.files("outputs/biomass_prediction", full.names = T)

load_outputs <- lapply(1:length(output_files), function(i) {
  
  load(output_files[i])
  assign(paste0("model_", i), sp_i)
  
})

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
