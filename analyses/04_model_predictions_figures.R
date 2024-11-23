# source functions ----
  
source("R/04_model_prediction_functions.R")

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
read_sp_eco <- read_sp_eco[-which(lapply(read_sp_eco, is.list) == FALSE)]

read_sp_eco <- pbmcapply::pbmclapply(1:length(read_sp_eco), function(i) {
  
  sp <- read_sp_eco[[i]]
  sp$survey_id <- as.numeric(sp$survey_id)
  to_remove <- which(sapply(sp$survey_id, is.na))
  if(length(to_remove) == 0){
    
    sp <- sp
    
  }else{
    
    sp <- sp[-to_remove,]
    
  }

  return(sp)
  
}, mc.cores = parallel::detectCores() - 1)

read_sp_eco <- do.call(rbind, read_sp_eco)

load("outputs/best_models.Rdata")

read_sp_eco_best <- read_sp_eco |> 
  dplyr::inner_join(best_models) |> 
  dplyr::filter(best_model == model)

# create plots 
observed_predicted_best_plot(input_data = read_sp_eco_best, 
                             nbins = 25)

# create plots 
observed_predicted_plot(input_data = read_sp_eco, 
                        nbins = 25, 
                        levels = c("glm", "gam", "spamm", "rf", "gbm", "sprf"))
 