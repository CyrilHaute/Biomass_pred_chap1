
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
read_sp_eco <- do.call(rbind, read_sp_eco)

read_sp_eco$survey_id <- as.numeric(read_sp_eco$survey_id)

remove_na <- which(sapply(read_sp_eco$survey_id, is.na))

if(length(remove_na) == 0){
  
  read_sp_eco <- read_sp_eco
  
}else{
  
  read_sp_eco <- read_sp_eco[-remove_na,]
  
}

read_sp_eco$validation_predict <- as.numeric(read_sp_eco$validation_predict)

remove_infinite <- which(sapply(read_sp_eco$validation_predict, is.infinite))

if(length(remove_infinite) == 0){
  
  read_sp_eco <- read_sp_eco
  
}else{
  
  read_sp_eco <- read_sp_eco[-remove_infinite,]
  
}

load("outputs/best_models.Rdata")

read_sp_eco <- read_sp_eco |> 
  dplyr::inner_join(best_models) |> 
  dplyr::filter(model == best_model)
read_sp_eco$survey_id <- as.character(read_sp_eco$survey_id)
read_sp_eco$validation_observed <- as.numeric(read_sp_eco$validation_observed)

load("data/new_derived_data/rls_covariates.RData")

species_name <- unique(read_sp_eco$species_name)

read_sp_eco$validation_predict <- log10(read_sp_eco$validation_predict) + 1

library(ggplot2)

dirpath <- "figures/response_dhw/"
dir.create(dirpath)

pbmcapply::pbmclapply(1:length(species_name), function(i) {
  
  test <- read_sp_eco[read_sp_eco$species_name == species_name[i],]
  
  test <- test |>
    dplyr::inner_join(rls_covariates)
  
  test <- ggplot(test, aes(x = validation_predict, y = max_5year_degree_heating_week)) +
    geom_line() +
    geom_smooth(se = FALSE) +
    labs(x = "predicted_biomass", title = species_name[i])
  
  ggsave(filename = paste0(dirpath, stringr::str_replace_all(species_name[i], pattern = " ", replacement = "_"), ".png"))
  
}, mc.cores = parallel::detectCores() - 1)
  