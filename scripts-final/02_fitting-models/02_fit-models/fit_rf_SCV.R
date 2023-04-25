# run random forest

source('scripts-final/00_functions/model-functions/rf_function_SCV.R')

print('rf biomass prediction')
rf_function(biomass = rls_biomass_SCV,
            covariates = covariates,
            species_name = colnames(rls_biomass_SCV[[1]]$fitting[,-1]),
            base_dir = base_dir)
