# run spatial random forest

source('scripts-final/00_functions/model-functions/spatialrf_function_SCV.R')

print('sprf biomass prediction')
spatialrf_function(biomass = rls_biomass_SCV,
                   covariates = spatial_covariates,
                   species_name = colnames(rls_biomass_SCV[[1]]$fitting[,-1]),
                   base_dir        = base_dir,
                   prediction_path = 'predictions_biomass')
