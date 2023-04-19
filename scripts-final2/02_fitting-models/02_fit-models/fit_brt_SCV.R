# run boosted regression trees

source('scripts-final/00_functions/model-functions/SCV/brt_function_SCV.R')

print('brt biomass prediction')
brt_function(biomass = rls_biomass_SCV,
             covariates = covariates,
             species_name = colnames(rls_biomass_SCV[[1]]$fitting[,-1]),
             n.cores=1,
             base_dir        = base_dir,
             prediction_path = 'predictions_biomass')