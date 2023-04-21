# run spamm (GLMM)

source('scripts-final/00_functions/model-functions/spamm_function_SCV.R')

print('spamm biomass prediction')
spamm_function(biomass = rls_biomass_SCV,
               covariates = spatial_covariates,
               species_name = colnames(rls_biomass_SCV[[1]]$fitting[,-1]),
               base_dir        = base_dir, 
               prediction_path = 'predictions_biomass')
