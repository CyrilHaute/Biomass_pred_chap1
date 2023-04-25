# run gam

source('scripts-final/00_functions/model-functions/gam_function_SCV.R')

print('gam biomass prediction')
gam_function(biomass = rls_biomass_SCV,
             covariates = covariates,
             species_name = colnames(rls_biomass_SCV[[1]]$fitting[,-1]),
             base_dir = base_dir)
