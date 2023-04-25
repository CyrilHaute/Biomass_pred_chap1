# run random forest

source('scripts-final/00_functions/contribution-functions/rf_cont_var.R')

print('rf biomass covariates contribution')
rf_function_cont(biomass = rls_biomass_cont,
                 covariates = covariates_cont,
                 species_name = names(rls_biomass_cont)[-1],
                 base_dir_cont   = base_dir_cont)
