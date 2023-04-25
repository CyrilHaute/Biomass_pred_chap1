# run gam

source('scripts-final/00_functions/contribution-functions/gam_cont_var.R')

print('gam biomass covariates contribution')
gam_function_cont(biomass = rls_biomass_cont,
                  covariates = covariates_cont,
                  species_name = names(rls_biomass_cont)[-1],
                  base_dir_cont = base_dir_cont)
