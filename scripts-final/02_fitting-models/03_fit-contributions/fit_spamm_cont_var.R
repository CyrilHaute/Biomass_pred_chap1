# run glmm (spamm)

source('scripts-final/00_functions/contribution-functions/spamm_cont_var.R')

print('glmm (spamm) biomass covariates contribution')
spamm_function_cont(biomass = rls_biomass_cont,
                    covariates = spatial_covariates_cont,
                    species_name = names(rls_biomass_cont)[-1],
                    base_dir_cont   = base_dir_cont)