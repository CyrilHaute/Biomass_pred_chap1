# run boosted regression trees

source('scripts-final/00_functions/contribution-functions/brt_cont_var.R')

print('brt biomass covariates contribution')
brt_function_cont(biomass = rls_biomass_cont,
                  covariates = covariates_cont,
                  species_name = names(rls_biomass_cont)[-1],
                  base_dir_cont   = base_dir_cont,
                  contribution_path = 'contributions_biomass')
