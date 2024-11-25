# source functions ----

library(patchwork)

source("R/06_species_traits_figures_function.R")
source("R/05_contributions_figures_functions.R")

pal_sp_trait <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1"), PNWColors::pnw_palette("Bay", 6, type = "continuous"))

load("outputs/best_models.Rdata")

#### Covariates contribution plot ####

# Merge all files together by species

glm <- list.files("outputs/glm_contribution3", full.names = TRUE)
glm <- lapply(1:length(glm), function(i) {
  
  load(glm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
glm <- do.call(rbind, glm)

gam <- list.files("outputs/gam_contribution3", full.names = TRUE)
gam <- lapply(1:length(gam), function(i) {
  
  load(gam[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gam <- do.call(rbind, gam)

spamm <- list.files("outputs/spamm_contribution3", full.names = TRUE)
spamm <- lapply(1:length(spamm), function(i) {
  
  load(spamm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
spamm <- do.call(rbind, spamm)

rf <- list.files("outputs/rf_contribution3", full.names = TRUE)
rf <- lapply(1:length(rf), function(i) {
  
  load(rf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
rf <- do.call(rbind, rf)

sprf <- list.files("outputs/sprf_contribution3", full.names = TRUE)
sprf <- lapply(1:length(sprf), function(i) {
  
  load(sprf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
sprf <- do.call(rbind, sprf)
colnames(sprf)[colnames(sprf) %in% c("global_dropout_loss", "global_sd_dropout_loss")] <- c("Dropout_loss", "sd_dropout_loss")

gbm <- list.files("outputs/brt_contribution3", full.names = TRUE)
gbm <- lapply(1:length(gbm), function(i) {
  
  load(gbm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gbm <- do.call(rbind, gbm)

bind_files <- list(glm, gam, spamm, rf, sprf, gbm)
bind_files <- lapply(1:length(bind_files), function(i) {
  
  model_i <- bind_files[[i]]
  
  sp <- unique(model_i$species_name)
  
  new_sp <- pbmcapply::pbmclapply(1:length(sp), function(j) {
    
    sp_j <- model_i |> 
      dplyr::filter(species_name == sp[[j]])
    
    sum_dr <- sum(sp_j$Dropout_loss)
    
    sp_j <- sp_j |> 
      dplyr::rowwise() |> 
      dplyr::mutate(Dropout_loss = (Dropout_loss * 100) / sum_dr)
    
  }, mc.cores = parallel::detectCores() - 1)
  new_sp <- do.call(rbind, new_sp)
  
})
bind_files <- purrr::reduce(bind_files, dplyr::full_join)



fitted_model <- unique(bind_files$fitted_model)

# Load species traits
sp_car <- read.csv("data/new_raw_data/Traits_tropical_spp_1906.csv", header = TRUE) |> 
  dplyr::rename(species_name = Species) |> 
  dplyr::filter(species_name %in% unique(bind_files$species_name))
sp_car$ML_cat <- NA
sp_car <- sp_car[which(is.na(sp_car$MaxLength) == FALSE),]
sp_car <- sp_car[which(is.na(sp_car$Trophic_guild_name) == FALSE),]

sp_car[sp_car$MaxLength > 0 & sp_car$MaxLength <= 20,]$ML_cat <- "0-20 cm"
sp_car[sp_car$MaxLength > 20 & sp_car$MaxLength <= 40,]$ML_cat <- "20-40 cm"
sp_car[sp_car$MaxLength > 40 & sp_car$MaxLength <= 60,]$ML_cat <- "40-60 cm"
sp_car[sp_car$MaxLength > 60 & sp_car$MaxLength <= 80,]$ML_cat <- "60-80 cm"
sp_car[sp_car$MaxLength > 80 & sp_car$MaxLength <= 300,]$ML_cat <- "80-300 cm"

sp_car[sp_car$Water.column == "Demersal",]$Water.column <- "demersal"
sp_car[sp_car$Water.column == "pelagic non-site attached",]$Water.column <- "pelagic"
sp_car[sp_car$Water.column == "pelagic site attached",]$Water.column <- "pelagic"

sp_car[sp_car$Habitat == "Coral",]$Habitat <- "coral"

sp_car[sp_car$Trophic_guild_name == "Herbivores Microvores Detritivores",]$Trophic_guild_name <- "herbivores"

load("outputs/performance_model.Rdata")
load("data/new_derived_data/species_count.Rdata")

performance_bind_best <- performance_bind |> 
  dplyr::inner_join(best_models) |> 
  dplyr::filter(model == best_model) |> 
  dplyr::select(-c(model)) |> 
  dplyr::filter(best_model %in% c("glm", "gam", "rf", "sprf"))

performance_bind_best <- performance_bind |> 
  dplyr::inner_join(best_models) |> 
  dplyr::filter(model == best_model) |> 
  dplyr::select(-c(model, best_model))



# species <- species_covariates_importance_function(plot_data = bind_files) |> 
#   dplyr::select(!sd) |> 
#   dplyr::filter(!species_name %in% "Pseudolabrus luculentus") |> 
#   tidyr::pivot_wider(names_from = VAR,
#                      values_from = value) |> 
#   dplyr::group_by(species_name) |> 
#   dplyr::mutate(varmax = names(dplyr::across(c("ENV", "HAB", "HUM", "BIOT")))[
#     max.col(dplyr::across(c("ENV", "HAB", "HUM", "BIOT")), ties.method = "first")
#   ]) |> 
#   dplyr::inner_join(sp_car[,colnames(sp_car) %in% c("species_name", "MaxLength", "Trophic_guild_name", "Trophic.Level")]) |> 
#   tidyr::drop_na()

library(ggplot2)

varnames <- unique(bind_files$variable)
species <- bind_files |> 
  dplyr::inner_join(best_models) |> 
  dplyr::filter(fitted_model == best_model) |> 
  dplyr::filter(!species_name %in% "Pseudolabrus luculentus") |> 
  dplyr::select(!c(sd_dropout_loss, best_model, fitted_model)) |> 
  tidyr::pivot_wider(names_from = variable,
                     values_from = Dropout_loss) |> 
  dplyr::group_by(species_name) |>
  dplyr::mutate(varmax = names(dplyr::across(varnames))[
    max.col(dplyr::across(varnames), ties.method = "first")
  ]) |> 
  dplyr::mutate(varmax = dplyr::case_when(varmax %in% c("max_1year_analysed_sst",
                                                          "max_5year_degree_heating_week",
                                                          "mean_1year_nppv",
                                                          "mean_1year_so_mean",
                                                          "min_1year_analysed_sst",
                                                          "min_5year_ph") ~ "ENV",
                                          varmax %in% c("protection_status2",
                                                          "gdp",
                                                          "marine_ecosystem_dependency",
                                                          "gravtot2",
                                                          "n_fishing_vessels",
                                                          "neartt") ~ "HUM",
                                          varmax %in% c("Rock_500m",
                                                          "Sand_500m",
                                                          "coral",
                                                          "coral_algae_500m",
                                                          "depth",
                                                          "reef_extent") ~ "HAB",
                                          varmax %in% c("delta_biomass",
                                                          "diversity",
                                                          "max_trophic",
                                                          "mean_biomass",
                                                          "mean_trophic",
                                                          "n_trophic") ~ "BIOT")) |> 
  dplyr::inner_join(sp_car) |> 
  dplyr::inner_join(performance_bind_best) |> 
  dplyr::inner_join(sp_count)
  
PCA <- FactoMineR::PCA(species[,c(2:23)], graph = FALSE, scale.unit = TRUE)
factoextra::fviz_pca_biplot(PCA,
                            geom.ind = "point",
                            pointshape = 19,
                            col.ind = log10(species$count),
                            # col.var = colnames(species)[2:23],
                            repel = TRUE,
                            ggtheme = theme_minimal()) +
  scale_color_continuous(type = "viridis")
factoextra::fviz_pca_biplot(PCA,
                            geom.ind = "point",
                            pointshape = 19,
                            col.ind = species$best_model,
                            # col.var = colnames(species)[2:23],
                            repel = TRUE,
                            ggtheme = theme_minimal())


factoextra::fviz_pca_biplot(PCA,
                            geom.ind = "point",
                            pointshape = 19,
                            col.ind = species$Trophic.Level,
                            # col.var = colnames(species)[2:23],
                            repel = TRUE,
                            ggtheme = theme_minimal()) +
  scale_color_continuous(type = "viridis") +
  scale_color_manual(values = c("max_1year_analysed_sst" = pal_sp_trait[2],
                                "max_5year_degree_heating_week" = pal_sp_trait[2],
                                "mean_1year_nppv" = pal_sp_trait[2],
                                "mean_1year_so_mean" = pal_sp_trait[2],
                                "min_1year_analysed_sst" = pal_sp_trait[2],
                                "min_5year_ph" = pal_sp_trait[2],
                                "protection_status2" = pal_sp_trait[1],
                                "gdp" = pal_sp_trait[1],
                                "marine_ecosystem_dependency" = pal_sp_trait[1],
                                "gravtot2" = pal_sp_trait[1],
                                "n_fishing_vessels" = pal_sp_trait[1],
                                "neartt" = pal_sp_trait[1],
                                "Rock_500m" = pal_sp_trait[5],
                                "Sand_500m" = pal_sp_trait[5],
                                "coral" = pal_sp_trait[5],
                                "coral_algae_500m" = pal_sp_trait[5],
                                "depth" = pal_sp_trait[5],
                                "reef_extent" = pal_sp_trait[5],
                                "delta_biomass" = pal_sp_trait[3],
                                "diversity" = pal_sp_trait[3],
                                "max_trophic" = pal_sp_trait[3],
                                "mean_biomass" = pal_sp_trait[3],
                                "mean_trophic" = pal_sp_trait[3],
                                "n_trophic" = pal_sp_trait[3],
                                "ENV" = pal_sp_trait[2],
                                "HUM" = pal_sp_trait[1],
                                "HAB" = pal_sp_trait[13],
                                "BIOT" = pal_sp_trait[3])) +
  theme(legend.position = "none")
factoextra::fviz_pca_biplot(PCA,
                            geom.ind = "point",
                            # gradient.cols = viridis::viridis(10),
                            col.ind = species$varmax) +
  scale_color_manual(values = c("ENV" = pal_sp_trait[2],
                                "HUM" = pal_sp_trait[1],
                                "HAB" = pal_sp_trait[6],
                                "BIOT" = pal_sp_trait[3]))
famd <- FactoMineR::FAMD(species[,!colnames(species) %in% c("varmax", "species_name")], graph = FALSE)
factoextra::fviz_famd_ind(famd,
                          geom = c("point"),
                          repel = TRUE,
                          ggtheme = theme_minimal())
factoextra::fviz_famd_var(famd,
                          c("quali.var"),
                          geom = c("arrow", "text"),
                          repel = TRUE,
                          ggtheme = theme_minimal())


plot_max.length <- species_traits_function(plot_data = bind_files,
                                           data_trait = sp_car,
                                           trait = "ML_cat",
                                           color = pal_sp_trait,
                                           labs_title = "A. Maximum length (cm)",
                                           x_text_angle = NULL,
                                           legend.position = "none")

plot_water.column <- species_traits_function(plot_data = bind_files,
                                             data_trait = sp_car,
                                             trait = "Water.column",
                                             color = pal_sp_trait,
                                             labs_title = "B. Water column position ",
                                             x_text_angle = NULL,
                                             legend.position = "none")

plot_habitat <- species_traits_function(plot_data = bind_files,
                                        data_trait = sp_car,
                                        trait = "Habitat",
                                        color = pal_sp_trait,
                                        labs_title = "C. Habitat",
                                        x_text_angle = NULL,
                                        legend.position = "none")

plot_trophic <- species_traits_function(plot_data = bind_files,
                                        data_trait = sp_car,
                                        trait = "Trophic_guild_name",
                                        color = pal_sp_trait,
                                        labs_title = "D. Trophic classes",
                                        x_text_angle = 55,
                                        legend.position = "none")

load("data/new_raw_data/RLS_actinopterygii_data.Rdata")
phylo <- RLS_actinopterygii_data |> 
  dplyr::select(species_name, order, family)
phylo <- unique(phylo)
n_phylo <- bind_files |> 
  dplyr::inner_join(phylo, multiple = "first") |> 
  dplyr::select(species_name, family) |> 
  unique() |> 
  dplyr::group_by(family) |> 
  dplyr::summarise(n = dplyr::n())
phylo <- phylo |> 
  dplyr::inner_join(n_phylo, by = "family") |> 
  dplyr::filter(n >= 10) |> 
  dplyr::select(-n)
bind_files <- bind_files |> 
  dplyr::filter(species_name %in% phylo$species_name)

plot_familly <- species_traits_function(plot_data = bind_files,
                                        data_trait = phylo,
                                        trait = "family",
                                        color = pal_sp_trait,
                                        labs_title = "E. Family",
                                        x_text_angle = 55,
                                        legend.position = c(0.85, 0.8))

plot_species_traits <- plot_max.length / plot_water.column / plot_habitat / plot_trophic / plot_familly

ggsave("figures/plot_species_traits.pdf", plot_species_traits, height = 25, width = 19)



scaridae <- bind_files |> 
  dplyr::inner_join(phylo, multiple = "first") |> 
  dplyr::filter(family == "Scaridae") |> 
  dplyr::ungroup() |> 
  dplyr::select(-c(order, family))

scaridae <- species_covariates_importance_function(plot_data = scaridae)

scaridae <- scaridae |> 
  dplyr::mutate(var_reordered = tidytext::reorder_within(VAR, value, species_name))

plot_scaridae <- ggplot(scaridae) +
  geom_col(aes(x = var_reordered, y = value, fill = VAR)) +
  geom_errorbar(aes(x = var_reordered, y = value, ymin = value - sd, ymax = value + sd), width = .1, position = position_dodge(.9)) +
  theme_bw() +
  coord_flip() +
  facet_wrap(~species_name, scales = "free_y", ncol = 5) +
  tidytext::scale_x_reordered() +
  scale_fill_manual(values = c("ENV" = pal_sp_trait[2],
                               "HUM" = pal_sp_trait[1],
                               "HAB" = pal_sp_trait[13],
                               "BIOT" = pal_sp_trait[3])) +
  labs(y = "Relative importance (RMSE)", x = "", fill = "") +
  theme(
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white", color = NA),
    title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 15),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = "grey50",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) 

ggsave("figures/plot_scaridae.png", plot_scaridae, height = 10, width = 16)


pomacentridae <- bind_files |> 
  dplyr::inner_join(phylo, multiple = "first") |> 
  dplyr::filter(family == "Pomacentridae") |> 
  dplyr::ungroup() |> 
  dplyr::select(-c(order, family))

pomacentridae <- species_covariates_importance_function(plot_data = pomacentridae)

pomacentridae <- pomacentridae |> 
  dplyr::mutate(var_reordered = tidytext::reorder_within(VAR, value, species_name))

plot_pomacentridae <- ggplot(pomacentridae) +
  geom_col(aes(x = var_reordered, y = value, fill = VAR)) +
  geom_errorbar(aes(x = var_reordered, y = value, ymin = value - sd, ymax = value + sd), width = .1, position = position_dodge(.9)) +
  theme_bw() +
  coord_flip() +
  facet_wrap(~species_name, scales = "free_y", ncol = 5) +
  tidytext::scale_x_reordered() +
  scale_fill_manual(values = c("ENV" = pal_sp_trait[2],
                               "HUM" = pal_sp_trait[1],
                               "HAB" = pal_sp_trait[13],
                               "BIOT" = pal_sp_trait[3])) +
  labs(y = "Relative importance (RMSE)", x = "", fill = "") +
  theme(
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white", color = NA),
    title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 15),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = "grey50",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) 

ggsave("figures/plot_pomacentridae.png", plot_pomacentridae, height = 30, width = 16)
