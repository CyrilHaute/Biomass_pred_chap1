# source functions ----

library(patchwork)
library(ggplot2)

source("R/05_contributions_figures_functions.R")
source("R/05_load_realm_contribution_function.R")

pal_contribution <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1"), PNWColors::pnw_palette("Bay", 6, type = "continuous"))

load("outputs/best_models.Rdata")

env_var <- c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_nppv", "mean_1year_so_mean", "min_1year_analysed_sst", "min_5year_ph")
hum_var <- c("protection_status2", "gdp", "gravtot2", "n_fishing_vessels", "neartt", "marine_ecosystem_dependency")
hab_var <- c("Rock_500m", "Sand_500m", "coral", "coral_algae_500m", "depth", "reef_extent")

#################### Global contribution ####################


#### Covariates contribution plot ####

glm <- list.files("outputs/glm_contribution", full.names = TRUE)
glm <- lapply(1:length(glm), function(i) {
  
  load(glm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
glm <- do.call(rbind, glm)
glm_inf <- list.files("outputs/glm_contribution_benthos_infered", full.names = TRUE)
glm_inf <- lapply(1:length(glm_inf), function(i) {
  
  load(glm_inf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
glm_inf <- do.call(rbind, glm_inf)

gam <- list.files("outputs/gam_contribution", full.names = TRUE)
gam <- lapply(1:length(gam), function(i) {
  
  load(gam[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gam <- do.call(rbind, gam)
gam_inf <- list.files("outputs/gam_contribution_benthos_infered", full.names = TRUE)
gam_inf <- lapply(1:length(gam_inf), function(i) {
  
  load(gam_inf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gam_inf <- do.call(rbind, gam_inf)

spamm <- list.files("outputs/spamm_contribution", full.names = TRUE)
spamm <- lapply(1:length(spamm), function(i) {
  
  load(spamm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
spamm <- do.call(rbind, spamm)
spamm_inf <- list.files("outputs/spamm_contribution_benthos_infered", full.names = TRUE)
spamm_inf <- lapply(1:length(spamm_inf), function(i) {
  
  load(spamm_inf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
spamm_inf <- do.call(rbind, spamm_inf)

rf <- list.files("outputs/rf_contribution", full.names = TRUE)
rf <- lapply(1:length(rf), function(i) {
  
  load(rf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
rf <- do.call(rbind, rf)
rf_inf <- list.files("outputs/rf_contribution_benthos_infered", full.names = TRUE)
rf_inf <- lapply(1:length(rf_inf), function(i) {
  
  load(rf_inf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
rf_inf <- do.call(rbind, rf_inf)

sprf <- list.files("outputs/sprf_contribution", full.names = TRUE)
sprf <- lapply(1:length(sprf), function(i) {
  
  load(sprf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
sprf <- do.call(rbind, sprf)
colnames(sprf)[colnames(sprf) %in% c("global_dropout_loss", "global_sd_dropout_loss")] <- c("Dropout_loss", "sd_dropout_loss")
sprf_inf <- list.files("outputs/sprf_contribution_benthos_infered", full.names = TRUE)
sprf_inf <- lapply(1:length(sprf_inf), function(i) {
  
  load(sprf_inf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
sprf_inf <- do.call(rbind, sprf_inf)
colnames(sprf_inf)[colnames(sprf_inf) %in% c("global_dropout_loss", "global_sd_dropout_loss")] <- c("Dropout_loss", "sd_dropout_loss")

gbm <- list.files("outputs/brt_contribution", full.names = TRUE)
gbm <- lapply(1:length(gbm), function(i) {
  
  load(gbm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gbm <- do.call(rbind, gbm)
gbm_inf <- list.files("outputs/brt_contribution_benthos_infered", full.names = TRUE)
gbm_inf <- lapply(1:length(gbm_inf), function(i) {
  
  load(gbm_inf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gbm_inf <- do.call(rbind, gbm_inf)

bind_files <- list(glm, gam, spamm, sprf, rf, gbm)
bind_files <- purrr::reduce(bind_files, dplyr::full_join)

bind_files_inf <- list(glm_inf, gam_inf, spamm_inf, sprf_inf, rf_inf, gbm_inf)
bind_files_inf <- purrr::reduce(bind_files_inf, dplyr::full_join)

covariates_importance_all <- covariates_importance_all_function(plot_data = bind_files,
                                                                title = "A.",
                                                                legend.position = "none",
                                                                title.size = 18,
                                                                axis.text.x = 15,
                                                                axis.text.y = 17,
                                                                axis.title = 21,
                                                                legend.text = 15,
                                                                strip.text.x = 20,
                                                                strip.text.y = 20,
                                                                fill = "")

merged_covariates_importance_all <- merged_covariates_importance_all_function(plot_data = bind_files,
                                                                              title = "B.",
                                                                              legend.position = "none",
                                                                              title.size = 18,
                                                                              axis.text.x = 15,
                                                                              axis.text.y = 17,
                                                                              axis.title = 21,
                                                                              legend.text = 15,
                                                                              strip.text.x = 20,
                                                                              strip.text.y = 20,
                                                                              geom.text.size = 5,
                                                                              fill = "")

covariates_importance_all_and_merged <- covariates_importance_all + merged_covariates_importance_all

ggsave("figures/covariates_importance_all_and_merged2.png", covariates_importance_all_and_merged, height = 7, width = 15)

covariates_importance_GLM <- plot_covariates_importance_function(plot_data = covariates_importance_function(plot_data = glm),
                                                                 color = pal_contribution,
                                                                 labs_y = "",
                                                                 labs_fill = "",
                                                                 ylim = c(0, 0.4),
                                                                 legend.position = "none")

covariates_importance_GAM <- plot_covariates_importance_function(plot_data = covariates_importance_function(plot_data = gam),
                                                                 color = pal_contribution,
                                                                 labs_y = "",
                                                                 labs_fill = "",
                                                                 ylim = c(0, 0.4),
                                                                 legend.position = "none")

covariates_importance_SPAMM <- plot_covariates_importance_function(plot_data = covariates_importance_function(plot_data = spamm),
                                                                   color = pal_contribution,
                                                                   labs_y = "",
                                                                   labs_fill = "",
                                                                   ylim = c(0, 0.4),
                                                                   legend.position = "none")

covariates_importance_RF <- plot_covariates_importance_function(plot_data = covariates_importance_function(plot_data = rf),
                                                                color = pal_contribution,
                                                                labs_y = "",
                                                                labs_fill = "",
                                                                ylim = c(0, 0.4),
                                                                legend.position = "none")

covariates_importance_GBM <- plot_covariates_importance_function(plot_data = covariates_importance_function(plot_data = gbm),
                                                                 color = pal_contribution,
                                                                 labs_y = "Relative importance (RMSE)",
                                                                 labs_fill = "",
                                                                 ylim = c(0, 0.4),
                                                                 legend.position = "none")

covariates_importance_SPRF <- plot_covariates_importance_function(plot_data = covariates_importance_function(plot_data = sprf),
                                                                  color = pal_contribution,
                                                                  labs_y = "Relative importance (RMSE)",
                                                                  labs_fill = "",
                                                                  ylim = c(0, 0.4),
                                                                  legend.position = c(0.75, 0.16))

covariates_importance_all <- (covariates_importance_GLM + covariates_importance_GAM) / (covariates_importance_SPAMM + covariates_importance_RF) / (covariates_importance_GBM + covariates_importance_SPRF)

ggsave("figures/covariates_importance_all.png", covariates_importance_all, height = 16, width = 13)

merged_covariates_importance_GLM <- plot_merged_covariates_importance_function(plot_data =  merged_covariates_importance_function(plot_data = glm,
                                                                                                                                  fitted_model = "glm"),
                                                                               color = pal_contribution,
                                                                               labs_y = "",
                                                                               labs_fill = "",
                                                                               legend.position = "none",
                                                                               mul = 2)

merged_covariates_importance_GAM <- plot_merged_covariates_importance_function(plot_data =  merged_covariates_importance_function(plot_data = gam,
                                                                                                                                  fitted_model = "gam"),
                                                                               color = pal_contribution,
                                                                               labs_y = "",
                                                                               labs_fill = "",
                                                                               legend.position = "none",
                                                                               mul = 2)

merged_covariates_importance_SPAMM <- plot_merged_covariates_importance_function(merged_covariates_importance_function(plot_data = spamm,
                                                                                                                       fitted_model = "spamm"),
                                                                                 color = pal_contribution,
                                                                                 labs_y = "",
                                                                                 labs_fill = "",
                                                                                 legend.position = "none",
                                                                                 mul = 2)

merged_covariates_importance_RF <- plot_merged_covariates_importance_function(merged_covariates_importance_function(plot_data = rf,
                                                                                                                    fitted_model = "rf"),
                                                                              color = pal_contribution,
                                                                              labs_y = "",
                                                                              labs_fill = "",
                                                                              legend.position = "none",
                                                                              mul = 3)

merged_covariates_importance_GBM <- plot_merged_covariates_importance_function(merged_covariates_importance_function(plot_data = gbm,
                                                                                                                     fitted_model = "gbm"),
                                                                               color = pal_contribution,
                                                                               labs_y = "Relative importance (RMSE)",
                                                                               labs_fill = "",
                                                                               legend.position = "none",
                                                                               mul = 3)

merged_covariates_importance_SPRF <- plot_merged_covariates_importance_function(merged_covariates_importance_function(plot_data = sprf,
                                                                                                                      fitted_model = "sprf"),
                                                                                color = pal_contribution,
                                                                                labs_y = "Relative importance (RMSE)",
                                                                                labs_fill = "",
                                                                                legend.position = c(0.8, 0.18),
                                                                                mul = 3)

merged_covariates_importance <- (merged_covariates_importance_GLM + merged_covariates_importance_GAM) / (merged_covariates_importance_SPAMM + merged_covariates_importance_RF) / (merged_covariates_importance_GBM + merged_covariates_importance_SPRF)

ggsave("figures/merged_covariates_importance.png", merged_covariates_importance, height = 15, width = 11)



#################### Per realm contribution ####################

glm_realm <- load_realm_cont_function(files_path = "outputs/glm_biomass_contribution_realm")
rf_realm <- load_realm_cont_function(files_path = "outputs/rf_biomass_contribution_realm")
gam_realm <- load_realm_cont_function(files_path = "outputs/gam_biomass_contribution_realm")
gbm_realm <- load_realm_cont_function(files_path = "outputs/gbm_biomass_contribution_realm")
spamm_realm <- load_realm_cont_function(files_path = "outputs/spamm_biomass_contribution_realm")
sprf_realm <- load_realm_cont_function(files_path = "outputs/sprf_biomass_contribution_realm")
sprf_realm <- do.call(rbind, sprf_realm)
colnames(sprf_realm)[colnames(sprf_realm) %in% c("global_dropout_loss", "global_sd_dropout_loss")] <- c("Dropout_loss", "sd_dropout_loss")
sprf_realm <- sprf_realm |> 
  dplyr::group_split(realm)

contribution_realm_data <- list(glm_realm, gam_realm, spamm_realm, rf_realm, gbm_realm, sprf_realm)
contribution_realm_data <- list(do.call(rbind, lapply(contribution_realm_data, '[[', 1)),
                                do.call(rbind, lapply(contribution_realm_data, '[[', 2)),
                                do.call(rbind, lapply(contribution_realm_data, '[[', 3)),
                                do.call(rbind, lapply(contribution_realm_data, '[[', 4)),
                                do.call(rbind, lapply(contribution_realm_data, '[[', 5)))
contribution_realm_data[[1]]$realm <- paste0(contribution_realm_data[[1]]$realm, " (n = 3527)")
contribution_realm_data[[2]]$realm <- paste0(contribution_realm_data[[2]]$realm, " (n = 205)")
contribution_realm_data[[3]]$realm <- paste0(contribution_realm_data[[3]]$realm, " (n = 186)")
contribution_realm_data[[4]]$realm <- paste0(contribution_realm_data[[4]]$realm, " (n = 246)")
contribution_realm_data[[5]]$realm <- paste0(contribution_realm_data[[5]]$realm, " (n = 338)")

merged_covariates_importance_all <- lapply(1:length(contribution_realm_data), function(i) {
  
  merged_covariates_importance_all_function(plot_data = contribution_realm_data[[i]],
                                            title = stringr::str_replace_all(unique(contribution_realm_data[[i]]$realm), c("_" = " ")),
                                            legend.position = "none",
                                            title.size = 8,
                                            axis.text.x = 9,
                                            axis.text.y = 9,
                                            axis.title = 10,
                                            legend.text = 15,
                                            strip.text.x = 9,
                                            strip.text.y = 9,
                                            geom.text.size = 3,
                                            fill = "")

})

covariates_importance_all <- lapply(1:length(contribution_realm_data), function(i) {
  
  covariates_importance_all_function(plot_data = contribution_realm_data[[i]],
                                     title = stringr::str_replace_all(unique(contribution_realm_data[[i]]$realm), "_", " "),
                                     legend.position = "none",
                                     title.size = 10,
                                     axis.text.x = 11,
                                     axis.text.y = 11,
                                     axis.title = 13,
                                     legend.text = 15,
                                     strip.text.x = 13,
                                     strip.text.y = 9,
                                     fill = "Covariate categories")
  
})

world <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")

rls_map <- ggplot(data = world) +
  geom_sf(fill = "white", color = "grey", size = 0.03) +
  scale_x_continuous(limits = c(-180, 180)) +
  scale_y_continuous(limits = c(-60, 60)) +
  theme_classic() +
  labs(title = "") +
  theme(legend.position = "none",
        axis.text = element_text(size = 15),
        title = element_text(size = 13))

map_contribution_merged_realm <- rls_map + 
  patchwork::inset_element(merged_covariates_importance_all[[3]], left = 0.35, bottom = 0.53, right = 0.52, top = 1) +
  patchwork::inset_element(merged_covariates_importance_all[[4]], left = 0.34, bottom = 0.07, right = 0.51, top = 0.52) +
  patchwork::inset_element(merged_covariates_importance_all[[1]], left = 0.80, bottom = 0.26, right = 0.96, top = 0.74) +
  patchwork::inset_element(merged_covariates_importance_all[[2]], left = 0.005, bottom = 0.14, right = 0.17, top = 0.61) +
  patchwork::inset_element(merged_covariates_importance_all[[5]], left = 0.16, bottom = 0.3, right = 0.32, top = 0.77) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")

covariates_importance_all_bind <- (covariates_importance_all[[1]] + covariates_importance_all[[2]] + covariates_importance_all[[3]] + covariates_importance_all[[4]] + covariates_importance_all[[5]])

ggplot2::ggsave("figures/covariates_importance_map.png", map_contribution_merged_realm, height = 6, width = 15)

ggplot2::ggsave("figures/contribution_realm.png", covariates_importance_all_bind, height = 10, width = 20)

##################### Contribution per species #####################

# Load species traits
sp_car <- read.csv("data/new_raw_data/Traits_tropical_spp_1906.csv", header = TRUE) |> 
  dplyr::rename(species_name = Species) |> 
  dplyr::filter(species_name %in% unique(bind_files$species_name))

sp_car <- sp_car[which(is.na(sp_car$Trophic_guild_name) == FALSE),]

sp_car[sp_car$Trophic_guild_name == "Herbivores Microvores Detritivores",]$Trophic_guild_name <- "Herbivores"
sp_car[sp_car$Trophic_guild_name == "planktivore",]$Trophic_guild_name <- "Planktivores"
sp_car[sp_car$Trophic_guild_name == "piscivore",]$Trophic_guild_name <- "Piscivores"
sp_car[sp_car$Trophic_guild_name == "microinvertivore",]$Trophic_guild_name <- "Microinvertivores"
sp_car[sp_car$Trophic_guild_name == "macroinvertivore",]$Trophic_guild_name <- "Macroinvertivores"
sp_car[sp_car$Trophic_guild_name == "sessile invertivores",]$Trophic_guild_name <- "Sessile invertivores"
sp_car[sp_car$Trophic_guild_name == "corallivore",]$Trophic_guild_name <- "Corallivores"
sp_car[sp_car$Trophic_guild_name == "crustacivore",]$Trophic_guild_name <- "Crustacivores"

trophic_group <- bind_files |> 
  dplyr::inner_join(sp_car[colnames(sp_car) %in% c("species_name", "Trophic_guild_name")])


trophic_group <- bind_files |> 
  dplyr::inner_join(sp_car[colnames(sp_car) %in% c("species_name", "Trophic_guild_name")]) |> 
  dplyr::mutate(VAR = dplyr::case_when(variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week",
                                                       "mean_1year_nppv", "mean_1year_so_mean", "min_1year_analysed_sst",
                                                       "min_5year_ph") ~ "ENV",
                                       variable %in% c("Rock_500m", "Sand_500m", "coral", "coral_algae_500m", "depth",
                                                       "reef_extent") ~ "HAB",
                                       variable %in% c("gdp", "gravtot2", "marine_ecosystem_dependency", "n_fishing_vessels",            
                                                       "neartt", "protection_status2") ~ "HUM"))
trophic_group <- trophic_group |> 
  dplyr::group_by(species_name, Trophic_guild_name, VAR) |> 
  dplyr::summarise(value = median(Dropout_loss),
                   sd = median(sd_dropout_loss))

trophic_group <- trophic_group |> 
  dplyr::filter(!species_name == "Chrysiptera notialis")

plot_trophic_group <- ggplot(trophic_group, aes(x = Trophic_guild_name, y = value, fill = VAR)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values = c("ENV" = pal_contribution[2],
                               "HUM" = pal_contribution[1],
                               "HAB" = pal_contribution[13])) +
  scale_y_continuous(limits = c(0, 0.15)) + 
  labs(y = "Relative importance (RMSE)", x = "", fill = "") +
  theme(
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white", color = NA),
    title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = "grey50",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

# Test ANOVA
anova_result <- aov(value ~ interaction(Trophic_guild_name, VAR), data = trophic_group)

# Test post-hoc Tukey
tukey_result <- TukeyHSD(anova_result)
tukey_table <- as.data.frame(tukey_result$`interaction(Trophic_guild_name, VAR)`)

diff <- tukey_table$`p adj` < 0.05
names(diff) <- rownames(tukey_table)

# Conversion en lettres avec multcompView
significance_letters <- multcompView::multcompLetters(diff,
                                                      Letters = letters)$Letters

y <- trophic_group |>
  dplyr::group_by(Trophic_guild_name, VAR) |>
  dplyr::summarise(y = quantile(value, probs = 0.75) + 0.05*quantile(value, probs = 0.75)) |>
  dplyr::ungroup() |>
  dplyr::mutate(interaction = paste0(Trophic_guild_name, ".", VAR))

label_data <- data.frame(
  interaction = names(significance_letters),
  Letter = significance_letters
) |>
  dplyr::inner_join(y)

label_data <- label_data |> 
  dplyr::arrange(interaction)

label_data$x <- gsub("\\..*", "", label_data$interaction)

label_data$x <- as.numeric(as.factor(label_data$x))

label_data <- label_data |> 
  dplyr::mutate(x = ifelse(VAR == "ENV", x - 0.35,
                           ifelse(VAR == "HUM", x + 0.35, x + 0.08)))

plot_trophic_group_letter <- plot_trophic_group + geom_text(data = label_data,
                                                            aes(x = x, y = y, label = Letter),
                                                            inherit.aes = FALSE,
                                                            size = 5)

ggsave("figures/plot_trophic_group_letter.png", plot_trophic_group_letter, height = 8, width = 18)



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
row_to_remove <- c(rownames(phylo[phylo$species_name == "Chlorurus microrhinos" & phylo$family == "Labridae",]),
                   rownames(phylo[phylo$species_name == "Oxycheilinus digramma" & phylo$order == "Eupercaria incertae sedis",]),
                   rownames(phylo[phylo$species_name == "Stegastes lacrymatus" & phylo$order == "Ovalentaria incertae sedis",]))
phylo <- phylo[-as.numeric(row_to_remove),]

phylo_bind <- bind_files |> 
  dplyr::inner_join(phylo)


phylo_bind <- phylo_bind |> 
  dplyr::mutate(VAR = dplyr::case_when(variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week",
                                                       "mean_1year_nppv", "mean_1year_so_mean", "min_1year_analysed_sst",
                                                       "min_5year_ph") ~ "ENV",
                                       variable %in% c("Rock_500m", "Sand_500m", "coral", "coral_algae_500m", "depth",
                                                       "reef_extent") ~ "HAB",
                                       variable %in% c("gdp", "gravtot2", "marine_ecosystem_dependency", "n_fishing_vessels",            
                                                       "neartt", "protection_status2") ~ "HUM"))
phylo_bind <- phylo_bind |> 
  dplyr::group_by(species_name, family, VAR) |> 
  dplyr::summarise(value = median(Dropout_loss),
                   sd = median(sd_dropout_loss))

phylo_bind <- phylo_bind |> 
  dplyr::filter(!species_name == "Chrysiptera notialis")

plot_family <- ggplot(phylo_bind, aes(x = family, y = value, fill = VAR)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values = c("ENV" = pal_contribution[2],
                               "HUM" = pal_contribution[1],
                               "HAB" = pal_contribution[13])) +
  scale_y_continuous(limits = c(0, 0.15)) + 
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

# Test ANOVA
anova_result <- aov(value ~ interaction(family, VAR), data = phylo_bind)

# Test post-hoc Tukey
tukey_result <- TukeyHSD(anova_result)
tukey_table <- as.data.frame(tukey_result$`interaction(family, VAR)`)

diff <- tukey_table$`p adj` < 0.05
names(diff) <- rownames(tukey_table)

# Conversion en lettres avec multcompView
significance_letters <- multcompView::multcompLetters(diff,
                                                      Letters = letters)$Letters

y <- phylo_bind |>
  dplyr::group_by(family, VAR) |>
  dplyr::summarise(y = quantile(value, probs = 0.75) + 0.05*quantile(value, probs = 0.75)) |>
  dplyr::ungroup() |>
  dplyr::mutate(interaction = paste0(family, ".", VAR))

label_data <- data.frame(
  interaction = names(significance_letters),
  Letter = significance_letters
) |>
  dplyr::inner_join(y)

label_data <- label_data |> 
  dplyr::arrange(interaction)

label_data$x <- gsub("\\..*", "", label_data$interaction)

label_data$x <- as.numeric(as.factor(label_data$x))

label_data <- label_data |> 
  dplyr::mutate(x = ifelse(VAR == "ENV", x - 0.35,
                           ifelse(VAR == "HUM", x + 0.35, x + 0.08)))

plot_family_letter <- plot_family + geom_text(data = label_data,
                                              aes(x = x, y = y, label = Letter),
                                              inherit.aes = FALSE,
                                              size = 5)

ggsave("figures/plot_trophic_group_letter.png", plot_trophic_group_letter, height = 8, width = 18)


##################### Test inferred benthos #####################

covariates_importance_all <- covariates_importance_all_function(plot_data = bind_files,
                                                                title = "With benthos inferred (n = 4684)",
                                                                legend.position = "none",
                                                                title.size = 16,
                                                                axis.text.x = 15,
                                                                axis.text.y = 17,
                                                                axis.title = 21,
                                                                legend.text = 15,
                                                                strip.text.x = 20,
                                                                strip.text.y = 20,
                                                                fill = "")
covariates_importance_all_inf <- covariates_importance_all_function(plot_data = bind_files_inf,
                                                                    title = "Without benthos inferred (n = 3363)",
                                                                    legend.position = "none",
                                                                    title.size = 16,
                                                                    axis.text.x = 15,
                                                                    axis.text.y = 17,
                                                                    axis.title = 21,
                                                                    legend.text = 15,
                                                                    strip.text.x = 20,
                                                                    strip.text.y = 20,
                                                                    fill = "")
test_benthos_infered <- covariates_importance_all + covariates_importance_all_inf
ggsave("figures/test_benthos_inferred.png", test_benthos_infered, height = 7, width = 20)



# load("data/new_raw_data/RLS_actinopterygii_data.Rdata")
# phylo <- RLS_actinopterygii_data |> 
#   dplyr::select(species_name, order, family)
# phylo <- unique(phylo)
# n_phylo <- bind_files |> 
#   dplyr::inner_join(phylo, multiple = "first") |> 
#   dplyr::select(species_name, family) |> 
#   unique() |> 
#   dplyr::group_by(family) |> 
#   dplyr::summarise(n = dplyr::n())
# phylo <- phylo |> 
#   dplyr::inner_join(n_phylo, by = "family") |> 
#   dplyr::filter(n >= 10) |> 
#   dplyr::select(-n)
# 
# scaridae <- bind_files |>
#   dplyr::inner_join(phylo, multiple = "first") |>
#   dplyr::filter(family == "Scaridae") |>
#   dplyr::ungroup() |>
#   dplyr::select(-c(order, family))
# 
# scaridae <- species_covariates_importance_function(plot_data = scaridae,
#                                                    group_by = "species_name")
# 
# scaridae <- scaridae |>
#   dplyr::mutate(var_reordered = tidytext::reorder_within(VAR, value, species_name))
# 
# plot_scaridae <- ggplot(scaridae) +
#   geom_col(aes(x = var_reordered, y = value, fill = VAR)) +
#   geom_errorbar(aes(x = var_reordered, y = value, ymin = value - sd, ymax = value + sd), width = .1, position = position_dodge(.9)) +
#   theme_bw() +
#   coord_flip() +
#   facet_wrap(~species_name, scales = "free_y", ncol = 4) +
#   tidytext::scale_x_reordered() +
#   scale_fill_manual(values = c("ENV" = pal_contribution[2],
#                                "HUM" = pal_contribution[1],
#                                "HAB" = pal_contribution[13])) +
#   labs(y = "Relative importance (RMSE)", x = "", fill = "") +
#   theme(
#     legend.direction = "vertical",
#     legend.background = element_rect(fill = "white"),
#     legend.key = element_rect(fill = "white", color = NA),
#     title = element_text(size = 15),
#     axis.text = element_text(size = 15),
#     axis.text.x = element_text(size = 15),
#     axis.text.y = element_text(size = 15),
#     axis.title = element_text(size = 15),
#     legend.text = element_text(size = 15),
#     legend.title = element_text(size = 15),
#     strip.text.x = element_text(size = 12, face = "italic"),
#     strip.text.y = element_text(size = 15),
#     strip.background = element_blank(),
#     panel.background = element_rect(fill = "white", colour = "grey50",
#                                     size = 1, linetype = "solid"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank())
# 
# ggsave("figures/plot_scaridae2.png", plot_scaridae, height = 10, width = 13)
