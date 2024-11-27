load("data/new_raw_data/RLS_actinopterygii_data.Rdata")
load("data/new_derived_data/rls_covariates.RData")
load("data/new_raw_data/00_rls_surveys.Rdata")

rls_surveys$survey_id <- as.character(rls_surveys$survey_id)
rls_fish_data <- dplyr::inner_join(RLS_actinopterygii_data, rls_surveys)

rls_coral_fish <- rls_fish_data |>
  dplyr::filter(survey_id %in% rls_covariates$survey_id)

rls_coral_fish_mean_biomass <- rls_coral_fish |> 
  dplyr::group_by(survey_id, site_code, species_name, latitude, longitude, survey_date, depth) |> 
  dplyr::summarise(biomass = sum(biomass))

rls_coral_fish_mean_biomass_count <- rls_coral_fish_mean_biomass |>
  dplyr::group_by(species_name) |>
  dplyr::mutate(count = dplyr::n()) |>
  dplyr::filter(count >= 50) |>
  dplyr::select(survey_id, species_name, biomass, latitude, longitude)

rls_coral_fish_mean_biomass_count <- rls_coral_fish_mean_biomass_count |> 
  dplyr::inner_join(rls_surveys)

sp_car <- read.csv("data/new_raw_data/Traits_tropical_spp_1906.csv", header = TRUE)
sp_car <- sp_car |> 
  dplyr::rename(species_name = Species)

sp_car <- sp_car[which(is.na(sp_car$MaxLength) == FALSE),]
sp_car <- sp_car[which(is.na(sp_car$Trophic_guild_name) == FALSE),]

sp_car$ML_cat <- NA
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

phylo <- RLS_actinopterygii_data |> 
  dplyr::select(species_name, order, family)
phylo <- unique(phylo)

n_phylo <- rls_coral_fish_mean_biomass_count |> 
  dplyr::inner_join(phylo, multiple = "first") |> 
  dplyr::select(species_name, family) |> 
  unique() |> 
  dplyr::group_by(family) |> 
  dplyr::summarise(n = dplyr::n())
phylo <- phylo |> 
  dplyr::inner_join(n_phylo, by = "family") |> 
  dplyr::filter(n >= 10) |> 
  dplyr::select(-n)

rls_coral_fish_mean_biomass_count <- rls_coral_fish_mean_biomass_count |> 
  dplyr::inner_join(sp_car) |> 
  dplyr::inner_join(rls_covariates) |> 
  dplyr::inner_join(phylo)

realm_selec <- rls_coral_fish_mean_biomass_count |> 
  dplyr::ungroup() |>
  dplyr::select(survey_id, realm) |> 
  unique() |> 
  dplyr::group_by(realm) |> 
  dplyr::summarise(n = dplyr::n()) |> 
  dplyr::filter(n > 100)

rls_coral_fish_mean_biomass_count <- rls_coral_fish_mean_biomass_count |> 
  dplyr::ungroup() |> 
  dplyr::filter(realm %in% realm_selec$realm)

mean_global <- mean(rls_coral_fish_mean_biomass_count$biomass)
med_global <- median(rls_coral_fish_mean_biomass_count$biomass)
min_global <- min(rls_coral_fish_mean_biomass_count$biomass)
max_global <- max(rls_coral_fish_mean_biomass_count$biomass)

stat_biomass_eco <- rls_coral_fish_mean_biomass_count |> 
  dplyr::group_by(ecoregion) |> 
  dplyr::summarise(med_biomass = median(biomass),
                   sd_biomass = sd(biomass),
                   max_biomass = max(biomass),
                   min_biomass = min(biomass),
                   inf_med = ifelse(med_biomass < med_global, TRUE, FALSE))
stat_biomass_realm <- rls_coral_fish_mean_biomass_count |> 
  dplyr::group_by(realm) |> 
  dplyr::summarise(med_biomass = median(biomass),
                   max_biomass = max(biomass),
                   q95_biomass = quantile(biomass, probs = 0.95),
                   q75_biomass = quantile(biomass, probs = 0.75),
                   inf_med = ifelse(med_biomass < med_global, TRUE, FALSE))
stat_biomass_realm_troph <- rls_coral_fish_mean_biomass_count |> 
  dplyr::group_by(realm, Trophic_guild_name) |> 
  dplyr::summarise(med_biomass = median(biomass),
                   max_biomass = max(biomass),
                   q95_biomass = quantile(biomass, probs = 0.95),
                   q75_biomass = quantile(biomass, probs = 0.75),
                   inf_med = ifelse(med_biomass < med_global, TRUE, FALSE))
stat_biomass_realm_fam <- rls_coral_fish_mean_biomass_count |> 
  dplyr::group_by(realm, family) |> 
  dplyr::summarise(med_biomass = median(biomass),
                   max_biomass = max(biomass),
                   q95_biomass = quantile(biomass, probs = 0.95),
                   q75_biomass = quantile(biomass, probs = 0.75),
                   inf_med = ifelse(med_biomass < med_global, TRUE, FALSE))

rls_coral_fish_mean_biomass_count$biomass <- log10(rls_coral_fish_mean_biomass_count$biomass + 1)

library(ggplot2)

pal_contribution <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1"), RColorBrewer::brewer.pal(n = 5, name = "Set3"))

troph_realm <- rls_coral_fish_mean_biomass_count |> 
  dplyr::mutate(troph_reordered = tidytext::reorder_within(Trophic_guild_name, biomass, realm, fun = median))

biomass_realm_troph_plot <- ggplot(troph_realm, aes(x = troph_reordered, y = biomass)) +
  geom_point(position = "jitter", alpha = 0.2, size = 0.5) +
  geom_boxplot(aes(fill = Trophic_guild_name), alpha = 0.9, outlier.shape = NA) +
  facet_wrap(~ realm, scales = "free_x") +
  theme_bw() +
  scale_fill_manual(values = pal_contribution) +
  tidytext::scale_x_reordered() +
  labs(y = "log10(Biomass + 1)", x = "", fill = "Trophic group", title = "A.") +
  theme(title = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        strip.background = element_blank())

fam_realm <- rls_coral_fish_mean_biomass_count |> 
  dplyr::mutate(fam_reordered = tidytext::reorder_within(family, biomass, realm, fun = median))

biomass_realm_fam_plot <- ggplot(fam_realm, aes(x = fam_reordered, y = biomass)) +
  geom_point(position = "jitter", alpha = 0.2, size = 0.5) +
  geom_boxplot(aes(fill = family), alpha = 0.9, outlier.shape = NA) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ realm, scales = "free_x") +
  theme_bw() +
  scale_fill_manual(values = pal_contribution) +
  tidytext::scale_x_reordered() +
  labs(y = "log10(Biomass + 1)", x = "", fill = "Family", title = "B.") +
  theme(title = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        strip.background = element_blank())

library(patchwork)

troph_fam_biomass_realm_plot <- biomass_realm_troph_plot / biomass_realm_fam_plot

ggsave("figures/troph_fam_biomass_realm_plot.png", troph_fam_biomass_realm_plot, width = 12, height = 16)

biomass_realm_ecoregion_plot <- rls_coral_fish_mean_biomass_count |> 
  dplyr::mutate(ecoregion = reorder(ecoregion, biomass, FUN = median)) |>
  ggplot(aes(x = ecoregion, y = biomass, fill = realm)) +
  geom_point(position = "jitter", alpha = 0.2, size = 0.5) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA) +
  geom_abline(slope = 0, intercept = median(rls_coral_fish_mean_biomass_count$biomass), color = "red", linewidth = 1.3, linetype = 2) +
  scale_fill_manual(values = pal_contribution) +
  theme_classic() +
  labs(y = "log10(Biomass + 1)", x = "Ecoregion", fill = "Realm") +
  theme(title = element_text(size = 18),
        axis.text.x = element_text(size = 10, angle = 55, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20))
# ggsave("figures/biomass_realm_ecoregion_plot.pdf", biomass_realm_ecoregion_plot, width = 20, height = 8)  
ggsave("figures/biomass_realm_ecoregion_plot.png", biomass_realm_ecoregion_plot, width = 20, height = 8)  


