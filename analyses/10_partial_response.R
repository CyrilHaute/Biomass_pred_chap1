source("R/05_load_realm_contribution_function.R")

partial_realm_rf <- load_realm_partial_function(files_path = "outputs/realm_partial_rf")
partial_realm_sprf <- load_realm_partial_function(files_path = "outputs/realm_partial_sprf")
partial_realm_gbm <- load_realm_partial_function(files_path = "outputs/realm_partial_gbm")
partial_realm_glm <- load_realm_partial_function(files_path = "outputs/realm_partial_glm")
partial_realm_gam <- load_realm_partial_function(files_path = "outputs/realm_partial_gam")
partial_realm_list <- c(partial_realm_rf, partial_realm_sprf, partial_realm_gbm, partial_realm_glm, partial_realm_gam)

partial_realm <- lapply(1:length(partial_realm_list), function(i) {
  
  realm_i <- partial_realm_list[[i]]
  
  partial_cov <- lapply(1:18, function(j) {
    
    partial_cov_j <- lapply(realm_i, '[[', j)
    partial_cov_j <- lapply(1:length(partial_cov_j), function(k) {
      
      data_k <- partial_cov_j[[k]]
      data_k$id <- 1:nrow(data_k)
      data_k
      
    })
    partial_cov_j <- do.call(rbind, partial_cov_j)
    partial_cov_j <- partial_cov_j |> 
      dplyr::mutate(var = colnames(partial_cov_j)[1])
    colnames(partial_cov_j)[1] <- "values"
    partial_cov_j
    
  })
  partial_cov <- do.call(rbind, partial_cov)
  
})
partial_realm <- do.call(rbind, partial_realm)

partial_realm <- partial_realm |> 
  dplyr::mutate(var_type = dplyr::case_when(var %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_nppv", "mean_1year_so_mean", "min_1year_analysed_sst", "min_5year_ph") ~ "ENV",
                                            var %in% c("protection_status2", "gdp", "gravtot2", "n_fishing_vessels", "neartt", "marine_ecosystem_dependency") ~ "HUM",
                                            var %in% c("Rock_500m", "Sand_500m", "coral", "coral_algae_500m", "depth", "reef_extent") ~ "HAB"))

partial_realm[partial_realm$var == "max_1year_analysed_sst",]$var <- "SST max"
partial_realm[partial_realm$var == "min_1year_analysed_sst",]$var <- "SST min"
partial_realm[partial_realm$var == "max_5year_degree_heating_week",]$var <- "DHW max"
partial_realm[partial_realm$var == "mean_1year_nppv",]$var <- "NPP mean"
partial_realm[partial_realm$var == "mean_1year_so_mean",]$var <- "SSS mean"
partial_realm[partial_realm$var == "min_5year_ph",]$var <- "pH min"

partial_realm[partial_realm$var == "protection_status2",]$var <- "MPA"
partial_realm[partial_realm$var == "gdp",]$var <- "GDP"
partial_realm[partial_realm$var == "gravtot2",]$var <- "Gravity"
partial_realm[partial_realm$var == "n_fishing_vessels",]$var <- "Fishing"
partial_realm[partial_realm$var == "neartt",]$var <- "Neartt"
partial_realm[partial_realm$var == "marine_ecosystem_dependency",]$var <- "MED"

partial_realm[partial_realm$var == "Rock_500m",]$var <- "Rock (%)"
partial_realm[partial_realm$var == "Sand_500m",]$var <- "Sand (%)"
partial_realm[partial_realm$var == "coral_algae_500m",]$var <- "Coral/Algae (%)"
partial_realm[partial_realm$var == "coral",]$var <- "Coral (RLS)"
partial_realm[partial_realm$var == "depth",]$var <- "Depth"
partial_realm[partial_realm$var == "reef_extent",]$var <- "Reef extent"

library(ggplot2)
library(patchwork)

partial_realm$values <- as.numeric(partial_realm$values)

pal_contribution <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1"), PNWColors::pnw_palette("Bay", 6, type = "continuous"))

partial_realm_mean <- partial_realm |> 
  dplyr::group_by(var, id, realm) |> 
  dplyr::summarise(mean_values = mean(values),
                   mean_biomass = mean(biomass),
                   sd_biomass = sd(biomass)) |> 
  dplyr::mutate(var_type = dplyr::case_when(var %in% c("SST max", "SST min", "DHW max", "NPP mean", "SSS mean", "pH min") ~ "ENV",
                                            var %in% c("MPA", "GDP", "Gravity", "Fishing", "Neartt", "MED") ~ "HUM",
                                            var %in% c("Rock (%)", "Sand (%)", "Coral/Algae (%)", "Coral (RLS)", "Depth", "Reef extent") ~ "HAB"))

sp_n_realm <- partial_realm |> 
  dplyr::select(realm, species_name) |> 
  unique() |> 
  dplyr::group_by(realm) |> 
  dplyr::summarise(n = dplyr::n())

partial_realm_mean <- partial_realm_mean |> 
  dplyr::inner_join(sp_n_realm) |> 
  dplyr::mutate(realm = paste0(realm, " (n = ", n, ")"))

central_indo_pacific_mean <- partial_realm_mean |> 
  dplyr::filter(var %in% c("NPP mean", "Neartt", "SST min", "Reef extent", "SST max", "Gravity"),
                realm == "Central Indo-Pacific (n = 21)") |> 
  dplyr::ungroup() |> 
  dplyr::select(-id) |> 
  dplyr::mutate(var = forcats::fct_relevel(var, c("NPP mean", "Neartt", "SST min", "Reef extent", "SST max", "Gravity"))) |> 
  ggplot() +
  geom_ribbon(aes(x = mean_values, y = mean_biomass, ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), fill = "grey", alpha = 0.5) +
  geom_line(aes(x = mean_values, y = mean_biomass, color = var_type), size = 1) +
  scale_color_manual(values = c("ENV" = pal_contribution[2],
                               "HUM" = pal_contribution[1],
                               "HAB" = pal_contribution[13],
                               "BIOT" = pal_contribution[3])) +
  facet_wrap(~var, scales = "free", nrow = 1) +
  theme_bw() +
  labs(title = "Central Indo-Pacific (n = 21)", x = "", y = "log10(Biomass+1)") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.position = "none")

eastern_indo_pacific_mean <- partial_realm_mean |> 
  dplyr::filter(var %in% c("MED", "GDP", "SSS mean", "SST min", "DHW max", "SST max"),
                realm == "Eastern Indo-Pacific (n = 68)") |> 
  dplyr::ungroup() |> 
  dplyr::select(-id) |> 
  dplyr::mutate(var = forcats::fct_relevel(var, c("MED", "GDP", "SSS mean", "SST min", "DHW max", "SST max"))) |> 
  ggplot() +
  geom_ribbon(aes(x = mean_values, y = mean_biomass, ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), fill = "grey", alpha = 0.5) +
  geom_line(aes(x = mean_values, y = mean_biomass, color = var_type), size = 1) +
  scale_color_manual(values = c("ENV" = pal_contribution[2],
                                "HUM" = pal_contribution[1],
                                "HAB" = pal_contribution[13],
                                "BIOT" = pal_contribution[3])) +
  facet_wrap(~var, scales = "free", nrow = 1) +
  theme_bw() +
  labs(title = "Eastern Indo-Pacific (n = 68)", x = "", y = "log10(Biomass+1)") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.position = "none")

tropical_atlantic_mean <- partial_realm_mean |> 
  dplyr::filter(var %in% c("MED", "SST min", "Depth", "Rock (%)", "Coral (RLS)", "pH min"),
                realm == "Tropical Atlantic (n = 46)") |> 
  dplyr::ungroup() |> 
  dplyr::select(-id) |> 
  dplyr::mutate(var = forcats::fct_relevel(var, c("MED", "SST min", "Depth", "Rock (%)", "Coral (RLS)", "pH min"))) |> 
  ggplot() +
  geom_ribbon(aes(x = mean_values, y = mean_biomass, ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), fill = "grey", alpha = 0.5) +
  geom_line(aes(x = mean_values, y = mean_biomass, color = var_type), size = 1) +
  scale_color_manual(values = c("ENV" = pal_contribution[2],
                                "HUM" = pal_contribution[1],
                                "HAB" = pal_contribution[13],
                                "BIOT" = pal_contribution[3])) +
  facet_wrap(~var, scales = "free", nrow = 1) +
  theme_bw() +
  labs(title = "Tropical Atlantic (n = 46)", x = "", y = "log10(Biomass+1)", color = " ") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.text = element_text(size = 12))

tropical_eastern_pacific_mean <- partial_realm_mean |> 
  dplyr::filter(var %in% c("Depth", "Gravity", "MED", "SSS mean", "NPP mean", "Neartt"),
                realm == "Tropical Eastern Pacific (n = 37)") |> 
  dplyr::ungroup() |> 
  dplyr::select(-id) |> 
  dplyr::mutate(var = forcats::fct_relevel(var, c("Depth", "Gravity", "MED", "SSS mean", "NPP mean", "Neartt"))) |> 
  ggplot() +
  geom_ribbon(aes(x = mean_values, y = mean_biomass, ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), fill = "grey", alpha = 0.5) +
  geom_line(aes(x = mean_values, y = mean_biomass, color = var_type), size = 1) +
  scale_color_manual(values = c("ENV" = pal_contribution[2],
                                "HUM" = pal_contribution[1],
                                "HAB" = pal_contribution[13])) +
  facet_wrap(~var, scales = "free", nrow = 1) +
  theme_bw() +
  labs(title = "Tropical Eastern Pacific (n = 37)", x = "", y = "log10(Biomass+1)") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.position = "none")

temperate_northern_atlantic_mean <- partial_realm_mean |> 
  dplyr::filter(var %in% c("Depth", "MED", "SSS mean", "Coral/Algae (%)", "Rock (%)", "Reef extent"),
                realm == "Temperate Northern Atlantic (n = 8)") |> 
  dplyr::ungroup() |> 
  dplyr::select(-id) |> 
  dplyr::mutate(var = forcats::fct_relevel(var, c("Depth", "MED", "SSS mean", "Coral/Algae (%)", "Rock (%)", "Reef extent"))) |> 
  ggplot() +
  geom_ribbon(aes(x = mean_values, y = mean_biomass, ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), fill = "grey", alpha = 0.5) +
  geom_line(aes(x = mean_values, y = mean_biomass, color = var_type), size = 1) +
  scale_color_manual(values = c("ENV" = pal_contribution[2],
                                "HUM" = pal_contribution[1],
                                "HAB" = pal_contribution[13])) +
  facet_wrap(~var, scales = "free", nrow = 1) +
  theme_bw() +
  labs(title = "Temperate Northern Atlantic (n = 8)", x = "", y = "log10(Biomass+1)") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.position = "none")

realm_partial_mean_plot <- central_indo_pacific_mean / eastern_indo_pacific_mean / tropical_atlantic_mean / tropical_eastern_pacific_mean / temperate_northern_atlantic_mean

ggsave("figures/realm_partial_mean_plot.png", realm_partial_mean_plot, width = 11, height = 13)

