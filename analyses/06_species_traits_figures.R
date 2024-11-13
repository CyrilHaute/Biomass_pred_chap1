# source functions ----

library(patchwork)

source("R/06_species_traits_figures_function.R")
source("R/05_contributions_figures_functions.R")

pal_sp_trait <- PNWColors::pnw_palette("Bay", 3, type = "discrete")

load("outputs/best_models.Rdata")

#### Covariates contribution plot ####

# Merge all files together by species

glm <- list.files("outputs/glm_contribution", full.names = TRUE)
glm <- lapply(1:length(glm), function(i) {
  
  load(glm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
glm <- do.call(rbind, glm)

gam <- list.files("outputs/gam_contribution", full.names = TRUE)
gam <- lapply(1:length(gam), function(i) {
  
  load(gam[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gam <- do.call(rbind, gam)

spamm <- list.files("outputs/spamm_contribution", full.names = TRUE)
spamm <- lapply(1:length(spamm), function(i) {
  
  load(spamm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
spamm <- do.call(rbind, spamm)

rf <- list.files("outputs/rf_contribution", full.names = TRUE)
rf <- lapply(1:length(rf), function(i) {
  
  load(rf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
rf <- do.call(rbind, rf)

sprf <- list.files("outputs/sprf_contribution", full.names = TRUE)
sprf <- lapply(1:length(sprf), function(i) {
  
  load(sprf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
sprf <- do.call(rbind, sprf)
colnames(sprf)[colnames(sprf) %in% c("global_dropout_loss", "global_sd_dropout_loss")] <- c("Dropout_loss", "sd_dropout_loss")

gbm <- list.files("outputs/brt_contribution", full.names = TRUE)
gbm <- lapply(1:length(gbm), function(i) {
  
  load(gbm[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
gbm <- do.call(rbind, gbm)

bind_files <- list(glm, gam, spamm, rf, sprf, gbm)
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
                                        x_text_angle = 65,
                                        legend.position = c(0.85, 0.8))

plot_species_traits <- plot_max.length / plot_water.column / plot_habitat / plot_trophic

ggsave("figures/plot_species_traits.pdf", plot_species_traits, height = 25, width = 19)



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
                                        labs_title = "",
                                        x_text_angle = 65,
                                        legend.position = c(0.85, 0.8))

ggsave("figures/plot_family.pdf", plot_familly, height = 10, width = 17)



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
  scale_fill_manual(values = c("ENV" = pal_sp_trait[1], 
                               "HUM" = pal_sp_trait[3], 
                               "HAB" = pal_sp_trait[2])) +
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
