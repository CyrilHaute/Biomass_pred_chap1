# source functions ----

library(patchwork)

source("R/06_species_traits_figures_function.R")

pal_sp_trait <- PNWColors::pnw_palette("Bay", 3, type = "discrete")

load("outputs/best_models.Rdata")

#### Covariates contribution plot ####

# Merge all files together by species

sprf <- list.files("outputs/biomass_contribution/sprf", full.names = T)

bind_files_sprf <- lapply(1:length(sprf), function(i) {
  
  load(sprf[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
bind_files_sprf <- do.call(rbind, bind_files_sprf)
bind_files_sprf$contributions_and_sd <- lapply(1:nrow(bind_files_sprf), function(i) { bind_files_sprf$contributions_and_sd[[i]] |> 
    dplyr::rename(Dropout_loss = global_dropout_loss,
                  sd_dropout_loss = global_sd_dropout_loss)
})
bind_files_sprf <- list(bind_files_sprf)

list_files_path <- list.files("outputs/biomass_contribution", full.names = T)
list_files_path <- list_files_path[which(grepl(pattern = "sprf", list_files_path) == FALSE)]
bind_files <- lapply(1:length(list_files_path), function(i) {
  
  load(list_files_path[i])
  assign(paste0("model_", i), extracted_contributions)
  
})
bind_files <- c(bind_files, bind_files_sprf)
bind_files <- do.call(rbind, bind_files)

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
                                           trait = "ML_cat",
                                           color = pal_sp_trait,
                                           labs_title = "A. Maximum length (cm)",
                                           aes_string_x = c(3.1, 4.85, 4.1, 2.9, 5.12, 2.12, 3.9, 1.9, 0.85, 1.12, 5.37, 4.37, 3.35, 2.3, 1.3),
                                           x_text_angle = NULL,
                                           legend.position = "none")

plot_water.column <- species_traits_function(plot_data = bind_files,
                                             trait = "Water.column",
                                             color = pal_sp_trait,
                                             labs_title = "B. Water column position ",
                                             aes_string_x = c(3.05, 2.8, 2.05, 1.8, 1.05, 3.3, 0.8, 2.3, 1.3),
                                             x_text_angle = NULL,
                                             legend.position = "none")

plot_habitat <- species_traits_function(plot_data = bind_files,
                                        trait = "Habitat",
                                        color = pal_sp_trait,
                                        labs_title = "C. Habitat",
                                        aes_string_x = c(4.05, 3.85, 2.1, 3.1, 1.825, 0.85, 1.1, 2.85, 1.3, 4.3, 3.3, 2.3),
                                        x_text_angle = NULL,
                                        legend.position = "none")

plot_trophic <- species_traits_function(plot_data = bind_files,
                                        trait = "Trophic_guild_name",
                                        color = pal_sp_trait,
                                        labs_title = "D. Trophic classes",
                                        aes_string_x = c(5.82, 2.82, 6.1, 3.1, 2.1, 4.05, 8.15, 3.83, 1.85, 6.85, 1.1, 7.87, 4.87, 5.13, 0.88, 7.15, 3.4, 7.45, 4.45, 6.45, 2.35, 8.4, 5.4, 1.3),
                                        x_text_angle = 65,
                                        legend.position = c(0.85, 0.8))

plot_species_traits <- plot_max.length / plot_water.column / plot_habitat / plot_trophic

ggsave("figures/plot_species_traits.pdf", plot_species_traits, height = 25, width = 19)
ggsave("figures/plot_species_traits_presentation.png", plot_species_traits, height = 25, width = 19)

plot_species_traits <- (plot_max.length / plot_water.column)
ggsave("figures/plot_species_traits_presentation.png", plot_species_traits, height = 13, width = 16)
plot_species_traits2 <- (plot_habitat / plot_trophic)
ggsave("figures/plot_species_traits_presentation2.png", plot_species_traits2, height = 13, width = 16)



bind_files_acp <- bind_files |> 
  dplyr::filter(fitted_model == "rf")

bind_files_acp <- pbmcapply::pbmclapply(1:nrow(bind_files_acp), function(i) {
  
  sp_model_i <- bind_files_acp[i,]
  
  pivot_cont <- sp_model_i$contributions_and_sd[[1]] |> 
    dplyr::select(-sd_dropout_loss) |> 
    tidyr::pivot_wider(names_from = "variable",
                     values_from = "Dropout_loss") |> 
    dplyr::mutate(species_name = sp_model_i$species_name,
                  fitted_model = sp_model_i$fitted_model)
  
}, mc.cores = parallel::detectCores() - 1)

bind_files_acp_bind <- do.call(rbind, bind_files_acp)

bind_files_acp_bind_traits <- dplyr::inner_join(bind_files_acp_bind, sp_car) |> 
  tidyr::drop_na()





ENV <- lapply(1:nrow(bind_files_acp_bind_traits), function(i) { bind_files_acp_bind_traits[i,][,which(colnames(bind_files_acp_bind_traits) %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"))]})
ENV <- do.call(rbind, ENV)
ENV <- dplyr::tibble(species_name = bind_files_acp_bind_traits$species_name,
                     value = rowMeans(ENV),
                     var = rep("ENV",nrow(ENV)))

SOC <- lapply(1:nrow(bind_files_acp_bind_traits), function(i) { bind_files_acp_bind_traits[i,][,which(colnames(bind_files_acp_bind_traits) %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"))]})
SOC <- do.call(rbind, SOC)
SOC <- dplyr::tibble(species_name = bind_files_acp_bind_traits$species_name,
                     value = rowMeans(SOC),
                     var = rep("HUM", nrow(SOC)))

HAB <- lapply(1:nrow(bind_files_acp_bind_traits), function(i) { bind_files_acp_bind_traits[i,][,which(colnames(bind_files_acp_bind_traits) %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"))]})
HAB <- do.call(rbind, HAB)
HAB <- dplyr::tibble(species_name = bind_files_acp_bind_traits$species_name,
                     value = rowMeans(HAB),
                     var = rep("HAB", nrow(HAB)))

cont <- ENV |>  
  dplyr::full_join(HAB) |>
  dplyr::full_join(SOC)

best <- dplyr::tibble(species_name = unique(cont$species_name),
                      ENV = cont[cont$var == "ENV",2],
                      SOC = cont[cont$var == "HUM",2],
                      HAB = cont[cont$var == "HAB",2])
best <- dplyr::inner_join(best, best_models, by = "species_name")
best <- as.matrix(best)

best <- lapply(1:nrow(best), function(i) {
  test <- best[i,]
  testt <- which.max(test[2:4])
  testtt <- dplyr::tibble(species_name = test[1],
                          varmax = testt)
})

best <- do.call(rbind, best)
best$varmax <- as.character(best$varmax)
best[best$varmax == "1" ,2] <- "ENV"
best[best$varmax == "2" ,2] <- "HUM"
best[best$varmax == "3" ,2] <- "HAB"


bind_files_acp_bind_traits <- bind_files_acp_bind_traits |> 
  dplyr::inner_join(best)









pca <- FactoMineR::FAMD(bind_files_acp_bind_traits[c(27:42)], graph = FALSE)

# bind_files_acp_bind_traits$ <- 

factoextra::fviz_famd(pca,
                      geom = c("point"),
                      repel = TRUE,
                      fill.ind = as.factor(bind_files_acp_bind_traits$varmax),
                      col.ind = as.factor(bind_files_acp_bind_traits$varmax),
                      palette = PNWColors::pnw_palette("Bay", 3, type = "discrete")
                      # pointshape = 21,
                      # pointsize = 3
                      )



pca <- FactoMineR::PCA(bind_files_acp_bind_traits[c(1:24)], graph = FALSE)

# bind_files_acp_bind_traits$ <- 

factoextra::fviz_pca_biplot(pca,
                            geom = "point",
                            repel = TRUE,
                            fill.ind = as.factor(bind_files_acp_bind_traits$Trophic_guild_name),
                            col.ind = as.factor(bind_files_acp_bind_traits$Trophic_guild_name),
                            pointshape = 21,
                            pointsize = 3)





cont <- lapply(1:length(fitted_model), function(i) {
  
  plot_level <- fitted_model[i]
  only_model <- bind_files |> 
    dplyr::filter(fitted_model == plot_level)
  only_model <- only_model |>
    dplyr::filter(species_name %in% best_models[best_models$best_model == plot_level,1]$species_name)
  
  ENV <- lapply(1:nrow(only_model), function(i) { only_model$contributions_and_sd[[i]][only_model$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$Dropout_loss})
  ENV_sd <- lapply(1:nrow(only_model), function(i) { only_model$contributions_and_sd[[i]][only_model$contributions_and_sd[[i]]$variable %in% c("max_1year_analysed_sst", "max_5year_degree_heating_week", "mean_1year_chl", "mean_1year_so_mean", "mean_7days_analysed_sst", "mean_7days_chl", "min_1year_analysed_sst", "min_5year_ph"),]$sd_dropout_loss})
  ENV <- do.call(rbind, ENV)
  ENV_sd <- do.call(rbind, ENV_sd)
  ENV <- dplyr::tibble(species_name = only_model$species_name,
                       value = matrixStats::rowMedians(ENV),
                       sd = matrixStats::rowMedians(ENV_sd),
                       var = rep("ENV",nrow(ENV)),
                       plot_level = rep(plot_level, nrow(ENV)))
  
  SOC <- lapply(1:nrow(only_model), function(i) { only_model$contributions_and_sd[[i]][only_model$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$Dropout_loss})
  SOC_sd <- lapply(1:nrow(only_model), function(i) { only_model$contributions_and_sd[[i]][only_model$contributions_and_sd[[i]]$variable %in% c("effectiveness", "gdp", "gravtot2", "hdi", "n_fishing_vessels", "natural_ressource_rent", "neartt", "ngo"),]$sd_dropout_loss})
  SOC <- do.call(rbind, SOC)
  SOC_sd <- do.call(rbind, SOC_sd)
  SOC <- dplyr::tibble(species_name = only_model$species_name,
                       value = matrixStats::rowMedians(SOC),
                       sd = matrixStats::rowMedians(SOC_sd),
                       var = rep("HUM",nrow(SOC)),
                       plot_level = rep(plot_level, nrow(SOC)))
  
  HAB <- lapply(1:nrow(only_model), function(i) { only_model$contributions_and_sd[[i]][only_model$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$Dropout_loss})
  HAB_sd <- lapply(1:nrow(only_model), function(i) { only_model$contributions_and_sd[[i]][only_model$contributions_and_sd[[i]]$variable %in% c("Rock_500m", "Rubble_500m", "Sand_500m", "coral", "coral_algae_500m", "coralline_algae", "depth", "reef_extent"),]$sd_dropout_loss})
  HAB <- do.call(rbind, HAB)
  HAB_sd <- do.call(rbind, HAB_sd)
  HAB <- dplyr::tibble(species_name = only_model$species_name,
                       value = matrixStats::rowMedians(HAB),
                       sd = matrixStats::rowMedians(HAB_sd),
                       var = rep("HAB",nrow(HAB)),
                       plot_level = rep(plot_level,nrow(HAB)))
  
  cont <- ENV |> 
    dplyr::full_join(HAB) |> 
    dplyr::full_join(SOC)
  cont$plot_level <- rep(plot_level, nrow(cont))
  cont
  
})

cont <- do.call(rbind, cont)

load("data/new_raw_data/RLS_actinopterygii_data.Rdata")
phylo <- RLS_actinopterygii_data |> 
  dplyr::select(species_name, order, family)
phylo <- unique(phylo)

cont <- cont |> 
  dplyr::inner_join(phylo, by = "species_name")

n_trait <- cont |>
  dplyr::group_by(order) |> 
  dplyr::summarise(n = dplyr::n() / 3)

cont <- cont |> 
  dplyr::inner_join(n_trait, by = "order") |> 
  dplyr::filter(n >= 10)

cont[colnames(cont) %in% "order"] <- sapply(1:nrow(cont[colnames(cont) %in% "order"]), function(i) {
  
  row_i <- cont[colnames(cont) %in% "order"][i,]
  
  which_row <- which(grepl(row_i, unlist(n_trait[colnames(n_trait) %in% "order"])) == TRUE)
  
  paste0(row_i, " (n = ", n_trait$n[which_row], ")")
  
})

library(ggplot2)

plot_trait <- cont |> 
  dplyr::mutate(var = forcats::fct_relevel(var, "ENV", "HAB", "HUM")) |> 
  ggplot(aes_string(x = "order", y = "value", fill = "var")) +
  geom_boxplot(aes(fill = factor(var)), width=0.6, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("ENV" = pal_sp_trait[1],
                               "HUM" = pal_sp_trait[3],
                               "HAB" = pal_sp_trait[2])) +
  theme_bw() +
  # geom_text(data = kruskal_test_trait$groups, aes_string(x = aes_string_x, y = "quant", label = "groups"), vjust=-0.48, size = 6) +
  coord_cartesian(ylim = c(0,0.25)) +
  labs(y = "Relative importance (RMSE)", x = "", fill = "") +
  theme(
        legend.direction = "vertical",
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA),
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("figures/plot_family.pdf", plot_trait, height = 10, width = 20)
