# script to produce figures that evaluate model performances

library(ggplot2)

# source functions ----

source("R/04_model_performance_functions.R")

# Set palette colors for performance figures

pal_best <- PNWColors::pnw_palette("Bay", 6, type = "continuous")
pal_perf <- PNWColors::pnw_palette("Bay", 6, type = "continuous")

# select best fitted model for each model type based on a concensus metrics ----

############### read data global ###############

load("outputs/performance_model.Rdata")

# estimate for each species the best model based on performance metrics  
best_models <- performance_bind |> 
  aggregate_metrics(., 
                    metrics = c("pearson", "spearman")) |>
  # find the best fitting model for each species within each fitted_model
  dplyr::group_by(species_name) |> 
  dplyr::do(best_model = .$model[which.max(.$discrimination)]) |> 
  tidyr::unnest(cols = c("best_model"))

save(best_models, file = "outputs/best_models.Rdata")

#### Best Model plot ####

best_models_pr <- best_models |>  
  dplyr::group_by(best_model) |> 
  dplyr::summarise(n = dplyr::n()) |> 
  dplyr::mutate(pr = (n*100)/sum(n))

best_models_pr$colors <- NA
best_models_pr[best_models_pr$best_model == "glm",4]$colors <- pal_best[1]
best_models_pr[best_models_pr$best_model == "gam",4]$colors <- pal_best[2]
best_models_pr[best_models_pr$best_model == "spamm",4]$colors <- pal_best[3]
best_models_pr[best_models_pr$best_model == "gbm",4]$colors <- pal_best[4]
best_models_pr[best_models_pr$best_model == "sprf",4]$colors <- pal_best[5]
best_models_pr[best_models_pr$best_model == "rf",4]$colors <- pal_best[6]

best_model <- best_models_pr |> 
  dplyr::mutate(best_model = forcats::fct_reorder(best_model, pr)) |> 
  ggplot(aes(x = best_model, y = pr, fill = best_model)) +
  geom_bar(width = 0.8, stat = 'identity') +
  scale_fill_manual(values = c("glm" = pal_best[1],
                               "gam" = pal_best[2],
                               "spamm" = pal_best[3],
                               "gbm" = pal_best[4],
                               "sprf" = pal_best[5],
                               "rf" = pal_best[6])) +
  labs(x = "Statistic methods", y = "Best model (%)", fill = "Method", title = "B") +
  theme(title = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 25),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

perf_corr <- performance_bind |> 
  dplyr::group_split(model)
perf_corr <- lapply(1:length(perf_corr), function(i) {
  
  model_i <- perf_corr[[i]]
  model_i <- model_i[colnames(model_i) %in% c("species_name", "pearson", "model")]
  colnames(model_i)[colnames(model_i) %in% "pearson"] <- unique(model_i$model)
  model_i <- model_i |> 
    dplyr::select(-model)
  
})
perf_corr <- purrr::reduce(perf_corr, dplyr::full_join)

library(ggplot2)

model_corr_plot <- GGally::ggcorr(perf_corr, label = TRUE)

ggsave(model_corr_plot, filename = "figures/model_corr_plot.png")

# produce histograms of model performance for best models ----

p_level <- unique(performance_bind$model)

perf_models_all <- performance_bind |> 
  dplyr::summarise_at(dplyr::vars(pearson, spearman), list(function(x) list(Q0.05 = round(quantile(x, 0.05, na.rm = T), 2),
                                                                           IQR0.25 = round(quantile(x, 0.25, na.rm = T), 2),
                                                                           median  = round(median(x, na.rm = T), 2),
                                                                           IQR0.75 = round(quantile(x, 0.75, na.rm = T), 2),
                                                                           Q0.95 = round(quantile(x, 0.95, na.rm = T), 2)))) |> 
  dplyr::mutate(summary_value = c('Q0.05','IQR0.25', 'median', 'IQR0.75', 'Q0.95')) |> 
  tidyr::unnest() |> 
  dplyr::group_by(summary_value) |> 
  dplyr::select(pearson, spearman) |> 
  t() |> 
  data.frame()

perf_models_all <- perf_models_all |> 
  dplyr::mutate(metric = rownames(perf_models_all)) |>
  dplyr::select(metric, X1:X5)

perf_models_details <- parallel::mclapply(1:length(p_level), function(i) {
  
  p_level <- p_level[i]
  
  perf_models <- performance_bind |> 
    dplyr::filter(model == p_level) |> 
    dplyr::summarise_at(dplyr::vars(pearson, spearman), list(function(x) list(Q0.05 = round(quantile(x, 0.05, na.rm = T), 2),
                                                                           IQR0.25 = round(quantile(x, 0.25, na.rm = T), 2),
                                                                           median  = round(median(x, na.rm = T), 2),
                                                                           IQR0.75 = round(quantile(x, 0.75, na.rm = T), 2),
                                                                           Q0.95 = round(quantile(x, 0.95, na.rm = T), 2)))) |> 
    dplyr::mutate(summary_value = c('Q0.05', 'IQR0.25', 'median', 'IQR0.75', 'Q0.95')) |> 
    tidyr::unnest() |> 
    dplyr::group_by(summary_value) |> 
    dplyr::select(pearson, spearman) |> 
    t() |> 
    data.frame() 
  
  perf_models <- perf_models |> 
    dplyr::mutate(metric = rownames(perf_models)) |> 
    dplyr::select(metric, X1:X5)
  
}, mc.cores = 1)
names(perf_models_details) <- p_level

# Select only the best model for each species

best_assessments_SCV <- dplyr::inner_join(performance_bind, best_models, by = "species_name")
best_assessments_SCV <- best_assessments_SCV[best_assessments_SCV$model == best_assessments_SCV$best_model,]

perf_models_all_best <- best_assessments_SCV |> 
  dplyr::summarise_at(dplyr::vars(pearson, spearman), list(function(x) list(Q0.05 = round(quantile(x, 0.05, na.rm = T), 2),
                                                                  IQR0.25 = round(quantile(x, 0.25, na.rm = T), 2),
                                                                  median  = round(median(x, na.rm = T), 2),
                                                                  IQR0.75 = round(quantile(x, 0.75, na.rm = T), 2),
                                                                  Q0.95 = round(quantile(x, 0.95, na.rm = T), 2)))) |> 
  dplyr::mutate(summary_value = c('Q0.05','IQR0.25', 'median', 'IQR0.75', 'Q0.95')) |> 
  tidyr::unnest() |> 
  dplyr::group_by(summary_value) |> 
  dplyr::select(pearson, spearman) |> 
  t() |> 
  data.frame() 

perf_models_all_best <- perf_models_all_best |> 
  dplyr::mutate(metric = rownames(perf_models_all_best)) |> 
  dplyr::select(metric, X1:X5)
max(best_assessments_SCV$pearson)
max(best_assessments_SCV$spearman)

# Manage data for performance plot

performance_all <- performance_bind |> 
  dplyr::select(-c(r2_sd, sd_pearson, sd_spearman)) |> 
  tidyr::pivot_longer(c(r2, spearman, pearson), names_to = "metric", values_to = "value")

performance_best <- best_assessments_SCV |> 
  dplyr::select(-c(r2_sd, sd_pearson, sd_spearman)) |> 
  tidyr::pivot_longer(c(r2, spearman, pearson), names_to = "metric", values_to = "value")

performance_all_best <- dplyr::full_join(performance_all, performance_best)
performance_all_best$cat <- NA

performance_all_best[which(performance_all_best$model == performance_all_best$best_model),6] <- "Best models"
performance_all_best[which(is.na(performance_all_best$cat) == TRUE),6] <- "All models"

plot_pearson <- performance_plot(performance_all_best,
                                 metrics_sel = "pearson",
                                 slope = 0,
                                 intercept = 1,
                                 color = pal_perf,
                                 ylim = c(-0.4, 1),
                                 legend.position = 'none',
                                 plot_title = "A.")

plot_spearman <- performance_plot(performance_all_best,
                                  metrics_sel = "spearman",
                                  slope = 0,
                                  intercept = 1,
                                  color = pal_perf,
                                  ylim = c(-0.4, 1),
                                  legend.position = c(0.5, 0.1),
                                  plot_title = "")

plot_perf <- patchwork::wrap_plots(plot_pearson, plot_spearman)

ggplot2::ggsave("figures/plot_perf.pdf", plot_perf, height = 6, width = 12)
ggplot2::ggsave("figures/plot_best.pdf", best_model, height = 6, width = 12)

################## Plot performance-traits relationship ################## 

load("data/new_derived_data/species_count.Rdata")
sp_count <- sp_count |> 
  dplyr::rename(occurence = count)
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

best_count <- best_assessments_SCV |> 
  dplyr::inner_join(sp_count) |> 
  dplyr::inner_join(sp_car[,colnames(sp_car) %in% c("species_name", "MaxLength", "IUCN_status")])

plot_spear_count <- ggplot(best_count, aes(x = log10(occurence), y = spearman)) +
  geom_point() +
  geom_smooth(method = "lm")

# ggsave(plot_spear_count, file = "figures/spearman_occurence.png")

summary(lm(formula = "spearman ~ occurence", data = best_count))

plot_pear_count <- ggplot(best_count, aes(x = log10(occurence), y = pearson)) +
  geom_point() +
  geom_smooth(method = "lm")

# ggsave(plot_pear_count, file = "figures/pearson_occurence.png")

summary(lm(formula = "pearson ~ occurence", data = best_count))

best_trait <- best_assessments_SCV |> 
  dplyr::inner_join(sp_car)

best_trait |> 
  dplyr::group_by(ML_cat) |> 
  dplyr::summarise(n = dplyr::n())

best_trait |> 
  dplyr::group_by(Water.column) |> 
  dplyr::summarise(n = dplyr::n())

best_trait |> 
  dplyr::group_by(Habitat) |> 
  dplyr::summarise(n = dplyr::n())

best_trait |> 
  dplyr::group_by(Trophic_guild_name) |> 
  dplyr::summarise(n = dplyr::n())

kruskal_ML <- kruskal_test_function_perf(best_trait,
                                         "ML_cat",
                                         "spearman")

if(kruskal_ML$statistics$p.chisq > 0.05){
  
  stop(print("No statistical differences among groups"))
  
}else{
  
  print("p.chisq < 0.05")
  
}

plot_MLclass_Spear <- ggplot(best_trait, aes(x = ML_cat, y = spearman), group = ML_cat) +
  geom_boxplot() +
  geom_text(data = kruskal_ML$groups, aes_string(x = "trait", y = "quant", label = "groups"), size = 5, vjust = -0.5, hjust = -0.55) +
  labs(x = "Size class")

plot_ML_Spear <- ggplot(best_count, aes(x = MaxLength, y = spearman), group = ML_cat) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Maximum length")

# ggsave(plot_MLclass_Spear, file = "figures/maximumlengthclass_spearman.png")
# ggsave(plot_ML_Spear, file = "figures/maximumlength_spearman.png")

kruskal_WC <- kruskal_test_function_perf(best_trait,
                                         "Water.column",
                                         "spearman")

if(kruskal_WC$statistics$p.chisq > 0.05){
  
  stop(print("No statistical differences among groups"))
  
}else{
  
  print("p.chisq < 0.05")
  
}

plot_WC_Spear <- ggplot(best_trait, aes(x = Water.column, y = spearman), group = Water.column) +
  geom_boxplot() +
  geom_text(data = kruskal_WC$groups, aes_string(x = "trait", y = "quant", label = "groups"), size = 5, vjust = -0.5, hjust = -0.55) +
  labs(x = "Water column")

# ggsave(plot_WC_Spear, file = "figures/watercolumn_spearman.png")

kruskal_HAB <- kruskal_test_function_perf(best_trait,
                                         "Habitat",
                                         "spearman")

if(kruskal_HAB$statistics$p.chisq > 0.05){
  
  stop(print("No statistical differences among groups"))
  
}else{
  
  print("p.chisq < 0.05")
  
}

plot_HAB_Spear <- ggplot(best_trait, aes(x = Habitat, y = spearman), group = Habitat) +
  geom_boxplot() +
  geom_text(data = kruskal_HAB$groups, aes_string(x = "trait", y = "quant", label = "groups"), size = 5, vjust = -0.5, hjust = -0.55)

# ggsave(plot_HAB_Spear, file = "figures/habitat_spearman.png")

kruskal_TR <- kruskal_test_function_perf(best_trait,
                                         "Trophic_guild_name",
                                         "spearman")

if(kruskal_TR$statistics$p.chisq > 0.05){
  
  stop(print("No statistical differences among groups"))
  
}else{
  
  print("p.chisq < 0.05")
  
}

plot_TR_Spear <- ggplot(best_trait, aes(x = Trophic_guild_name, y = spearman), group = Trophic_guild_name) +
  geom_boxplot() +
  geom_text(data = kruskal_TR$groups, aes_string(x = "trait", y = "quant", label = "groups"), size = 5, vjust = -0.5, hjust = -0.55)

# ggsave(plot_TR_Spear, file = "figures/trophicgroup_spearman.png", width = 10)


occ_ml_sp <- ggplot(best_count, aes(x = log10(occurence), y = MaxLength, color = spearman)) + 
  geom_point() +
  geom_point(data = best_count[best_count$IUCN_status %in% "Thr",], colour = "red", pch = 21, size = 4, stroke = 1) +
  scale_colour_viridis_c()

best_count_log <- best_count |> 
  dplyr::mutate(occurence = log10(occurence))

glm_occ_ml_sp <- glm(formula = "spearman ~ occurence * MaxLength", data = best_count_log)
summary(glm_occ_ml_sp)
with(summary(glm_occ_ml_sp), 1 - deviance/null.deviance)

vis_plot <- visreg::visreg2d(glm_occ_ml_sp, "occurence", "MaxLength", scale = "response", plot.type = "gg") + 
  geom_contour(aes(z = z), color = "black") +
  labs(y = "Maximum length (cm)", x = "log10(Occurrence)", fill = "spearman") +
  scale_fill_viridis_c(option = "viridis") +
  theme(title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(vis_plot, file = "figures/visreg_occ_ml_sp.pdf", width = 8)
