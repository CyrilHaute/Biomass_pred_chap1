# script to produce figures that evaluate model performances

library(ggplot2)
library(GGally)

# source functions ----

source("R/04_model_performance_functions.R")
source("R/04_model_prediction_functions.R")
source("R/04_model_correlation_function.R")

# Set palette colors for performance figures

pal_best <- PNWColors::pnw_palette("Bay", 6, type = "continuous")
pal_perf <- PNWColors::pnw_palette("Bay", 6, type = "continuous")

# select best fitted model for each model type based on a concensus metrics ----

############### read data global ###############

load("outputs/performance_model.Rdata")
load("data/new_derived_data/species_count.Rdata")

# estimate for each species the best model based on performance metrics  
best_models <- performance_bind |> 
  aggregate_metrics(metrics = c("pearson", "spearman")) |>
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

best_models_pr[best_models_pr$best_model == "gam",]$best_model <- "GAM"
best_models_pr[best_models_pr$best_model == "gbm",]$best_model <- "GBM"
best_models_pr[best_models_pr$best_model == "glm",]$best_model <- "GLM"
best_models_pr[best_models_pr$best_model == "rf",]$best_model <- "RF"
best_models_pr[best_models_pr$best_model == "spamm",]$best_model <- "SPAMM"
best_models_pr[best_models_pr$best_model == "sprf",]$best_model <- "SPRF"

best_model <- best_models_pr |> 
  dplyr::mutate(best_model = forcats::fct_reorder(best_model, pr)) |> 
  ggplot(aes(x = best_model, y = pr, fill = best_model)) +
  geom_bar(width = 0.8, stat = "identity") +
  scale_fill_manual(values = c("GLM" = pal_best[1],
                               "GAM" = pal_best[2],
                               "SPAMM" = pal_best[3],
                               "GBM" = pal_best[4],
                               "SPRF" = pal_best[5],
                               "RF" = pal_best[6])) +
  labs(x = "Statistical model", y = "Best model (%)", fill = "Method", title = "A.") +
  theme_minimal() + 
  theme(title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Pearon and Spearman plot
# Select only the best model for each species

best_assessments_SCV <- dplyr::inner_join(performance_bind, best_models, by = "species_name")
best_assessments_SCV <- best_assessments_SCV[best_assessments_SCV$model == best_assessments_SCV$best_model,]

perf_model_gouped_best <- best_assessments_SCV |> 
  tidyr::drop_na() |> 
  dplyr::group_by(model) |>
  dplyr::summarise(medianpearson = median(pearson, na.rm = TRUE),
                   q95_pearson = quantile(pearson, probs = 0.95, na.rm = TRUE),
                   q05_pearson = quantile(pearson, probs = 0.05, na.rm = TRUE),
                   maxpearson = max(pearson, na.rm = TRUE),
                   minpearson = min(pearson, na.rm = TRUE),
                   sdpearson = sd(pearson, na.rm = TRUE),
                   medianspearman = median(spearman, na.rm = TRUE),
                   q95_spearman = quantile(spearman, probs = 0.95, na.rm = TRUE),
                   q05_spearman = quantile(spearman, probs = 0.05, na.rm = TRUE),
                   maxspearman = max(spearman, na.rm = TRUE),
                   minspearman = min(spearman, na.rm = TRUE),
                   sdspearman = sd(spearman, na.rm = TRUE))
perf_model_best <- best_assessments_SCV |> 
  tidyr::drop_na() |>
  dplyr::summarise(medianpearson = median(pearson, na.rm = TRUE),
                   q95_pearson = quantile(pearson, probs = 0.95, na.rm = TRUE),
                   q05_pearson = quantile(pearson, probs = 0.05, na.rm = TRUE),
                   maxpearson = max(pearson, na.rm = TRUE),
                   minpearson = min(pearson, na.rm = TRUE),
                   sdpearson = sd(pearson, na.rm = TRUE),
                   medianspearman = median(spearman, na.rm = TRUE),
                   q95_spearman = quantile(spearman, probs = 0.95, na.rm = TRUE),
                   q05_spearman = quantile(spearman, probs = 0.05, na.rm = TRUE),
                   maxspearman = max(spearman, na.rm = TRUE),
                   minspearman = min(spearman, na.rm = TRUE),
                   sdspearman = sd(spearman, na.rm = TRUE))
max(best_assessments_SCV$pearson)
max(best_assessments_SCV$spearman)

# Manage data for performance plot

performance_all <- performance_bind |> 
  tidyr::drop_na() |> 
  dplyr::select(-c(r2_sd, sd_pearson, sd_spearman)) |> 
  tidyr::pivot_longer(c(r2, spearman, pearson), names_to = "metric", values_to = "value")

performance_best <- best_assessments_SCV |> 
  dplyr::select(-c(r2_sd, sd_pearson, sd_spearman)) |> 
  tidyr::pivot_longer(c(r2, spearman, pearson), names_to = "metric", values_to = "value")

performance_all_best <- dplyr::full_join(performance_all, performance_best)
performance_all_best$cat <- NA

performance_all_best[which(performance_all_best$model == performance_all_best$best_model),6] <- "Best models"
performance_all_best[which(is.na(performance_all_best$cat) == TRUE),6] <- "All models"

performance_all_best[performance_all_best$model == "gam",]$model <- "GAM"
performance_all_best[performance_all_best$model == "gbm",]$model <- "GBM"
performance_all_best[performance_all_best$model == "glm",]$model <- "GLM"
performance_all_best[performance_all_best$model == "rf",]$model <- "RF"
performance_all_best[performance_all_best$model == "spamm",]$model <- "SPAMM"
performance_all_best[performance_all_best$model == "sprf",]$model <- "SPRF"

performance_all_best[performance_all_best$metric == "pearson",]$metric <- "Pearson"
performance_all_best[performance_all_best$metric == "spearman",]$metric <- "Spearman"

performance_all_best <- performance_all_best |> 
  dplyr::mutate(model2 = paste0(model, "_", metric))

performance_all_best <- performance_all_best |> 
  dplyr::filter(metric != c("r2"))

performance_all_best$model2 <- factor(performance_all_best$model2, levels = c("GLM_Pearson", "GAM_Pearson", "SPAMM_Pearson", "GBM_Pearson", "SPRF_Pearson", "RF_Pearson", "",
                                                                              "GLM_Spearman", "GAM_Spearman", "SPAMM_Spearman", "GBM_Spearman", "SPRF_Spearman", "RF_Spearman"))

mean_pearson <- mean(performance_all_best[performance_all_best$metric == "Pearson" & performance_all_best$cat == "All models",]$value)
mean_spearman <- mean(performance_all_best[performance_all_best$metric == "Spearman" & performance_all_best$cat == "All models",]$value)

plot_perf <- performance_all_best |> 
  dplyr::mutate(model2 = forcats::fct_relevel(model2, "GLM_Pearson", "GAM_Pearson", "SPAMM_Pearson", "GBM_Pearson", "SPRF_Pearson", "RF_Pearson", "",
                                              "GLM_Spearman", "GAM_Spearman", "SPAMM_Spearman", "GBM_Spearman", "SPRF_Spearman", "RF_Spearman")) |>
  ggplot() +
  geom_boxplot(aes(x = model2, y = value, fill = cat), 
               alpha = 0.8, 
               outlier.shape = NA, 
               size = 0.4, 
               position = position_dodge(width = 0.9)) +
  geom_segment(
    aes(x = 0, xend = 6.5, y = mean_pearson, yend = mean_pearson),
    color = "red",
    size = 2,
    linetype = "dashed"
  ) +
  geom_segment(
    aes(x = 7.5, xend = 13.5, y = mean_spearman, yend = mean_spearman),
    color = "red",
    size = 2,
    linetype = "dashed"
  ) +
  scale_x_discrete(drop = FALSE,
                   labels = c("GLM_Pearson" = "GLM",
                              "GAM_Pearson" = "GAM",
                              "SPAMM_Pearson" = "SPAMM",
                              "GBM_Pearson" = "GBM", 
                              "SPRF_Pearson" = "SPRF",
                              "RF_Pearson" = "RF",
                              "GLM_Spearman" = "GLM",
                              "GAM_Spearman" = "GAM",
                              "SPAMM_Spearman" = "SPAMM",
                              "GBM_Spearman" = "GBM",
                              "SPRF_Spearman" = "SPRF",
                              "RF_Spearman" = "RF")) +
  annotate("text", x = 3.5, y = 0.95, label = "Pearson", size = 5) +
  annotate("text", x = 10.5, y = 0.95, label = "Spearman", size = 5) +
  scale_fill_manual(values = c("Best models" = pal_perf[6],
                               "All models" = pal_perf[2])) +
  coord_cartesian(ylim = c(-0.3, 1)) + 
  labs(y = "Performance values", x = "", title = "C.", fill = "") +
  theme_minimal() + 
  theme(
    legend.position = c(0.8, 0.1),
    legend.direction = "horizontal",
    title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    panel.background = element_rect(fill = "white", colour = "grey50",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())


# Load data to plot model predictions

glm <- list.files("outputs/glm_prediction", full.names = TRUE)
rf <- list.files("outputs/rf_prediction", full.names = TRUE)
brt <- list.files("outputs/brt_prediction", full.names = TRUE)
gam <- list.files("outputs/gam_prediction", full.names = TRUE)
spamm <- list.files("outputs/spamm_prediction", full.names = TRUE)
sprf <- list.files("outputs/sprf_prediction", full.names = TRUE)

file_model <- c(glm, rf, brt, gam, spamm, sprf)

read_sp_eco <- lapply(1:length(file_model), function(i) {
  
  load(file_model[i])
  assign(as.character(i), cv_j_bind)
  
})

read_sp_eco <- pbmcapply::pbmclapply(1:length(read_sp_eco), function(i) {
  
  sp <- read_sp_eco[[i]]
  sp$survey_id <- as.numeric(sp$survey_id)
  to_remove <- which(sapply(sp$survey_id, is.na))
  if(length(to_remove) == 0){
    
    sp <- sp
    
  }else{
    
    sp <- sp[-to_remove,]
    
  }
  
  return(sp)
  
}, mc.cores = parallel::detectCores() - 1)

read_sp_eco <- do.call(rbind, read_sp_eco)

load("outputs/best_models.Rdata")

read_sp_eco_best <- read_sp_eco |> 
  dplyr::inner_join(best_models) |> 
  dplyr::filter(best_model == model)

# create plots 
observed_predicted_best_plot <- observed_predicted_best_plot(input_data = read_sp_eco_best, 
                                                             nbins = 25)

# create plots 
observed_predicted_all_plot <- observed_predicted_plot(input_data = read_sp_eco, 
                                                       nbins = 25, 
                                                       levels = c("glm", "gam", "spamm", "rf", "gbm", "sprf"))

ggsave("figures/all_predictions_pres.png", observed_predicted_all_plot, width = 9, height = 11)


# Plot Figure 1

plot_perf_best <- (best_model + observed_predicted_best_plot) / (plot_perf)

ggplot2::ggsave("figures/plot_perf_best.png", plot_perf_best, height = 13, width = 12)


################ Performance correlation figure ################

performance_bind[performance_bind$model == "gam",]$model <- "GAM"
performance_bind[performance_bind$model == "gbm",]$model <- "GBM"
performance_bind[performance_bind$model == "glm",]$model <- "GLM"
performance_bind[performance_bind$model == "rf",]$model <- "RF"
performance_bind[performance_bind$model == "spamm",]$model <- "SPAMM"
performance_bind[performance_bind$model == "sprf",]$model <- "SPRF"

perf_corr <- performance_bind |> 
  tidyr::drop_na() |> 
  dplyr::group_split(model)
perf_corr <- lapply(1:length(perf_corr), function(i) {
  
  model_i <- perf_corr[[i]]
  model_i <- model_i[colnames(model_i) %in% c("species_name", "pearson", "model")]
  colnames(model_i)[colnames(model_i) %in% "pearson"] <- unique(model_i$model)
  model_i <- model_i |> 
    dplyr::select(-model)
  
})
perf_corr <- purrr::reduce(perf_corr, dplyr::full_join) |> 
  tidyr::drop_na()
perf_corr <- perf_corr[, c("species_name", "GLM", "GAM", "SPAMM", "GBM", "RF", "SPRF")]

legend_plot <- ggplot(data.frame(x = 0, y = 0, cor = c(0.35, 0.85)), aes(fill = cor)) +
  geom_tile(aes(x = 1, y = cor)) +
  scale_fill_gradient2(
    low = "#3B9AB2", mid = "#EEEEEE", high = "#F21A00",
    midpoint = 0.6, limits = c(0.35, 0.85),
    name = "Pearson"
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.key.height = unit(1, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

# Extraire la légende
legend <- grab_legend(legend_plot)

corr_plot <- ggpairs(perf_corr[,-1], 
                     lower = list(continuous = custom_ggally_scat),
                     upper = list(continuous = custom_ggally_cor),
                     diag = list(continuous = custom_diag_text),
                     legend = legend)  +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(-0.2, "cm")
  )

ggsave(corr_plot, filename = "figures/model_corr_plot.png", width = 10, height = 7)
