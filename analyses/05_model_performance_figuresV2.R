# script to produce figures that evaluate model performances

library(ggplot2)
library(GGally)

# source functions ----

source("R/04_model_performance_functions.R")
source("R/10_model_correlation_function.R")

# Set palette colors for performance figures

pal_best <- PNWColors::pnw_palette("Bay", 6, type = "continuous")
pal_perf <- PNWColors::pnw_palette("Bay", 6, type = "continuous")

# select best fitted model for each model type based on a concensus metrics ----

############### read data global ###############

load("outputs/performance_model.Rdata")
load("data/new_derived_data/species_count.Rdata")

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

best_models_pr[best_models_pr$best_model == "gam",]$best_model <- "GAM"
best_models_pr[best_models_pr$best_model == "gbm",]$best_model <- "GBM"
best_models_pr[best_models_pr$best_model == "glm",]$best_model <- "GLM"
best_models_pr[best_models_pr$best_model == "rf",]$best_model <- "RF"
best_models_pr[best_models_pr$best_model == "spamm",]$best_model <- "SPAMM"
best_models_pr[best_models_pr$best_model == "sprf",]$best_model <- "SPRF"

best_model <- best_models_pr |> 
  dplyr::mutate(best_model = forcats::fct_reorder(best_model, pr)) |> 
  ggplot(aes(x = best_model, y = pr, fill = best_model)) +
  geom_bar(width = 0.8, stat = 'identity') +
  scale_fill_manual(values = c("GLM" = pal_best[1],
                               "GAM" = pal_best[2],
                               "SPAMM" = pal_best[3],
                               "GBM" = pal_best[4],
                               "SPRF" = pal_best[5],
                               "RF" = pal_best[6])) +
  labs(x = "Statistical model", y = "Best model (%)", fill = "Method", title = "A.") +
  theme_minimal() + 
  theme(title = element_text(size = 15),
        axis.text=element_text(size = 15),
        axis.title=element_text(size = 20),
        legend.text=element_text(size = 10),
        legend.title=element_text(size = 15),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

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

plot_perf_pearson <- performance_all_best |> 
  dplyr::filter(metric %in% c("Pearson")) |> 
  dplyr::mutate(model = forcats::fct_relevel(model, "GLM", "GAM", "SPAMM", "GBM", "SPRF", "RF")) |> 
  ggplot() +
  geom_boxplot(aes(x = model, y = value, fill = cat), 
               alpha = 0.8, 
               outlier.shape = NA, 
               size = 0.6, 
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Best models" = pal_perf[6],
                               "All models" = pal_perf[2])) +
  coord_cartesian(ylim = c(-0.3, 1)) + 
  labs(y = "Pearson", x = "", title = "C.", fill = "", linetype = "") +
  theme_minimal() + 
  theme(
    legend.position = "none",
    title = element_text(size = 15),
    axis.text=element_text(size = 15),
    axis.title=element_text(size = 20),
    legend.text=element_text(size = 10),
    legend.title=element_text(size = 15),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    panel.background = element_rect(fill = "white", colour = "grey50",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

plot_perf_spearman <- performance_all_best |> 
  dplyr::filter(metric %in% c("Spearman")) |> 
  dplyr::mutate(model = forcats::fct_relevel(model, "GLM", "GAM", "SPAMM", "GBM", "SPRF", "RF")) |> 
  ggplot() +
  geom_boxplot(aes(x = model, y = value, fill = cat), 
               alpha = 0.8, 
               outlier.shape = NA, 
               size = 0.6, 
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Best models" = pal_perf[6],
                               "All models" = pal_perf[2])) +
  coord_cartesian(ylim = c(-0.3, 1)) + 
  labs(y = "Spearman", x = "", title = "", fill = "", linetype = "") +
  theme_minimal() + 
  theme(
    legend.position = c(0.7, 0.1),
    legend.direction = "horizontal",
    title = element_text(size = 15),
    axis.text=element_text(size = 15),
    axis.title=element_text(size = 20),
    legend.text=element_text(size = 10),
    legend.title=element_text(size = 15),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    panel.background = element_rect(fill = "white", colour = "grey50",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

library(patchwork)

plot_perf_best <- (best_model + observed_predicted_best_plot) / (plot_perf_pearson + plot_perf_spearman)

ggplot2::ggsave("figures/test.png", plot_perf_best, height = 13, width = 12)



















plot_perf_all <- performance_all_best |> 
  dplyr::filter(cat %in% c("All models"),
                metric %in% c("Pearson", "Spearman")) |> 
  dplyr::mutate(model = forcats::fct_relevel(model, "GLM", "GAM", "SPAMM", "GBM", "SPRF", "RF")) |> 
  ggplot() +
  geom_boxplot(aes(x = model, y = value, fill = metric), 
               alpha = 0.8, 
               outlier.shape = NA, 
               size = 0.6, 
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Pearson" = pal_perf[6],
                               "Spearman" = pal_perf[2])) +
  coord_cartesian(ylim = c(-0.3, 1)) + 
  labs(y = "Performance value", x = "", title = "C.    All models", fill = "", linetype = "") +
  theme_minimal() + 
  theme(
    legend.position = "none",
    title = element_text(size = 15),
    axis.text=element_text(size = 15),
    axis.title=element_text(size = 20),
    legend.text=element_text(size = 10),
    legend.title=element_text(size = 15),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    panel.background = element_rect(fill = "white", colour = "grey50",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

plot_perf_best <- performance_all_best |> 
  dplyr::filter(cat %in% c("Best models"),
                metric %in% c("Pearson", "Spearman")) |> 
  dplyr::mutate(model = forcats::fct_relevel(model, "GLM", "GAM", "SPAMM", "GBM", "SPRF", "RF")) |> 
  ggplot() +
  geom_boxplot(aes(x = model, y = value, fill = metric), 
               alpha = 0.8, 
               outlier.shape = NA, 
               size = 0.6, 
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Pearson" = pal_perf[6],
                               "Spearman" = pal_perf[2])) +
  coord_cartesian(ylim = c(-0.3, 1)) + 
  labs(y = "Performance value", x = "", title = "Best models", fill = "", linetype = "") +
  theme_minimal() + 
  theme(
    legend.position = c(0.7, 0.1),
    legend.direction = "horizontal",
    title = element_text(size = 15),
    axis.text=element_text(size = 15),
    axis.title=element_text(size = 20),
    legend.text=element_text(size = 10),
    legend.title=element_text(size = 15),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    panel.background = element_rect(fill = "white", colour = "grey50",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

plot_perf_best2 <- (best_model + observed_predicted_best_plot) / (plot_perf_all + plot_perf_best)

ggplot2::ggsave("figures/test2.png", plot_perf_best2, height = 13, width = 12)

