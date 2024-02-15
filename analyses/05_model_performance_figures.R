# script to produce figures that evaluate model performances

# source functions ----

source("R/05_model_performance_functions.R")

# Set palette colors for performance figures

pal_best = PNWColors::pnw_palette("Bay", 6 , type = "continuous")
pal_perf = PNWColors::pnw_palette("Bay",6, type = "continuous")

# select best fitted model for each model type based on a concensus metrics ----

metrics = c('Intercept', 'Slope', 'Pearson', 'Spearman')

# read data

load("outputs/model_assessment_validation/validation.Rdata")
all_assessments_SCV <- model_assessment
all_assessments_SCV <- do.call(rbind, all_assessments_SCV)

# select only the columns to be used later 
all_assessments_SCV <- all_assessments_SCV |> 
  
  dplyr::select(fitted_model, species_name, metrics)

# estimate for each species the best model based on performance metrics  
best_models <- all_assessments_SCV |> 
  aggregate_metrics(., 
                    metrics = c("Intercept", "Slope", "Pearson", "Spearman")) |>
  # find the best fitting model for each species within each fitted_model
  dplyr::group_by(species_name) |> 
  dplyr::do(best_model = .$fitted_model[which.max(.$discrimination)]) |> 
  tidyr::unnest(cols = c('best_model'))

#### Best Model plot ####

best_models_pr <- best_models |>  
  dplyr::group_by(best_model) |> 
  dplyr::summarise(n = dplyr::n()) |> 
  dplyr::mutate(pr = (n*100)/sum(n))

library(ggplot2)

best_model <- best_models_pr |> 
  # mutate(best_model = fct_relevel(best_model, "GLM", "GAM", "SPAMM", "RF", "GBM", "SPRF")) %>%
  dplyr::mutate(best_model = forcats::fct_relevel(best_model, "GLM", "SPAMM", "RF", "SPRF")) |> 
  ggplot(aes(x = best_model, y = pr, fill = best_model)) +
  geom_bar(width = 0.8, stat = 'identity') +
  scale_fill_manual(values = pal_best) +
  # scale_y_continuous(limits=c(0, 80)) +
  labs(x = "Statistic methods", y = "Best model (%)", fill = "Method", title = "B") +
  theme(title = element_text(size=20),
        axis.text=element_text(size=15),
        axis.title=element_text(size=25),
        legend.text=element_text(size=20), 
        legend.title=element_text(size=25),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# produce histograms of model performance for best models ----

p_level <- unique(all_assessments_SCV$fitted_model)

perf_models_all <- all_assessments_SCV |> 
  dplyr::summarise_at(vars(Intercept:Spearman), list(function(x) list(Q0.05 = round(quantile(x, 0.05, na.rm = T), 2),
                                                                         IQR0.25 = round(quantile(x, 0.25, na.rm = T), 2),
                                                                         median  = round(median(x, na.rm = T), 2),
                                                                         IQR0.75 = round(quantile(x, 0.75, na.rm = T), 2),
                                                                         Q0.95 = round(quantile(x, 0.95, na.rm = T), 2)))) |> 
  dplyr::mutate(summary_value = c('Q0.05','IQR0.25', 'median', 'IQR0.75', 'Q0.95')) |> 
  tidyr::unnest() |> 
  dplyr::group_by(summary_value) |> 
  dplyr::select(Intercept:Spearman) |> 
  t() |> 
  data.frame()

perf_models_all <- perf_models_all |> 
  dplyr::mutate(metric = rownames(perf_models_all)) |>
  dplyr::select(metric, X1:X5)

perf_models_details <- parallel::mclapply(1:length(unique(all_assessments_SCV$fitted_model)), function(i) {
  
  p_level <- p_level[i]
  
  perf_models <- all_assessments_SCV |> 
    dplyr::filter(fitted_model == p_level) |> 
    dplyr::summarise_at(vars(Intercept:Spearman), list(function(x) list(Q0.05 = round(quantile(x, 0.05, na.rm = T), 2),
                                                                           IQR0.25 = round(quantile(x, 0.25, na.rm = T), 2),
                                                                           median  = round(median(x, na.rm = T), 2),
                                                                           IQR0.75 = round(quantile(x, 0.75, na.rm = T), 2),
                                                                           Q0.95 = round(quantile(x, 0.95, na.rm = T), 2)))) |> 
    dplyr::mutate(summary_value = c('Q0.05', 'IQR0.25', 'median', 'IQR0.75', 'Q0.95')) |> 
    tidyr::unnest() |> 
    dplyr::group_by(summary_value) |> 
    dplyr::select(Intercept:Spearman) |> 
    t() |> 
    data.frame() 
  
  perf_models <- perf_models |> 
    dplyr::mutate(metric = rownames(perf_models)) |> 
    dplyr::select(metric, X1:X5)
  
}, mc.cores = 1)
names(perf_models_details) <- p_level

# Select only the best model for each species

best_assessments_SCV <- dplyr::inner_join(all_assessments_SCV, best_models, by = "species_name")
best_assessments_SCV <- best_assessments_SCV[best_assessments_SCV$fitted_model == best_assessments_SCV$best_model,]

perf_models_all_best <- best_assessments_SCV |> 
  dplyr::summarise_at(vars(Intercept:Spearman), list(function(x) list(Q0.05 = round(quantile(x, 0.05, na.rm = T), 2),
                                                                  IQR0.25 = round(quantile(x, 0.25, na.rm = T), 2),
                                                                  median  = round(median(x, na.rm = T), 2),
                                                                  IQR0.75 = round(quantile(x, 0.75, na.rm = T), 2),
                                                                  Q0.95 = round(quantile(x, 0.95, na.rm = T), 2)))) |> 
  dplyr::mutate(summary_value = c('Q0.05','IQR0.25', 'median', 'IQR0.75', 'Q0.95')) |> 
  tidyr::unnest() |> 
  dplyr::group_by(summary_value) |> 
  dplyr::select(Intercept:Spearman) |> 
  t() |> 
  data.frame() 

perf_models_all_best <- perf_models_all_best |> 
  dplyr::mutate(metric = rownames(perf_models_all_best)) |> 
  dplyr::select(metric, X1:X5)

# Manage data for performance plot

# performance_all <- all_assessments_SCV[,c(1,2,4:7)]
performance_all <- dplyr::tibble(species_name = rep(all_assessments_SCV$species_name, 4),
                                 value = c(all_assessments_SCV$Intercept, all_assessments_SCV$Slope, all_assessments_SCV$Pearson, all_assessments_SCV$Spearman),
                          metrics = c(rep("Intercept", nrow(all_assessments_SCV)), 
                                      rep("Slope", nrow(all_assessments_SCV)),
                                      rep("Pearson", nrow(all_assessments_SCV)), 
                                      rep("Spearman", nrow(all_assessments_SCV))),
                          model = rep(all_assessments_SCV$fitted_model, 4))

performance_all[performance_all$metrics == "Intercept",2] <- log10(performance_all[performance_all$metrics == "Intercept",2]$value + 1)
performance_all[performance_all$metrics == "Slope",2] <- log10(performance_all[performance_all$metrics == "Slope",2]$value + 1)

performance_best <- dplyr::tibble(species_name = rep(best_assessments_SCV$species_name, 4),
                          value = c(best_assessments_SCV$Intercept, best_assessments_SCV$Slope, best_assessments_SCV$Pearson, best_assessments_SCV$Spearman),
                          metrics = c(rep("Intercept", nrow(best_assessments_SCV)), 
                                      rep("Slope", nrow(best_assessments_SCV)),
                                      rep("Pearson", nrow(best_assessments_SCV)), 
                                      rep("Spearman", nrow(best_assessments_SCV))),
                          best_model = rep(best_assessments_SCV$fitted_model, 4))

performance_best[performance_best$metrics == "Intercept",2] <- log10(performance_best[performance_best$metrics == "Intercept",2]$value + 1)
performance_best[performance_best$metrics == "Slope",2] <- log10(performance_best[performance_best$metrics == "Slope",2]$value + 1)

performance_all_best <- dplyr::full_join(performance_all, performance_best)
performance_all_best$cat <- NA

performance_all_best[which(performance_all_best$model == performance_all_best$best_model),6] <- "Best models"
performance_all_best[which(is.na(performance_all_best$cat) == TRUE),6] <- "All models"

plot_pearson <- performance_plot(performance_all_best,
                                 metrics_sel = "Pearson",
                                 slope = 0,
                                 intercept = 1,
                                 color = pal_perf,
                                 ylim = c(-0.5, 1),
                                 legend.position = 'none',
                                 plot_title = "")

plot_spearman <- performance_plot(performance_all_best,
                                  metrics_sel = "Spearman",
                                  slope = 0,
                                  intercept = 1,
                                  color = pal_perf,
                                  ylim = c(-0.5, 1),
                                  legend.position = c(0.5, 0.1),
                                  plot_title = "")

plot_slope <- performance_plot(performance_all_best,
                               metrics_sel = "Slope",
                               slope = 0,
                               intercept = 0.30103,
                               color = pal_perf,
                               ylim = c(-0.05, 0.1),
                               legend.position = "none",
                               plot_title = "")

plot_intercept <- performance_plot(performance_all_best,
                                   metrics_sel = "Intercept",
                                   slope = 0,
                                   intercept = 0,
                                   color = pal_perf,
                                   ylim = c(-1,10),
                                   legend.position = 'none',
                                   plot_title = "A")

all_plots <- patchwork::wrap_plots(plot_intercept, plot_slope, plot_pearson, plot_spearman)
all_plots <- all_plots / best_model

ggplot2::ggsave("figures/plot_perf_best.pdf", all_plots, height = 10, width = 10)
