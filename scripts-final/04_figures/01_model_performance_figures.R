# script to produce figures that evaluate model performances

# load in packages ---- 

libs <- c('tidyverse', 'gridExtra', 'ggplot2', 'patchwork', 'matrixStats', 'parallel','PNWColors', 'agricolae')
lapply(libs, library, character.only = T, lib.loc = '/home/marbec/R/x86_64-pc-linux-gnu-library/4.1')

# check all packages are loaded
if(sum(libs %in% (.packages())) != length(libs)){
  stop('packages not loaded correctly')}

# source functions ----

source('scripts-final/00_functions/model_performance_functions.R')

pal_best = pnw_palette("Bay", 6 , type = "continuous")
pal_perf = pnw_palette("Bay",2, type = "continuous")

# select best fitted model for each model type based on a concensus metrics ----

metrics = c('Amae_rel_mean', 'Intercept', 'Slope', 'Pearson', 'Spearman', 'Pdispersion')

# read and clean assessment data as in the script 01-model-performance-figures.R

all_assessments_SCV <- readRDS("results/model_assessment_all_R3/SCV/validation.rds")
all_assessments_SCV <- do.call(rbind, all_assessments_SCV)

# select only the columns to be used later 
all_assessments_SCV <- all_assessments_SCV %>% 
  
  dplyr::select(fitted_model, species_name, metrics, Evaluation_number, Evaluation_message)

# estimate for each species the best model based on performance metrics  
best_models <- all_assessments_SCV %>% 
  # estimate the relative metric performance within a cross validation and dataset
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(., 
                                                                  metrics = c('Amae_rel_mean', 'Intercept', 'Slope', 
                                                                              'Pearson', 'Spearman', 'Pdispersion')))) %>% 
  .$metric_aggregation %>% 
  do.call(rbind, .) %>% 
  # find the best fitting model for each species within each fitted_model
  group_by(species_name) %>% 
  do(best_model = .$fitted_model[which.max(.$discrimination)]) %>% 
  unnest(cols = c('best_model'))

saveRDS(best_models, file = 'results/overall_best_models_R3.rds')

#### Best Model plot ####

best_models <- readRDS("results/overall_best_models_R3.rds")

best_models_pr <- best_models %>% 
  group_by(best_model) %>% 
  summarise(n = n()) %>% 
  mutate(pr = (n*100)/sum(n))

best_model <- best_models_pr %>%
  mutate(best_model = fct_relevel(best_model, "GLM", "GAM", "SPAMM", "RF", "GBM", "SPRF")) %>%
  ggplot(aes(x = best_model, y = pr, fill = best_model)) +
  geom_bar(width = 0.8, stat = 'identity') +
  scale_fill_manual(values = pal_best) +
  scale_y_continuous(limits=c(0, 30)) +
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

all_assessments_SCV <- readRDS("results/model_assessment_all_R3/SCV/validation.rds")
all_assessments_SCV <- do.call(rbind, all_assessments_SCV)

all_assessments_SCV <- all_assessments_SCV %>% 
  
  dplyr::select(fitted_model, species_name, metrics, Evaluation_number, Evaluation_message)

p_level <- unique(all_assessments_SCV$fitted_model)

perf_models_all <- all_assessments_SCV %>%
  summarise_at(., vars(Amae_rel_mean:Pdispersion), list(function(x) list(Q0.05 = round(quantile(x, 0.05, na.rm = T), 2),
                                                                         IQR0.25 = round(quantile(x, 0.25, na.rm = T), 2),
                                                                         median  = round(median(x, na.rm = T), 2),
                                                                         IQR0.75 = round(quantile(x, 0.75, na.rm = T), 2),
                                                                         Q0.95 = round(quantile(x, 0.95, na.rm = T), 2)))) %>%
  mutate(summary_value = c('Q0.05','IQR0.25', 'median', 'IQR0.75', 'Q0.95')) %>%
  unnest() %>%
  group_by(summary_value) %>%
  dplyr::select(Amae_rel_mean:Pdispersion) %>%
  t() %>%
  data.frame() %>%
  mutate(metric = rownames(.)) %>%
  dplyr::select(metric, X1:X5)
perf_models_all <- perf_models_all[-c(2,7),]

perf_models_details <- mclapply(1:length(unique(all_assessments_SCV$fitted_model)), function(i) {
  p_level <- p_level[i]
  perf_models <- all_assessments_SCV %>%
    filter(fitted_model == p_level) %>% 
    summarise_at(., vars(Amae_rel_mean:Pdispersion), list(function(x) list(Q0.05 = round(quantile(x, 0.05, na.rm = T), 2),
                                                                           IQR0.25 = round(quantile(x, 0.25, na.rm = T), 2),
                                                                           median  = round(median(x, na.rm = T), 2),
                                                                           IQR0.75 = round(quantile(x, 0.75, na.rm = T), 2),
                                                                           Q0.95 = round(quantile(x, 0.95, na.rm = T), 2)))) %>%
    mutate(summary_value = c('Q0.05', 'IQR0.25', 'median', 'IQR0.75', 'Q0.95')) %>%
    unnest() %>%
    group_by(summary_value) %>%
    dplyr::select(Amae_rel_mean:Pdispersion) %>%
    t() %>%
    data.frame() %>%
    mutate(metric = rownames(.)) %>%
    dplyr::select(metric, X1:X5)
  perf_models <- perf_models[-c(2,7),]
},mc.cores = 1)
names(perf_models_details) <- p_level

best_assessments_SCV <- inner_join(all_assessments_SCV, best_models, by = "species_name")
best_assessments_SCV <- best_assessments_SCV[best_assessments_SCV$fitted_model == best_assessments_SCV$best_model,]

# Manage data for performance plot

performance_all <- all_assessments_SCV[,c(1,2,4:7)]
performance_all <- tibble(species_name = rep(performance_all$species_name, 4),
                          value = c(performance_all$Intercept, performance_all$Slope, performance_all$Pearson, performance_all$Spearman),
                          metrics = c(rep("Intercept", nrow(performance_all)), 
                                      rep("Slope", nrow(performance_all)),
                                      rep("Pearson", nrow(performance_all)), 
                                      rep("Spearman", nrow(performance_all))),
                          model = rep(performance_all$fitted_model, 4))

performance_all[performance_all$metrics == "Intercept",2] <- log10(performance_all[performance_all$metrics == "Intercept",2]$value + 1)
performance_all[performance_all$metrics == "Slope",2] <- log10(performance_all[performance_all$metrics == "Slope",2]$value + 1)

performance_best <- best_assessments_SCV[,c(1,2,4:7)]
performance_best <- tibble(species_name = rep(performance_best$species_name, 4),
                          value = c(performance_best$Intercept, performance_best$Slope, performance_best$Pearson, performance_best$Spearman),
                          metrics = c(rep("Intercept", nrow(performance_best)), 
                                      rep("Slope", nrow(performance_best)),
                                      rep("Pearson", nrow(performance_best)), 
                                      rep("Spearman", nrow(performance_best))),
                          best_model = rep(performance_best$fitted_model, 4))

performance_best[performance_best$metrics == "Intercept",2] <- log10(performance_best[performance_best$metrics == "Intercept",2]$value + 1)
performance_best[performance_best$metrics == "Slope",2] <- log10(performance_best[performance_best$metrics == "Slope",2]$value + 1)

performance_best <- performance_best[,c(1,4)]

performance_all_best <- inner_join(performance_all, performance_best, by = "species_name")
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
                                  legend.position = 'none',
                                  plot_title = "")

plot_slope <- performance_plot(performance_all_best,
                               metrics_sel = "Slope",
                               slope = 0,
                               intercept = 0.30103,
                               color = pal_perf,
                               ylim = c(-0.058, 0.302),
                               legend.position = c(0.75, 0.73),
                               plot_title = "")

plot_intercept <- performance_plot(performance_all_best,
                                   metrics_sel = "Intercept",
                                   slope = 0,
                                   intercept = 0,
                                   color = pal_perf,
                                   ylim = c(-1,8.5),
                                   legend.position = 'none',
                                   plot_title = "A")

all_plots <- wrap_plots(plot_intercept, plot_slope, plot_pearson, plot_spearman)
all_plots <- all_plots / best_model

ggsave("figures-R3/plot_perf_best.png", all_plots, height = 15, width = 11)
