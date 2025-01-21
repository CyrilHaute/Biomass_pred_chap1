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