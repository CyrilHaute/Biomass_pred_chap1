load("outputs/performance_model.Rdata")
load("data/new_derived_data/species_count.Rdata")
load("data/new_raw_data/RLS_actinopterygii_data.Rdata")

# Load species traits
sp_car <- read.csv("data/new_raw_data/Traits_tropical_spp_1906.csv", header = TRUE) |> 
  dplyr::rename(species_name = Species)
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

sp_car <- sp_car |> 
  dplyr::select(species_name, MaxLength, Water.column, Habitat, Trophic_guild_name)

phylo <- RLS_actinopterygii_data |> 
  dplyr::select(species_name, order, family)
phylo <- unique(phylo)

performance_bind <- performance_bind |> 
  dplyr::select(species_name, pearson, spearman, model) |> 
  dplyr::inner_join(sp_car) |> 
  dplyr::inner_join(sp_count) |> 
  dplyr::inner_join(phylo)

pearson_model <- ranger::ranger(x = performance_bind[!colnames(performance_bind) %in% c("pearson", "spearman")],
                                y = unlist(performance_bind[colnames(performance_bind) %in% "pearson"]),
                                num.trees = 1000,
                                importance = "impurity")
spearman_model <- ranger::ranger(x = performance_bind[!colnames(performance_bind) %in% c("pearson", "spearman")],
                                 y = unlist(performance_bind[colnames(performance_bind) %in% "spearman"]),
                                 num.trees = 1000,
                                 importance = "impurity")
var_imp <- data.frame(pearson = pearson_model$variable.importance,
                      spearman = spearman_model$variable.importance,
                      covariates = unique(names(c(pearson_model$variable.importance, spearman_model$variable.importance))))
var_imp[var_imp$covariates == "species_name",]$covariates <- "Species"
var_imp[var_imp$covariates == "model",]$covariates <- "Model"
var_imp[var_imp$covariates == "MaxLength",]$covariates <- "Maximum Length (cm)"
var_imp[var_imp$covariates == "Trophic_guild_name",]$covariates <- "Trophic classes"
var_imp[var_imp$covariates == "count",]$covariates <- "Occurrence"
var_imp[var_imp$covariates == "order",]$covariates <- "Order"
var_imp[var_imp$covariates == "family",]$covariates <- "Family"
var_imp <- var_imp |> 
  tidyr::pivot_longer(c("pearson", "spearman"),
                      names_to = "peformance")

var_imp <- var_imp |> 
  dplyr::mutate(covariates = tidytext::reorder_within(covariates, value, peformance))

library(ggplot2)

perf_imp_var <- ggplot(var_imp) +
  geom_col(aes(x = reorder(covariates, value), y = value)) +
  theme_minimal() +
  coord_flip() +
  facet_wrap(~ peformance, scales = "free_y") +
  tidytext::scale_x_reordered() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(12),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10)) + 
  labs(x = "", y = "Importance")

ggsave("figures/perf_imp_var.pdf", perf_imp_var, width = 9, height = 5)

patial_plot_pear_count <- autoplot(pdp::partial(pearson_model, train = performance_bind, pred.var = c("count"))) +
  labs(y = "pearson", x = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(12),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))
patial_plot_pear_ml <- autoplot(pdp::partial(pearson_model, train = performance_bind, pred.var = c("MaxLength"))) +
  labs(y = "pearson", x = "Maximum Length (cm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(12),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))

patial_plot_spear_count <- autoplot(pdp::partial(spearman_model, train = performance_bind, pred.var = c("count"))) +
  labs(y = "spearman", x = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(12),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))
patial_plot_spear_ml <- autoplot(pdp::partial(spearman_model, train = performance_bind, pred.var = c("MaxLength"))) +
  labs(y = "spearman", x = "Maximum Length (cm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(12),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))

library(patchwork)

patial_plot <- (patial_plot_pear_count + patial_plot_pear_ml) / (patial_plot_spear_count + patial_plot_spear_ml)

ggsave("figures/patial_plot.png", patial_plot, width = 10, height = 7)
