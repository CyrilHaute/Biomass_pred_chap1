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
  tidyr::drop_na() |> 
  dplyr::select(species_name, pearson, spearman, model) |> 
  dplyr::inner_join(sp_car) |> 
  dplyr::inner_join(sp_count) |> 
  dplyr::inner_join(phylo) |> 
  dplyr::select(-species_name)

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
var_imp[var_imp$covariates == "model",]$covariates <- "Model"
var_imp[var_imp$covariates == "MaxLength",]$covariates <- "Maximum body length"
var_imp[var_imp$covariates == "Trophic_guild_name",]$covariates <- "Trophic class"
var_imp[var_imp$covariates == "count",]$covariates <- "Number of occurrences"
var_imp[var_imp$covariates == "order",]$covariates <- "Order"
var_imp[var_imp$covariates == "family",]$covariates <- "Family"
var_imp[var_imp$covariates == "Water.column",]$covariates <- "Water column position"
colnames(var_imp)[colnames(var_imp) == "pearson"] <- "Pearson"
colnames(var_imp)[colnames(var_imp) == "spearman"] <- "Spearman"
var_imp$cov_type <- NA
var_imp <- var_imp |> 
  dplyr::mutate(cov_type = dplyr::case_when(covariates %in% c("Order", "Family") ~ "Taxonomy",
                                            covariates %in% c("Maximum body length", "Trophic class", "Number of occurrences", "Habitat", "Water column position") ~ "Traits",
                                            covariates %in% c("Model") ~ "Model"))
var_imp <- var_imp |> 
  tidyr::pivot_longer(c("Pearson", "Spearman"),
                      names_to = "peformance")

var_imp <- var_imp |> 
  dplyr::mutate(covariates = tidytext::reorder_within(covariates, value, peformance))

library(ggplot2)

pal <- PNWColors::pnw_palette("Bay", 6, type = "continuous")

perf_imp_var_pear <- ggplot(var_imp |> 
                              dplyr::filter(peformance == "Pearson")) +
  geom_col(aes(x = reorder(covariates, value), y = value, fill = cov_type)) +
  scale_fill_manual(values = c("Traits" = pal[2],
                               "Taxonomy" = pal[4],
                               "Model" = pal[6])) + 
  theme_classic() +
  coord_flip() +
  tidytext::scale_x_reordered() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        title = element_text(size = 18),
        legend.position = "none") + 
  labs(x = "", y = "Importance", title = "A.     Pearson r² = 0.39", fill = "")
perf_imp_var_spear <- ggplot(var_imp |> 
                               dplyr::filter(peformance == "Spearman")) +
  geom_col(aes(x = reorder(covariates, value), y = value, fill = cov_type)) +
  scale_fill_manual(values = c("Traits" = pal[2],
                               "Taxonomy" = pal[4],
                               "Model" = pal[6])) + 
  theme_classic() +
  coord_flip() +
  tidytext::scale_x_reordered() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        title = element_text(size = 18)) + 
  labs(x = "", y = "Importance", title = "B.     Spearman r² = 0.75", fill = "")

patial_pear_count <- pdp::partial(pearson_model, train = performance_bind, pred.var = c("count")) |> 
  dplyr::rename(pearson = yhat)
patial_spear_count <- pdp::partial(spearman_model, train = performance_bind, pred.var = c("count")) |> 
  dplyr::rename(spearman = yhat)
patial_pear_ml <- pdp::partial(pearson_model, train = performance_bind, pred.var = c("MaxLength")) |> 
  dplyr::rename(pearson = yhat)
patial_spear_ml <- pdp::partial(spearman_model, train = performance_bind, pred.var = c("MaxLength")) |> 
  dplyr::rename(spearman = yhat)
patial_count <- patial_pear_count |> 
  dplyr::inner_join(patial_spear_count) |>
  dplyr::rename(Pearson = pearson,
                Spearman = spearman) |> 
  tidyr::pivot_longer(c("Pearson", "Spearman"),
                      names_to = "performance")
patial_ml <- patial_pear_ml |> 
  dplyr::inner_join(patial_spear_ml) |>
  dplyr::rename(Pearson = pearson,
                Spearman = spearman) |> 
  tidyr::pivot_longer(c("Pearson", "Spearman"),
                      names_to = "performance")

patial_count_plot_pear <- patial_count |> 
  dplyr::filter(performance == "Pearson") |> 
  ggplot() +
  geom_line(aes(x = count, y = value, color = performance), size = 1.2) +
  scale_color_manual(values = c("Pearson" = pal[1])) +
  labs(y = "Pearson", x = "Number of occurrences", color = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 10),
        strip.text.x = element_text(size = 6),
        strip.text.y = element_text(size = 6),
        legend.position = "none")
patial_count_plot_spear <- patial_count |> 
  dplyr::filter(performance == "Spearman") |> 
  ggplot() +
  geom_line(aes(x = count, y = value, color = performance), size = 1.2) +
  scale_color_manual(values = c("Spearman" = pal[1])) +
  labs(y = "Spearman", x = "Number of occurrences", color = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 10),
        strip.text.x = element_text(size = 6),
        strip.text.y = element_text(size = 6),
        legend.position = "none")

patial_ml_plot_pear <- patial_ml |> 
  dplyr::filter(performance == "Pearson") |> 
  ggplot() +
  geom_line(aes(x = MaxLength, y = value, color = performance), size = 1.2) +
  scale_color_manual(values = c("Pearson" = pal[1])) +
  labs(y = "Pearson", x = "Maximum body length", color = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 10),
        strip.text.x = element_text(size = 6),
        strip.text.y = element_text(size = 6),
        legend.position = "none")
patial_ml_plot_spear <- patial_ml |> 
  dplyr::filter(performance == "Spearman") |> 
  ggplot() +
  geom_line(aes(x = MaxLength, y = value, color = performance), size = 1.2) +
  scale_color_manual(values = c("Spearman" = pal[1])) +
  labs(y = "Spearman", x = "Maximum body length", color = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 10),
        strip.text.x = element_text(size = 6),
        strip.text.y = element_text(size = 6),
        legend.position = "none")

library(patchwork)

partial_pear <- patial_count_plot_pear + patial_ml_plot_pear
partial_spear <- patial_count_plot_spear + patial_ml_plot_spear

imp_patial_plot_pear <- perf_imp_var_pear + 
  patchwork::inset_element(partial_pear, left = 0.2, bottom = 0.00, right = 1, top = 0.37)
imp_patial_plot_spear <- perf_imp_var_spear + 
  patchwork::inset_element(partial_spear, left = 0.2, bottom = 0.00, right = 1, top = 0.37)


imp_patial_plot <- imp_patial_plot_pear + imp_patial_plot_spear
ggsave("figures/imp_patial_plot_pear2.png", imp_patial_plot, width = 17.5, height = 10)
