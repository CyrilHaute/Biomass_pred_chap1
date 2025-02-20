# This script generate partial dependence plot of model performance

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

pearson_imp_var <- ggplot(var_imp |> 
                         dplyr::filter(peformance == "Pearson")) +
  geom_col(aes(x = reorder(covariates, value), y = value, fill = cov_type)) +
  scale_fill_manual(values = c("Traits" = pal[2],
                               "Taxonomy" = pal[4],
                               "Model" = pal[6])) + 
  theme_classic() +
  coord_flip() +
  tidytext::scale_x_reordered() +
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 6),
        strip.text.x = element_text(size = 5),
        strip.text.y = element_text(size = 5),
        legend.position = "none",
        title = element_text(size = 6)) + 
  labs(x = "", y = "Importance", title = "A.     Pearson r² = 0.33", fill = "")
spearman_imp_var <- ggplot(var_imp |> 
                            dplyr::filter(peformance == "Spearman")) +
  geom_col(aes(x = reorder(covariates, value), y = value, fill = cov_type)) +
  scale_fill_manual(values = c("Traits" = pal[2],
                               "Taxonomy" = pal[4],
                               "Model" = pal[6])) + 
  theme_classic() +
  coord_flip() +
  tidytext::scale_x_reordered() +
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 6),
        strip.text.x = element_text(size = 5),
        strip.text.y = element_text(size = 5),
        legend.position = "none",
        title = element_text(size = 6)) + 
  labs(x = "", y = "Importance", title = "B.     Spearman r² = 0.59", fill = "")

# yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata)$predictions) 
# ALEPlot::ALEPlot(as.data.frame(performance_bind), pearson_model, pred.fun=yhat, J = 8, K = 100,
#                  NA.plot = TRUE)
# https://christophm.github.io/interpretable-ml-book/pdp.html#disadvantages-5

patial_pear_count <- pdp::partial(pearson_model, train = performance_bind, pred.var = c("count")) |> 
  dplyr::rename(pearson = yhat)
patial_pear_ml <- pdp::partial(pearson_model, train = performance_bind, pred.var = c("MaxLength")) |> 
  dplyr::rename(pearson = yhat)

patial_pear_count_plot <- patial_pear_count |> 
  ggplot() +
  geom_line(aes(x = count, y = pearson), size = 0.8, color = pal[1]) +
  labs(y = "Pearson", x = "Number of occurrences") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 4),
        strip.text.x = element_text(size = 4),
        strip.text.y = element_text(size = 4),
        legend.position = "none")
patial_pear_ml_plot <- patial_pear_ml |> 
  ggplot() +
  geom_line(aes(x = MaxLength, y = pearson), size = 0.8, color = pal[1]) +
  labs(y = "Pearson", x = "Maximum body length") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 4),
        strip.text.x = element_text(size = 4),
        strip.text.y = element_text(size = 4),
        legend.position = "none")
patial_pear_plot <- patchwork::wrap_plots(patial_pear_count_plot / patial_pear_ml_plot)

pearson_plot <- pearson_imp_var + 
  patchwork::inset_element(patial_pear_plot, left = 0.31, bottom = 0.01, right = 1, top = 0.5)


patial_spear_count <- pdp::partial(spearman_model, train = performance_bind, pred.var = c("count")) |> 
  dplyr::rename(spearman = yhat)
patial_spear_ml <- pdp::partial(spearman_model, train = performance_bind, pred.var = c("MaxLength")) |> 
  dplyr::rename(spearman = yhat)

patial_spear_count_plot <- patial_spear_count |> 
  ggplot() +
  geom_line(aes(x = count, y = spearman), size = 0.8, color = pal[1]) +
  labs(y = "Spearman", x = "Number of occurrences") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 4),
        strip.text.x = element_text(size = 4),
        strip.text.y = element_text(size = 4),
        legend.position = "none")
patial_spear_ml_plot <- patial_spear_ml |> 
  ggplot() +
  geom_line(aes(x = MaxLength, y = spearman), size = 0.8, color = pal[1]) +
  labs(y = "Spearman", x = "Maximum body length") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 4),
        strip.text.x = element_text(size = 4),
        strip.text.y = element_text(size = 4),
        legend.position = "none")
patial_spear_plot <- patchwork::wrap_plots(patial_spear_count_plot / patial_spear_ml_plot)

spearman_plot <- spearman_imp_var + 
  patchwork::inset_element(patial_spear_plot, left = 0.31, bottom = 0.01, right = 1, top = 0.5)

pearson_spearman <- patchwork::wrap_plots(pearson_plot + spearman_plot)

ggsave("figures/imp_patial_plot.pdf", pearson_spearman, width = 15, height = 11, units = "cm")
