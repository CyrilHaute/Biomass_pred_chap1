
library(ggplot2)

load("data/new_raw_data/00_rls_surveys.Rdata")
rls_surveys$survey_id <- as.character(rls_surveys$survey_id)
load("data/new_derived_data/rls_covariates.RData")

rls_points <- rls_surveys |> 
  dplyr::inner_join(rls_covariates, by = "survey_id") |> 
  dplyr::select(latitude, longitude)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

RLS_Map <- ggplot(data = world) +
  geom_sf() +
  scale_x_continuous(limits=c(-180, 180)) +
  scale_y_continuous(limits=c(-70, 70)) +
  geom_sf(fill= "white") +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "blank", size = 0.5), panel.background = element_rect(fill = "azure1", colour = "azure1")) +
  geom_point(data = rls_points, aes(x = longitude, y = latitude, fill="RLS points"), size = 2, color="blue") +
  coord_sf(expand = FALSE) +
  annotate(geom = "text", x = c(-140,-10,80), y = c(10,-45,-35), label = c("Pacific Ocean", "Atlantic Ocean", "Indian Ocean"), color = "grey22", size = 6) +
  theme(legend.position = "none",
        axis.text=element_text(size=10),
        axis.title=element_text(size=15)) +  
  labs(fill = "")

ggsave("figures/RLS_Map.png", RLS_Map, width = 12, height = 5)
