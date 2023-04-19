library(RColorBrewer)
library(dplyr)
library(viridis)
library(raster)
library(sf)
library(readxl)
library(lwgeom)
library(rgeos)
library(geojsonsf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)


### pts_RLS_global = fichier excel contenant l'ensemble des points RLS
pts_RLS_global <- readRDS("../Biomass_prediction/data/Cyril_data/RLS_sitesInfos.rds")
### ERboudaries = Shapefile, couche QGIS avec l'ensemble des limites des écorégions
ERboundaries <- st_read ("../hab_geo/boundaries/Boudaries.shp")

  ### fonction st_as_sf : va transformer ton fichier excel pts_RLS_global en objet SIG (points)
pts_RLS_global<- st_as_sf(x = pts_RLS_global, coords = c("SiteLongitude" , "SiteLatitude"), crs = "WGS84")
pts_RLS_global<- pts_RLS_global %>% filter(SurveyID %in% Fish_RLS$SurveyID)

world <- ne_countries(scale = "medium", returnclass = "sf")

### fonction st_intersection = intersection points RLS x boundaries écorégions.
sf::sf_use_s2(FALSE)
pts_RLS <- st_intersection(pts_RLS_global, ERboundaries)
st_write (pts_RLS, "../Biomass_prediction/data/Cyril_data",driver= "ESRI Shapefile")
pts_RLS$SurveyID <- as.factor(pts_RLS$SurveyID)
pts_RLS$Longitude <- sapply(1:nrow(pts_RLS), function(i){ pts_RLS$geometry[[i]][1]})
pts_RLS$Latitude <- sapply(1:nrow(pts_RLS), function(i){ pts_RLS$geometry[[i]][2]})


world <- ne_countries(scale = "medium", returnclass = "sf")
RLS_Map <- ggplot(data = world) +
             geom_sf() +
             scale_x_continuous(limits=c(-180, 180)) +
             scale_y_continuous(limits=c(-70, 70)) +
             geom_sf(fill= "white") +
             theme(panel.grid.major = element_line(color = gray(.5), linetype = "blank", size = 0.5), panel.background = element_rect(fill = "azure1", colour = "azure1")) +
             geom_point(data = pts_RLS, aes(x=Longitude, y=Latitude, fill="RLS points"), size = 2, color="blue") +
             coord_sf(expand = FALSE) +
             annotate(geom = "text", x = c(-140,-10,80), y = c(10,-45,-35), label = c("Pacific Ocean", "Atlantic Ocean", "Indian Ocean"), color = "grey22", size = 6) +
             theme(legend.position = "none",
                   axis.text=element_text(size=10),
                   axis.title=element_text(size=15)) +  
  labs(fill = "")


ggsave("RLS_Map.png", RLS_Map, width = 10, height = 10)
  
