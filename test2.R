# boat_csv <- read.csv("industrial_vessels_v20240102.csv", header = TRUE, sep = ",")
# save(boat_csv, file = "vessels_rdata.RData")
load("vessels_rdata.RData")

# Keep only fishing boat
fishing_total <- vessels_rdata |> 
  dplyr::mutate(category = ifelse(matched_category == "unmatched" & fishing_score >= 0.8, "unmatched_fishing",
                                  #If the match is unknown, we consider it as a fishing boat if fishing score > 0.8
                                  ifelse(matched_category == "matched_unknown" & fishing_score >= 0.8, "matched_fishing",
                                         # #If it is fishing boat but the length is higher than the longest fishing boat in the world then matched_nonfishing
                                         ifelse(matched_category == "matched_fishing" & length_m > 228, "matched_unknown", matched_category)))) |>
  #Also if unmatched_fishing and length higher than 80% quantile then delete it
  dplyr::mutate(category = ifelse(category == "unmatched_fishing" & length_m > quantile(length_m, 0.8),"matched_unknown",category)) |>
  dplyr::filter(category %in% c("unmatched_fishing", "matched_fishing"),
                lon > 135 & lon < 175,
                lat > -35 & lat < 0)

test <- terra::vect(fishing_total, geom = c("lon", "lat"), crs = "WGS84")

terra::writeVector(test, file = "fishing_shapefile.shp")
