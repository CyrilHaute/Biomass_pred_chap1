
test <- arrow::read_feather("raster_5th_degree.feather") |> 
  as.data.frame() |> 
  dplyr::mutate(lon_index = lon_index/5,
                lat_index = lat_index/5) |> 
  dplyr::rename(x = lon_index,
                y = lat_index) |> 
  dplyr::relocate(x, .before = y) |> 
  dplyr::select(x, y, AIS_fishing, dark_fishing) |> 
  terra::vect(geom = c("x", "y"), crs = "WGS84")

terra::writeVector(test, "fishing_boat.shp")


test <- arrow::read_feather("raster_5th_degree.feather") |> 
  as.data.frame() |> 
  dplyr::mutate(lon_index = lon_index / 5,
                lat_index = lat_index / 5) |> 
  dplyr::rename(x = lon_index,
                y = lat_index) |>
  dplyr::relocate(x, .before = y) |> 
  dplyr::mutate(detection_type = ifelse(dark_fishing < 0.7, "not_publicy_tracked", "publicy_tracked")) |> 
  terra::rast(type = "xyz", crs = "WGS84")

terra::writeRaster(test, "fishing_boat.tif", overwrite = TRUE)
test |> dplyr::group_by(detection_type) |> dplyr::summarise(n=dplyr::n())
