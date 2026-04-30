# Curated workflow for extracting seasonal polygon means from monthly NetCDF files.

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(stringr)

netcdf_sources <- data.frame(
  variable_name = c("chla", "sst", "sss", "vo", "uo"),
  source_name = c("chla", "thetao", "so", "vo", "uo"),
  file_name = c(
    "chla_monthly_1993-2023-herring.nc",
    "sst_sss_monthly_1993-2023-herring.nc",
    "sst_sss_monthly_1993-2023-herring.nc",
    "vo_uo_monthly_1993-2023-herring.nc",
    "vo_uo_monthly_1993-2023-herring.nc"
  ),
  stringsAsFactors = FALSE
)

polygon_files <- paste0("pol", 1:4, "b.kml")
output_file <- "ocean_polygons_seasonal_1993_2023.xlsx"

read_polygons <- function(kml_files) {
  map2_dfr(
    kml_files,
    paste0("pol", seq_along(kml_files)),
    function(file_name, polygon_id) {
      st_read(file_name, quiet = TRUE) |>
        st_make_valid() |>
        mutate(polygon_id = polygon_id)
    }
  ) |>
    st_transform(4326) |>
    st_collection_extract("POLYGON")
}

read_netcdf_variable <- function(netcdf_file, source_name) {
  raster_stack <- try(
    terra::rast(netcdf_file, subds = source_name),
    silent = TRUE
  )

  if (inherits(raster_stack, "try-error")) {
    raster_stack <- terra::rast(paste0("NETCDF:", netcdf_file, ":", source_name))
  }

  names(raster_stack) <- paste0(source_name, "_", seq_len(nlyr(raster_stack)))
  raster_stack
}

get_layer_dates <- function(raster_stack) {
  layer_dates <- terra::time(raster_stack)

  if (is.null(layer_dates) || all(is.na(layer_dates))) {
    layer_dates <- seq(
      from = as.Date("1993-01-01"),
      to = as.Date("2023-12-01"),
      by = "month"
    )
  }

  if (length(layer_dates) != nlyr(raster_stack)) {
    stop(
      "The number of dates does not match the number of raster layers. ",
      "Check whether the NetCDF file contains an additional dimension."
    )
  }

  as.Date(layer_dates)
}

extract_monthly_polygon_means <- function(netcdf_file, source_name, variable_name, polygons_sf) {
  message("Processing variable: ", variable_name)

  raster_stack <- read_netcdf_variable(netcdf_file, source_name)

  if (terra::xmin(raster_stack) >= 0 && terra::xmax(raster_stack) > 180) {
    raster_stack <- terra::rotate(raster_stack)
  }

  layer_dates <- get_layer_dates(raster_stack)
  polygon_vectors <- terra::vect(st_transform(polygons_sf, terra::crs(raster_stack)))

  extracted_values <- terra::extract(
    raster_stack,
    polygon_vectors,
    fun = mean,
    na.rm = TRUE,
    weights = TRUE
  )

  extracted_values$polygon_id <- polygons_sf$polygon_id[extracted_values$ID]

  extracted_values |>
    dplyr::select(-ID) |>
    pivot_longer(
      cols = starts_with(source_name),
      names_to = "layer_name",
      values_to = "mean_value"
    ) |>
    mutate(
      layer_id = as.integer(str_extract(layer_name, "\\d+$")),
      date = layer_dates[layer_id],
      year = lubridate::year(date),
      month = lubridate::month(date),
      variable = variable_name
    ) |>
    dplyr::select(polygon_id, variable, date, year, month, mean_value)
}

polygons_sf <- read_polygons(polygon_files)

monthly_covariates <- pmap_dfr(
  netcdf_sources,
  function(variable_name, source_name, file_name) {
    extract_monthly_polygon_means(
      netcdf_file = file_name,
      source_name = source_name,
      variable_name = variable_name,
      polygons_sf = polygons_sf
    )
  }
)

seasonal_covariates_long <- monthly_covariates |>
  mutate(
    season = case_when(
      month %in% 4:6 ~ "spring",
      month %in% 10:12 ~ "fall",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(season)) |>
  group_by(year, polygon_id, season, variable) |>
  summarise(
    mean_value = mean(mean_value, na.rm = TRUE),
    n_months = sum(!is.na(mean_value)),
    .groups = "drop"
  )

seasonal_covariates_wide <- seasonal_covariates_long |>
  mutate(column_name = paste(variable, polygon_id, season, sep = "_")) |>
  dplyr::select(year, column_name, mean_value) |>
  pivot_wider(
    names_from = column_name,
    values_from = mean_value
  ) |>
  arrange(year)

ordered_columns <- expand.grid(
  variable = c("chla", "sst", "sss", "vo", "uo"),
  polygon_id = paste0("pol", 1:4),
  season = c("spring", "fall"),
  stringsAsFactors = FALSE
) |>
  mutate(column_name = paste(variable, polygon_id, season, sep = "_")) |>
  pull(column_name)

seasonal_covariates_wide <- seasonal_covariates_wide |>
  dplyr::select(year, any_of(ordered_columns))

openxlsx::write.xlsx(seasonal_covariates_wide, output_file)
message("Saved seasonal polygon covariates to: ", output_file)