# Extract annual oceanographic extremes at stock-specific FST cells

# Packages

library(terra)
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(stringr)


# Input files

coords_file <- "herring_stock_FST_extreme_locations.xlsx"
archive_file <- "ocean_monthly.7z"

stopifnot(file.exists(coords_file), file.exists(archive_file))

archive_path <- normalizePath(
  archive_file, winslash = "/", mustWork = TRUE
)

sevenzip <- Sys.which("7z")

if (!nzchar(sevenzip)) {
  stop("7-Zip executable not found.")
}


# Oceanographic variables

variables <- tibble::tribble(
  ~variable, ~nc_var,  ~nc_file,
  "chla",    "chla",   "chla_monthly_1993-2026-herring.nc",
  "sss",     "so",     "sss_monthly_1993-2026-herring.nc",
  "sst",     "thetao", "sst_monthly_1993-2026-herring.nc",
  "uo",      "uo",     "uo_monthly_1993-2026-herring.nc",
  "vo",      "vo",     "vo_monthly_1993-2026-herring.nc"
)


# Read minimum- and maximum-FST locations

min_fst_points <- readxl::read_excel(
  coords_file, sheet = "min_FST"
) |>
  dplyr::transmute(
    stock = as.character(stock),
    fst_location = "minFST",
    predicted_fst = as.numeric(fst),
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
  )

max_fst_points <- readxl::read_excel(
  coords_file, sheet = "max_FST"
) |>
  dplyr::transmute(
    stock = as.character(stock),
    fst_location = "maxFST",
    predicted_fst = as.numeric(fst),
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
  )

fst_points <- dplyr::bind_rows(min_fst_points, max_fst_points)

stopifnot(
  nrow(min_fst_points) == 9L,
  nrow(max_fst_points) == 9L,
  nrow(fst_points) == 18L
)

if (any(!is.finite(fst_points$longitude)) ||
    any(!is.finite(fst_points$latitude))) {
  stop("Invalid coordinates in FST location file.")
}


# Create spatial points

points_ll <- terra::vect(
  fst_points, geom = c("longitude", "latitude"), crs = "EPSG:4326"
)


# Extract one NetCDF temporarily from the 7z archive

extract_nc_from_archive <- function(nc_file) {
  
  temp_nc_dir <- file.path(tempdir(), "herring_ocean_nc")
  
  if (!dir.exists(temp_nc_dir)) {
    dir.create(temp_nc_dir, recursive = TRUE)
  }
  
  nc_path <- file.path(temp_nc_dir, nc_file)
  
  if (file.exists(nc_path)) {
    unlink(nc_path)
  }
  
  message("Extracting: ", nc_file)
  
  status <- system2(
    sevenzip,
    args = c(
      "e",
      "-y",
      shQuote(paste0("-o", temp_nc_dir)),
      shQuote(archive_path),
      shQuote(nc_file)
    )
  )
  
  if (status != 0L || !file.exists(nc_path)) {
    stop("Unable to extract ", nc_file, " from ", archive_file, ".")
  }
  
  nc_path
}


# Read one NetCDF variable

read_nc_variable <- function(nc_path, varname) {
  
  message("Opening ", basename(nc_path), " [", varname, "]")
  
  r <- try(
    terra::rast(nc_path, subds = varname),
    silent = TRUE
  )
  
  if (inherits(r, "try-error")) {
    
    nc_gdal <- paste0('NETCDF:"', nc_path, '":', varname)
    
    r <- try(
      terra::rast(nc_gdal),
      silent = TRUE
    )
  }
  
  if (inherits(r, "try-error")) {
    stop(
      "Unable to read variable ", varname,
      " from ", basename(nc_path), "."
    )
  }
  
  names(r) <- paste0("layer_", seq_len(terra::nlyr(r)))
  
  r
}


# Retrieve NetCDF dates

get_nc_dates <- function(r) {
  
  dates <- terra::time(r)
  
  if (is.null(dates) || length(dates) == 0L || all(is.na(dates))) {
    dates <- seq(
      from = as.Date("1993-01-01"),
      to = as.Date("2026-12-01"),
      by = "month"
    )
  }
  
  if (length(dates) != terra::nlyr(r)) {
    stop(
      "Number of dates does not match raster layers: ",
      length(dates), " dates and ", terra::nlyr(r), " layers."
    )
  }
  
  as.Date(dates)
}


# Extract monthly values for one variable

extract_monthly_variable <- function(variable, nc_var, nc_file) {
  
  message("\nProcessing ", variable)
  
  nc_path <- extract_nc_from_archive(nc_file)
  
  on.exit(
    {
      if (file.exists(nc_path)) {
        unlink(nc_path)
      }
    },
    add = TRUE
  )
  
  r <- read_nc_variable(
    nc_path = nc_path, varname = nc_var
  )
  
  dates <- get_nc_dates(r)
  
  message(
    "  ", terra::nlyr(r), " layers: ",
    min(dates), " to ", max(dates)
  )
  
  if (terra::xmin(r) >= 0 && terra::xmax(r) > 180) {
    r <- terra::rotate(r)
  }
  
  raster_crs <- terra::crs(r)
  
  if (is.na(raster_crs) || !nzchar(raster_crs)) {
    stop("Raster for ", variable, " has no CRS.")
  }
  
  points_r <- terra::project(points_ll, raster_crs)
  
  ext <- terra::extract(
    r, points_r, method = "simple", ID = TRUE
  )
  
  ext$stock <- fst_points$stock[ext$ID]
  ext$fst_location <- fst_points$fst_location[ext$ID]
  ext$predicted_fst <- fst_points$predicted_fst[ext$ID]
  
  ext_long <- ext |>
    tidyr::pivot_longer(
      cols = dplyr::starts_with("layer_"),
      names_to = "layer",
      values_to = "value"
    ) |>
    dplyr::mutate(
      layer_id = as.integer(stringr::str_extract(layer, "\\d+$")),
      date = dates[layer_id],
      year = lubridate::year(date),
      month = lubridate::month(date),
      variable = variable
    ) |>
    dplyr::filter(year >= 1993, year <= 2025) |>
    dplyr::select(
      stock, fst_location, predicted_fst, variable,
      date, year, month, value
    )
  
  rm(r)
  gc()
  
  ext_long
}

# Extract monthly values for all variables

monthly_values <- purrr::pmap_dfr(
  variables,
  function(variable, nc_var, nc_file) {
    extract_monthly_variable(
      variable = variable,
      nc_var = nc_var,
      nc_file = nc_file
    )
  }
)


# Calculate annual minima and maxima

annual_extremes <- monthly_values |>
  dplyr::group_by(
    stock, fst_location, predicted_fst, variable, year
  ) |>
  dplyr::summarise(
    n_months = sum(is.finite(value)),
    annual_min = if (all(is.na(value))) {
      NA_real_
    } else {
      min(value, na.rm = TRUE)
    },
    annual_max = if (all(is.na(value))) {
      NA_real_
    } else {
      max(value, na.rm = TRUE)
    },
    .groups = "drop"
  )


# Check years with incomplete monthly data

missing_months <- annual_extremes |>
  dplyr::filter(n_months < 12L) |>
  dplyr::arrange(stock, fst_location, variable, year)

if (nrow(missing_months) > 0L) {
  warning(
    nrow(missing_months),
    " stock-location-variable-year combinations contain fewer ",
    "than 12 valid months."
  )
} else {
  message("All annual summaries contain 12 valid monthly observations.")
}


# Convert annual extrema to long format

annual_extremes_long <- annual_extremes |>
  tidyr::pivot_longer(
    cols = c(annual_min, annual_max),
    names_to = "annual_summary",
    values_to = "value"
  ) |>
  dplyr::mutate(
    annual_summary = dplyr::recode(
      annual_summary,
      annual_min = "min",
      annual_max = "max"
    ),
    predictor = paste(
      annual_summary, variable, fst_location, sep = "_"
    )
  )


# Create wide predictor table

annual_wide <- annual_extremes_long |>
  dplyr::select(stock, year, predictor, value) |>
  tidyr::pivot_wider(
    names_from = predictor,
    values_from = value
  ) |>
  dplyr::arrange(stock, year)


# Order predictors

variable_order <- c("chla", "sst", "sss", "vo", "uo")
summary_order <- c("min", "max")
fst_location_order <- c("minFST", "maxFST")

desired_predictors <- tidyr::expand_grid(
  fst_location = fst_location_order,
  variable = variable_order,
  annual_summary = summary_order
) |>
  dplyr::mutate(
    predictor = paste(
      annual_summary, variable, fst_location, sep = "_"
    )
  ) |>
  dplyr::pull(predictor)

annual_wide <- annual_wide |>
  dplyr::select(
    stock, year, dplyr::any_of(desired_predictors)
  )


# Print diagnostics

cat(
  "\nStocks:", dplyr::n_distinct(annual_wide$stock), "\n"
)

cat(
  "Years:",
  min(annual_wide$year, na.rm = TRUE), "-",
  max(annual_wide$year, na.rm = TRUE), "\n"
)

cat(
  "Predictors per stock-year:",
  ncol(annual_wide) - 2L, "\n"
)


# Create one worksheet per stock

stock_sheets <- split(
  annual_wide,
  annual_wide$stock
)

stock_sheets <- lapply(
  stock_sheets,
  function(x) {
    x |>
      dplyr::select(-stock)
  }
)

names(stock_sheets) <- substr(
  names(stock_sheets), 1, 31
)


# Add supporting worksheets

output_sheets <- c(
  stock_sheets,
  list(
    coordinates = fst_points,
    missing_months = missing_months,
    annual_long = annual_extremes_long
  )
)


# Save results

output_file <- "annual_oceanographic_extremes_1993_2025.xlsx"

openxlsx::write.xlsx(
  output_sheets,
  file = output_file,
  overwrite = TRUE
)

message("\nOutput saved as: ", output_file)