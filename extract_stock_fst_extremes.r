# ==============================================================================
# Identify minimum and maximum predicted FST cells within herring stock polygons
# ==============================================================================

library(terra)
library(sf)
library(openxlsx)


# Stock polygons
stock_files <- c(
  "her.27.1-24a514a.kml",
  "her.27.20-24.kml",
  "her.27.25-2932.kml",
  "her.27.28.kml",
  "her.27.3031.kml",
  "her.27.3a47d.kml",
  "her.27.5a.kml",
  "her.27.irls.kml",
  "her.27.nirs.kml"
)


# Predicted FST surface from the iteration with the highest OOB R2
fst_surface <- terra::rast("Herring_RF_final_resistance_surface.tif")

names(fst_surface) <- "FST"


# Output tables
min_fst_locations <- data.frame()
max_fst_locations <- data.frame()


# Process each stock polygon
for (stock_file in stock_files) {
  
  message("Processing: ", stock_file)
  
  
  # Read and project polygon
  stock_polygon <- sf::st_read(
    paste0("stocks/",stock_file),
    quiet = TRUE
  )
  
  stock_polygon <- sf::st_make_valid(
    stock_polygon
  )
  
  stock_polygon <- sf::st_transform(
    stock_polygon,
    terra::crs(fst_surface)
  )
  
  stock_vect <- terra::vect(
    stock_polygon
  )
  
  
  # Crop and mask the FST surface
  fst_stock <- terra::crop(
    fst_surface,
    stock_vect
  )
  
  fst_stock <- terra::mask(
    fst_stock,
    stock_vect
  )
  
  
  # Extract valid raster values
  fst_values <- terra::values(
    fst_stock,
    mat = FALSE
  )
  
  valid_cells <- which(
    is.finite(fst_values)
  )
  
  
  if (length(valid_cells) == 0L) {
    
    warning(
      "No valid raster cells found inside ",
      stock_file
    )
    
    next
  }
  
  
  # Minimum FST cell
  min_cell_local <- valid_cells[
    which.min(
      fst_values[valid_cells]
    )
  ]
  
  min_xy <- terra::xyFromCell(
    fst_stock,
    min_cell_local
  )
  
  min_cell_global <- terra::cellFromXY(
    fst_surface,
    min_xy
  )
  
  
  min_fst_locations <- rbind(
    min_fst_locations,
    data.frame(
      stock = sub("\\.kml$", "", stock_file),
      fst = fst_values[min_cell_local],
      cell = min_cell_global,
      longitude = min_xy[1],
      latitude = min_xy[2]
    )
  )
  
  
  # Maximum FST cell
  max_cell_local <- valid_cells[
    which.max(
      fst_values[valid_cells]
    )
  ]
  
  max_xy <- terra::xyFromCell(
    fst_stock,
    max_cell_local
  )
  
  max_cell_global <- terra::cellFromXY(
    fst_surface,
    max_xy
  )
  
  
  max_fst_locations <- rbind(
    max_fst_locations,
    data.frame(
      stock = sub("\\.kml$", "", stock_file),
      fst = fst_values[max_cell_local],
      cell = max_cell_global,
      longitude = max_xy[1],
      latitude = max_xy[2]
    )
  )
}


# Inspect results
print(min_fst_locations)
print(max_fst_locations)


# Save minimum and maximum FST locations
openxlsx::write.xlsx(
  list(
    min_FST = min_fst_locations,
    max_FST = max_fst_locations
  ),
  file = "herring_stock_FST_extreme_locations.xlsx",
  overwrite = TRUE
)
