# Evaluate isolation by sea distance for Atlantic herring.
# Uses pair-wise FST estimates and shortest marine distances.

library(sf)
library(terra)
library(raster)
library(gdistance)
library(rnaturalearth)
library(sp)
library(data.table)
library(vegan)

# Read pair-wise FST estimates and linearize them.
fst_matrix_file <- "paired_Fst.csv"
fst_matrix <- fread(fst_matrix_file)

linearized_fst <- as.dist(fst_matrix / (1 - fst_matrix))

# Read sampling localities.
localities <- fread("SNP_samples_herring.csv")

# Convert localities from WGS84 to ETRS89 / LAEA Europe.
localities_sf <- sf::st_as_sf(
  localities,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

localities_3035 <- sf::st_transform(
  localities_sf,
  3035
)

# Define a study area extending 250 km beyond the locality bbox.
study_area <- sf::st_as_sfc(
  sf::st_bbox(localities_3035)
)

study_area <- sf::st_buffer(
  study_area,
  dist = 250000
)

# Obtain Natural Earth land polygons and repair invalid geometries.
s2_previous <- sf::sf_use_s2()
sf::sf_use_s2(FALSE)

land <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf"
)

land <- sf::st_make_valid(land)

sf::sf_use_s2(s2_previous)

# Project land polygons and crop them to the study area.
land_3035 <- sf::st_transform(
  land,
  3035
)

land_crop <- suppressWarnings(
  sf::st_crop(
    land_3035,
    sf::st_bbox(study_area)
  )
)

land_crop <- sf::st_union(land_crop)

# Build a 2-km raster template in a metric projection.
resolution_m <- 2000

bbox_projected <- sf::st_bbox(study_area)

template <- terra::rast(
  xmin = bbox_projected["xmin"],
  xmax = bbox_projected["xmax"],
  ymin = bbox_projected["ymin"],
  ymax = bbox_projected["ymax"],
  resolution = resolution_m,
  crs = "EPSG:3035"
)

# Rasterize land and define sea cells as 1 and land cells as NA.
land_raster <- terra::rasterize(
  terra::vect(land_crop),
  template,
  field = 1,
  background = NA
)

sea_raster <- terra::ifel(
  is.na(land_raster),
  1,
  NA
)

# Convert to RasterLayer for gdistance.
sea <- raster::raster(sea_raster)

# Extract projected coordinates of sampling localities.
xy <- sf::st_coordinates(localities_3035)

points_sp <- sp::SpatialPoints(
  xy,
  proj4string = raster::crs(sea)
)

# Detect points falling on land because of raster discretization.
idx_na <- which(
  is.na(
    raster::extract(
      sea,
      points_sp
    )
  )
)

# Move affected points to the nearest valid sea cell.
if (length(idx_na) > 0) {
  for (i in idx_na) {
    distance_raster <- raster::distanceFromPoints(
      sea,
      xy[i, , drop = FALSE]
    )
    
    distance_raster[
      is.na(sea[])
    ] <- NA
    
    nearest_cell <- which.min(
      distance_raster[]
    )
    
    xy[i, ] <- raster::xyFromCell(
      sea,
      nearest_cell
    )
  }
}

# Rebuild spatial points using adjusted coordinates.
points_sp <- sp::SpatialPoints(
  xy,
  proj4string = raster::crs(sea)
)

# Build a uniform marine resistance network.
transition_function <- function(x) {
  1 / mean(x, na.rm = TRUE)
}

sea_transition <- gdistance::transition(
  sea,
  transitionFunction = transition_function,
  directions = 8
)

# Correct transition costs for distance between cell centers.
sea_transition <- gdistance::geoCorrection(
  sea_transition,
  type = "c",
  multpl = FALSE,
  scl = FALSE
)

# Calculate shortest marine distances among all locality pairs.
sea_distance_m <- as.matrix(
  gdistance::costDistance(
    sea_transition,
    points_sp,
    points_sp
  )
)

sea_distance_km <- sea_distance_m / 1000

# Assign sample names while preserving locality order.
rownames(sea_distance_km) <- localities$sample
colnames(sea_distance_km) <- localities$sample

sea_distance_km <- as.dist(sea_distance_km)

# Test isolation by sea distance with a Pearson Mantel test.
mantel_result <- vegan::mantel(
  sea_distance_km,
  linearized_fst,
  method = "pearson",
  permutations = 9999
)

# Calculate a two-tailed permutation P value.
p_two_tailed <- (
  sum(
    abs(mantel_result$perm) >=
      abs(mantel_result$statistic)
  ) + 1
) / (
  length(mantel_result$perm) + 1
)

print(mantel_result$statistic)
print(p_two_tailed)
