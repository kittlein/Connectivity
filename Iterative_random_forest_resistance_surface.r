# Iteratively estimate an oceanographic resistance surface with Random Forest.
# Uses linearized pair-wise FST and shortest marine paths among localities.

library(data.table)
library(gdistance)
library(randomForestSRC)
library(raster)
library(terra)

# Define inputs and the number of Random Forest iterations.
predictor_file <- "layers.tif"
locality_file <- "SNP_samples_herring.csv"
fst_matrix_file <- "paired_Fst.csv"
n_iterations <- 30L

# Read oceanographic predictor layers.
predictor_stack <- raster::stack(
  predictor_file
)

predictor_raster <- terra::rast(
  predictor_stack
)

valid_cells <- which(
  !is.na(raster::calc(predictor_stack, sum)[])
)

prediction_data <- data.frame(
  raster::getValues(predictor_stack)
)

# Read sampling localities and assign the raster coordinate reference system.
localities <- read.csv(
  locality_file
)

sp::coordinates(localities) <- c(
  "longitude",
  "latitude"
)

sp::proj4string(localities) <- raster::crs(
  predictor_stack
)

# Move localities on NA cells to the nearest valid raster cell.
missing_cells <- which(
  is.na(
    raster::extract(
      predictor_stack[[1]],
      sp::coordinates(localities)
    )
  )
)

original_coordinates <- sp::coordinates(
  localities
)

adjusted_coordinates <- original_coordinates

for (i in missing_cells) {
  distance_raster <- raster::distanceFromPoints(
    predictor_stack[[1]],
    original_coordinates[i, , drop = FALSE]
  )

  distance_raster[
    is.na(predictor_stack[[1]][])
  ] <- NA

  nearest_cell <- which.min(
    distance_raster[]
  )

  adjusted_coordinates[i, ] <- raster::xyFromCell(
    predictor_stack[[1]],
    nearest_cell
  )
}

localities <- sp::SpatialPointsDataFrame(
  coords = adjusted_coordinates,
  data = localities@data,
  proj4string = sp::CRS(
    sp::proj4string(localities)
  )
)

# Read pair-wise FST and linearize it as FST / (1 - FST).
fst_matrix <- data.table::fread(
  fst_matrix_file
)

population_names <- localities$sample

rownames(fst_matrix) <- colnames(fst_matrix)

linearized_fst_matrix <- fst_matrix / (
  1 - fst_matrix
)

diag(linearized_fst_matrix) <- 0

linearized_fst <- as.matrix(linearized_fst_matrix)[upper.tri(linearized_fst_matrix)]

# Convert the resistance surface to transition conductance.
transition_function <- function(x) {
  1 / mean(x, na.rm = TRUE)
}

# Calculate coefficient of determination from observed and predicted values.
calculate_r2 <- function(observed, predicted) {
  1 - sum((observed - predicted)^2) /
    sum((observed - mean(observed))^2)
}

# Calculate shortest paths between all pairs of sampling localities.
calculate_shortest_paths <- function(
  resistance_surface,
  localities
) {
  transition_layer <- gdistance::transition(
    x = resistance_surface,
    transitionFunction = transition_function,
    directions = 16
  )

  transition_layer <- gdistance::geoCorrection(
    transition_layer,
    type = "c"
  )

  n_localities <- nrow(localities)

  paths <- gdistance::shortestPath(
    x = transition_layer,
    origin = sp::coordinates(localities)[1, ],
    goal = sp::coordinates(localities)[-1, ],
    output = "SpatialLines"
  )

  if (n_localities > 2L) {
    for (k in 2:(n_localities - 1L)) {
      new_paths <- gdistance::shortestPath(
        x = transition_layer,
        origin = sp::coordinates(localities)[k, ],
        goal = sp::coordinates(localities)[
          (k + 1L):n_localities,
        ],
        output = "SpatialLines"
      )

      paths <- rbind(
        paths,
        new_paths
      )
    }
  }

  sp::proj4string(paths) <- raster::crs(
    resistance_surface
  )

  paths
}

# Extract mean oceanographic values along each shortest path.
extract_path_predictors <- function(
  paths,
  predictor_raster,
  linearized_fst
) {
  paths_terra <- terra::vect(
    paths
  )

  path_values <- terra::extract(
    x = predictor_raster,
    y = paths_terra,
    fun = mean,
    na.rm = TRUE,
    touches = TRUE
  )

  path_values$ID <- NULL
  path_values$linearized_fst <- as.vector(
    linearized_fst
  )

  path_values
}

# Tune and fit the Random Forest model.
fit_random_forest <- function(path_values) {
  tuning <- randomForestSRC::tune.rfsrc(
    linearized_fst ~ .,
    data = path_values,
    trace = FALSE
  )

  model <- randomForestSRC::rfsrc(
    formula = linearized_fst ~ .,
    data = path_values,
    nodesize = tuning$optimal["nodesize"],
    mtry = tuning$optimal["mtry"]
  )

  model
}

# Initialize a spatially uniform resistance surface.
mean_linearized_fst <- mean(
  linearized_fst,
  na.rm = TRUE
)

resistance_surface <- predictor_stack[[1]]
resistance_surface[] <- mean_linearized_fst

resistance_surface[
  is.na(predictor_stack[[1]][])
] <- NA

names(resistance_surface) <- "linearized_fst"

# Calculate initial shortest paths and fit the first Random Forest.
shortest_paths <- calculate_shortest_paths(
  resistance_surface,
  localities
)

path_values <- extract_path_predictors(
  shortest_paths,
  predictor_raster,
  linearized_fst
)

rf_model <- fit_random_forest(
  path_values
)

rf_prediction <- predict(
  object = rf_model,
  newdata = prediction_data
)

resistance_surface[
  valid_cells
] <- rf_prediction$predicted

r2 <- calculate_r2(
  observed = rf_model$yvar,
  predicted = rf_model$predicted
)

r2_oob <- calculate_r2(
  observed = rf_model$yvar,
  predicted = rf_model$predicted.oob
)

# Store initial results and allocate objects for subsequent iterations.
resistance_surfaces <- vector(
  "list",
  n_iterations + 1L
)

path_predictors <- vector(
  "list",
  n_iterations + 1L
)

shortest_path_list <- vector(
  "list",
  n_iterations + 1L
)

model_performance <- vector(
  "list",
  n_iterations + 1L
)

resistance_surfaces[[1]] <- resistance_surface
path_predictors[[1]] <- path_values
shortest_path_list[[1]] <- shortest_paths

model_performance[[1]] <- data.frame(
  iteration = 0L,
  r2 = r2,
  r2_oob = r2_oob
)

# Recalculate paths and update the resistance surface iteratively.
for (iteration in seq_len(n_iterations)) {
  shortest_paths <- calculate_shortest_paths(
    resistance_surface,
    localities
  )

  path_values <- extract_path_predictors(
    shortest_paths,
    predictor_raster,
    linearized_fst
  )

  rf_model <- fit_random_forest(
    path_values
  )

  rf_prediction <- predict(
    object = rf_model,
    newdata = prediction_data
  )

  resistance_surface[
    valid_cells
  ] <- rf_prediction$predicted

  r2 <- calculate_r2(
    observed = rf_model$yvar,
    predicted = rf_model$predicted
  )

  r2_oob <- calculate_r2(
    observed = rf_model$yvar,
    predicted = rf_model$predicted.oob
  )

  result_index <- iteration + 1L

  resistance_surfaces[[result_index]] <- resistance_surface
  path_predictors[[result_index]] <- path_values
  shortest_path_list[[result_index]] <- shortest_paths

  model_performance[[result_index]] <- data.frame(
    iteration = iteration,
    r2 = r2,
    r2_oob = r2_oob
  )
  
  cat("iter = ", iteration, "rsq = ", r2, ", rsq.oob = ", r2_oob, "\n")
}

# Save all iterations, model performance, and the final resistance surface.
model_performance <- data.table::rbindlist(
  model_performance
)

analysis_results <- list(
  resistance_surfaces = resistance_surfaces,
  path_predictors = path_predictors,
  shortest_paths = shortest_path_list,
  model_performance = model_performance
)

saveRDS(
  analysis_results,
  "Herring_RF_resistance_iterations.rds"
)

data.table::fwrite(
  model_performance,
  "Herring_RF_resistance_model_performance.csv"
)

imax = which.max(model_performance$r2_oob)

raster::writeRaster(
  resistance_surfaces[[imax]],
  filename = "Herring_RF_final_resistance_surface.tif",
  overwrite = TRUE
)

