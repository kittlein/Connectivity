# Curated workflow for iterative connectivity-surface modelling.

library(raster)
library(gdistance)
library(randomForestSRC)
library(sp)
library(sf)

predictor_layers <- stack("connectivity_predictors.tif")

non_missing_cells <- which(is.na(calc(predictor_layers, sum)[]) == FALSE)
predictor_table <- data.frame(getValues(predictor_layers)[])

sampling_sites <- read.table("baltic.txt", header = TRUE)
coordinates(sampling_sites) <- c("x", "y")
proj4string(sampling_sites) <- proj4string(predictor_layers)

connectivity_matrix <- readRDS("conectNeutralLoci.rds")

mean_connectivity <- mean(connectivity_matrix, na.rm = TRUE)

connectivity_surface <- predictor_layers[[1]]
connectivity_surface[] <- mean_connectivity
connectivity_surface[which(is.na(predictor_layers[[1]][]) == TRUE)] <- NA

names(connectivity_surface) <- "Connectivity"

mean_transition_function <- function(x) {
  mean(x, na.rm = TRUE)
}

transition_matrix <- transition(
  x = connectivity_surface,
  transitionFunction = mean_transition_function,
  directions = 16
)

transition_matrix <- geoCorrection(transition_matrix, type = "c")

site_index <- 1

least_cost_lines <- shortestPath(
  x = transition_matrix,
  origin = coordinates(sampling_sites)[site_index, ],
  goal = coordinates(sampling_sites)[-site_index, ],
  output = "SpatialLines"
)

for (site_index in 2:length(sampling_sites)) {
  least_cost_lines <- rbind(
    least_cost_lines,
    shortestPath(
      x = transition_matrix,
      origin = coordinates(sampling_sites)[site_index, ],
      goal = coordinates(sampling_sites)[-site_index, ],
      output = "SpatialLines"
    )
  )
}

proj4string(least_cost_lines) <- proj4string(predictor_layers)

least_cost_lines_sf <- sf::st_as_sf(least_cost_lines)

sinusoidal_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +units=m +datum=WGS84 +no_defs"
least_cost_lines_sinusoidal <- sf::st_transform(least_cost_lines_sf, crs = sinusoidal_crs)

valid_connectivity_index <- which(!is.na(connectivity_matrix) == TRUE)

connectivity_values <- as.vector(connectivity_matrix[valid_connectivity_index])
path_distances_km <- sf::st_length(least_cost_lines_sinusoidal) / 1000

sinusoidal_distance_matrix <- matrix(1000, 18, 18)
diag(sinusoidal_distance_matrix) <- NA

valid_distance_index <- which(!is.na(sinusoidal_distance_matrix))
sinusoidal_distance_matrix[valid_distance_index] <- path_distances_km

mantel_test <- vegan::mantel(
  as.dist(sinusoidal_distance_matrix),
  as.dist(connectivity_matrix),
  permutations = 9999
)

path_predictor_values <- raster::extract(
  x = predictor_layers,
  y = least_cost_lines,
  df = TRUE,
  fun = mean,
  na.rm = TRUE
)

path_predictor_values$ID <- NULL
path_predictor_values$Connectivity <- connectivity_values

random_forest_tuning <- tune.rfsrc(Connectivity ~ ., data = path_predictor_values)

random_forest_model <- rfsrc(
  formula = Connectivity ~ .,
  data = path_predictor_values,
  nodesize = random_forest_tuning$optimal["nodesize"],
  mtry = random_forest_tuning$optimal["mtry"]
)

surface_prediction <- predict(object = random_forest_model, newdata = predictor_table)
connectivity_surface[non_missing_cells] <- surface_prediction$predicted

r_squared <- round(MLmetrics::R2_Score(random_forest_model$predicted, random_forest_model$yvar), 3)
r_squared_oob <- round(MLmetrics::R2_Score(random_forest_model$predicted.oob, random_forest_model$yvar), 3)

connectivity_surface_list <- list()
r_squared_list <- list()
least_cost_lines_list <- list()
path_predictor_values_list <- list()

r_squared_list[[1]] <- data.frame(r_squared, r_squared_oob)
connectivity_surface_list[[1]] <- connectivity_surface
least_cost_lines_list[[1]] <- least_cost_lines
path_predictor_values_list[[1]] <- path_predictor_values

for (iteration in 1:5) {
  
  transition_matrix <- transition(
    x = connectivity_surface,
    transitionFunction = mean_transition_function,
    directions = 16
  )
  
  transition_matrix <- geoCorrection(transition_matrix, type = "c")
  
  site_index <- 1
  
  least_cost_lines <- shortestPath(
    x = transition_matrix,
    origin = coordinates(sampling_sites)[site_index, ],
    goal = coordinates(sampling_sites)[-site_index, ],
    output = "SpatialLines"
  )
  
  for (site_index in 2:length(sampling_sites)) {
    least_cost_lines <- rbind(
      least_cost_lines,
      shortestPath(
        x = transition_matrix,
        origin = coordinates(sampling_sites)[site_index, ],
        goal = coordinates(sampling_sites)[-site_index, ],
        output = "SpatialLines"
      )
    )
  }
  
  proj4string(least_cost_lines) <- proj4string(predictor_layers)
  
  path_predictor_values <- raster::extract(
    x = predictor_layers,
    y = least_cost_lines,
    df = TRUE,
    fun = mean,
    na.rm = TRUE
  )
  
  path_predictor_values$ID <- NULL
  path_predictor_values$Connectivity <- connectivity_values
  
  random_forest_tuning <- tune.rfsrc(Connectivity ~ ., data = path_predictor_values)
  
  random_forest_model <- rfsrc(
    formula = Connectivity ~ .,
    data = path_predictor_values,
    importance = TRUE,
    nodesize = random_forest_tuning$optimal["nodesize"],
    mtry = random_forest_tuning$optimal["mtry"]
  )
  
  surface_prediction <- predict(object = random_forest_model, newdata = predictor_table)
  connectivity_surface[non_missing_cells] <- surface_prediction$predicted
  
  r_squared <- round(MLmetrics::R2_Score(random_forest_model$predicted, random_forest_model$yvar), 3)
  r_squared_oob <- round(MLmetrics::R2_Score(random_forest_model$predicted.oob, random_forest_model$yvar), 3)
  
  r_squared_list[[iteration + 1]] <- data.frame(r_squared, r_squared_oob)
  connectivity_surface_list[[iteration + 1]] <- connectivity_surface
  least_cost_lines_list[[iteration + 1]] <- least_cost_lines
  path_predictor_values_list[[iteration + 1]] <- path_predictor_values
}

oob_r_squared_index <- grep("_oob", names(unlist(r_squared_list)))

best_iteration_index <- which(
  unlist(r_squared_list) == max(unlist(r_squared_list)[oob_r_squared_index])
) / 2

best_connectivity_surface <- connectivity_surface_list[[best_iteration_index]]
best_least_cost_lines <- least_cost_lines_list[[best_iteration_index]]
best_path_predictor_values <- path_predictor_values_list[[best_iteration_index]]
best_r_squared <- r_squared_list[[best_iteration_index]]
