# Curated workflow for iterative connectivity-surface modelling.

library(terra)
library(sf)
library(dplyr)
library(purrr)
library(tidyr)
library(randomForestSRC)
library(gdistance)
library(raster)

predictor_raster_file <- "connectivity_predictors.tif"
sampling_sites_file <- "sampling_sites.gpkg"
sampling_sites_layer <- "sampling_sites"
sampling_site_id_column <- "site_id"
connectivity_matrix_file <- "neutral_connectivity_matrix.rds"
number_of_iterations <- 6
number_of_trees <- 2000

iteration_summary_file <- "connectivity_surface_iteration_summary.csv"
best_surface_file <- "best_connectivity_surface.tif"
best_paths_file <- "best_connectivity_paths.gpkg"
best_importance_file <- "best_connectivity_predictor_importance.csv"

calculate_apparent_r2 <- function(observed, predicted) {
  valid_index <- is.finite(observed) & is.finite(predicted)
  observed <- observed[valid_index]
  predicted <- predicted[valid_index]

  if (!length(observed) || stats::var(observed) == 0) {
    return(NA_real_)
  }

  1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2)
}

calculate_oob_r2 <- function(forest, observed) {
  if (is.null(forest$err.rate) || !length(forest$err.rate) || stats::var(observed) == 0) {
    return(NA_real_)
  }

  final_oob_mse <- tail(forest$err.rate, 1)
  1 - final_oob_mse / stats::var(observed)
}

load_sampling_sites <- function(file_name, layer_name, id_column, target_crs) {
  sites_sf <- st_read(file_name, layer = layer_name, quiet = TRUE) |>
    st_transform(target_crs)

  if (!id_column %in% names(sites_sf)) {
    sites_sf[[id_column]] <- paste0("site_", seq_len(nrow(sites_sf)))
  }

  sites_sf
}

ordered_site_pairs <- function(site_ids) {
  expand.grid(
    from_site = site_ids,
    to_site = site_ids,
    stringsAsFactors = FALSE
  ) |>
    filter(from_site != to_site)
}

create_initial_surface <- function(template_raster, initial_value) {
  initial_surface <- template_raster[[1]]
  initial_surface[] <- ifelse(is.na(initial_surface[]), NA, initial_value)
  initial_surface
}

surface_to_transition <- function(surface_raster) {
  conductance_raster <- raster::raster(surface_raster)
  transition_matrix <- gdistance::transition(conductance_raster, mean, directions = 8)
  gdistance::geoCorrection(transition_matrix, type = "c")
}

extract_path_covariates <- function(path_sf, predictor_stack) {
  path_values <- terra::extract(
    predictor_stack,
    terra::vect(path_sf),
    fun = mean,
    na.rm = TRUE
  )

  bind_cols(path_sf |> st_drop_geometry() |> select(from_site, to_site), path_values |> select(-ID))
}

build_directional_paths <- function(surface_raster, sites_sf, id_column) {
  transition_matrix <- surface_to_transition(surface_raster)
  site_ids <- sites_sf[[id_column]]
  site_pairs <- ordered_site_pairs(site_ids)

  path_list <- pmap(site_pairs, function(from_site, to_site) {
    from_geometry <- sites_sf |> filter(.data[[id_column]] == from_site) |> st_geometry()
    to_geometry <- sites_sf |> filter(.data[[id_column]] == to_site) |> st_geometry()

    path_sp <- gdistance::shortestPath(
      transition_matrix,
      origin = as(from_geometry, "Spatial"),
      goal = as(to_geometry, "Spatial"),
      output = "SpatialLines"
    )

    st_as_sf(path_sp) |>
      mutate(from_site = from_site, to_site = to_site)
  })

  bind_rows(path_list)
}

connectivity_matrix_to_long <- function(connectivity_matrix) {
  as.data.frame(as.table(connectivity_matrix)) |>
    rename(from_site = Var1, to_site = Var2, connectivity = Freq) |>
    filter(from_site != to_site)
}

fit_connectivity_model <- function(model_data) {
  tuning <- tune.rfsrc(connectivity ~ ., data = model_data)
  tuned_nodesize <- as.numeric(tuning$optimal["nodesize"])
  tuned_mtry <- as.numeric(tuning$optimal["mtry"])

  forest <- rfsrc(
    connectivity ~ .,
    data = model_data,
    ntree = number_of_trees,
    nodesize = tuned_nodesize,
    mtry = tuned_mtry,
    importance = TRUE
  )

  apparent_predictions <- predict(forest, newdata = model_data)$predicted

  list(
    forest = forest,
    nodesize = tuned_nodesize,
    mtry = tuned_mtry,
    apparent_r2 = calculate_apparent_r2(model_data$connectivity, apparent_predictions),
    oob_r2 = calculate_oob_r2(forest, model_data$connectivity)
  )
}

predict_connectivity_surface <- function(forest, predictor_stack) {
  predictor_table <- terra::as.data.frame(predictor_stack, na.rm = FALSE)
  complete_index <- complete.cases(predictor_table)
  predicted_values <- rep(NA_real_, nrow(predictor_table))

  predicted_values[complete_index] <- predict(
    forest,
    newdata = predictor_table[complete_index, , drop = FALSE]
  )$predicted

  predicted_surface <- predictor_stack[[1]]
  predicted_surface[] <- predicted_values
  predicted_surface
}

predictor_stack <- terra::rast(predictor_raster_file)
sampling_sites <- load_sampling_sites(
  file_name = sampling_sites_file,
  layer_name = sampling_sites_layer,
  id_column = sampling_site_id_column,
  target_crs = terra::crs(predictor_stack)
)

neutral_connectivity_matrix <- readRDS(connectivity_matrix_file)
connectivity_response <- connectivity_matrix_to_long(neutral_connectivity_matrix)
initial_connectivity_value <- mean(connectivity_response$connectivity, na.rm = TRUE)
current_surface <- create_initial_surface(predictor_stack, initial_connectivity_value)

iteration_results <- vector("list", number_of_iterations)

for (iteration_id in seq_len(number_of_iterations)) {
  message("Running connectivity-surface iteration: ", iteration_id)

  path_sf <- build_directional_paths(current_surface, sampling_sites, sampling_site_id_column)
  directional_covariates <- extract_path_covariates(path_sf, predictor_stack)

  model_data <- directional_covariates |>
    left_join(connectivity_response, by = c("from_site", "to_site")) |>
    drop_na()

  fitted_model <- fit_connectivity_model(model_data)
  predicted_surface <- predict_connectivity_surface(fitted_model$forest, predictor_stack)

  iteration_results[[iteration_id]] <- list(
    iteration = iteration_id,
    path_sf = path_sf,
    model_data = model_data,
    forest = fitted_model$forest,
    surface = predicted_surface,
    apparent_r2 = fitted_model$apparent_r2,
    oob_r2 = fitted_model$oob_r2,
    nodesize = fitted_model$nodesize,
    mtry = fitted_model$mtry
  )

  current_surface <- predicted_surface
}

iteration_summary <- map_dfr(iteration_results, function(result) {
  tibble(
    iteration = result$iteration,
    nodesize = result$nodesize,
    mtry = result$mtry,
    apparent_r2 = result$apparent_r2,
    oob_r2 = result$oob_r2
  )
})

best_iteration_id <- iteration_summary |>
  filter(oob_r2 == max(oob_r2, na.rm = TRUE)) |>
  slice(1) |>
  pull(iteration)

best_result <- iteration_results[[best_iteration_id]]

predictor_importance <- tibble(
  variable = names(best_result$forest$importance),
  importance = as.numeric(best_result$forest$importance)
) |>
  arrange(desc(importance))

write.csv(iteration_summary, iteration_summary_file, row.names = FALSE)
terra::writeRaster(best_result$surface, best_surface_file, overwrite = TRUE)
st_write(best_result$path_sf, best_paths_file, delete_dsn = TRUE, quiet = TRUE)
write.csv(predictor_importance, best_importance_file, row.names = FALSE)

message("Saved iteration summary to: ", iteration_summary_file)
message("Saved best connectivity surface to: ", best_surface_file)
message("Saved best least-cost paths to: ", best_paths_file)
message("Saved predictor importance to: ", best_importance_file)