# Curated workflow for stock-specific recruitment models and variable importance.

library(randomForestSRC)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

input_workbook <- "Recruitment data Herring stocks chla sst sss vo uo.xlsx"
recruitment_sheet <- "Recruitment"
catch_sheet <- "Catches"
environment_sheet <- "Environment"

recruitment_years <- 1993:2022
catch_years <- 1992:2021
number_of_trees <- 2000
number_of_importance_replicates <- 200

calculate_apparent_r2 <- function(observed, predicted) {
  valid_index <- is.finite(observed) & is.finite(predicted)
  observed <- observed[valid_index]
  predicted <- predicted[valid_index]

  if (!length(observed) || stats::var(observed) == 0) {
    return(NA_real_)
  }

  1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2)
}

standardize_environment_names <- function(environment_data) {
  names(environment_data) <- names(environment_data) |>
    str_replace("pol", "p") |>
    str_replace("spring", "S") |>
    str_replace("fall", "F")

  environment_data
}

build_stock_dataset <- function(stock_name, recruitment_data, catch_data, environment_data) {
  recruitment_index <- recruitment_data$Years %in% recruitment_years
  catch_index <- catch_data$Years %in% catch_years

  stock_data <- data.frame(
    Year = recruitment_data$Years[recruitment_index],
    Recruitment = recruitment_data[recruitment_index, stock_name][[1]],
    LaggedCatch = catch_data[catch_index, stock_name][[1]],
    environment_data[seq_len(sum(recruitment_index)), -1],
    check.names = FALSE
  )

  stock_data
}

fit_random_forest_model <- function(stock_data) {
  model_data <- stock_data |>
    dplyr::select(-Year)

  tuning <- tune.rfsrc(
    Recruitment ~ ., 
    data = model_data
  )

  tuned_nodesize <- as.numeric(tuning$optimal["nodesize"])
  tuned_mtry <- as.numeric(tuning$optimal["mtry"])

  forest <- rfsrc(
    Recruitment ~ ., 
    data = model_data,
    ntree = number_of_trees,
    nodesize = tuned_nodesize,
    mtry = tuned_mtry
  )

  list(
    forest = forest,
    nodesize = tuned_nodesize,
    mtry = tuned_mtry,
    apparent_r2 = calculate_apparent_r2(model_data$Recruitment, forest$predicted)
  )
}

summarize_variable_importance <- function(stock_data, tuned_nodesize, tuned_mtry, n_replicates) {
  model_data <- stock_data |>
    dplyr::select(-Year)

  importance_replicates <- map_dfr(seq_len(n_replicates), function(replicate_id) {
    forest <- rfsrc(
      Recruitment ~ ., 
      data = model_data,
      ntree = number_of_trees,
      nodesize = tuned_nodesize,
      mtry = tuned_mtry,
      importance = "permute",
      block.size = 1
    )

    variable_importance <- vimp(forest)$importance

    tibble(
      variable = names(variable_importance),
      importance = as.numeric(variable_importance),
      replicate = replicate_id
    )
  })

  importance_replicates |>
    group_by(variable) |>
    summarise(
      mean_importance = mean(importance, na.rm = TRUE),
      sd_importance = sd(importance, na.rm = TRUE),
      q025 = quantile(importance, 0.025, na.rm = TRUE),
      q50 = quantile(importance, 0.50, na.rm = TRUE),
      q975 = quantile(importance, 0.975, na.rm = TRUE),
      prob_positive = mean(importance > 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      mean_importance_positive = pmax(mean_importance, 0),
      relative_importance = if_else(
        sum(mean_importance_positive, na.rm = TRUE) > 0,
        100 * mean_importance_positive / sum(mean_importance_positive, na.rm = TRUE),
        0
      )
    ) |>
    arrange(desc(relative_importance))
}

recruitment_data <- readxl::read_excel(input_workbook, sheet = recruitment_sheet, skip = 1)
catch_data <- readxl::read_excel(input_workbook, sheet = catch_sheet)
environment_data <- readxl::read_excel(input_workbook, sheet = environment_sheet) |>
  standardize_environment_names()

stock_names <- names(recruitment_data)[-1]

model_results <- map(stock_names, function(stock_name) {
  stock_data <- build_stock_dataset(stock_name, recruitment_data, catch_data, environment_data)
  fitted_model <- fit_random_forest_model(stock_data)
  importance_summary <- summarize_variable_importance(
    stock_data = stock_data,
    tuned_nodesize = fitted_model$nodesize,
    tuned_mtry = fitted_model$mtry,
    n_replicates = number_of_importance_replicates
  )

  list(
    stock_name = stock_name,
    stock_data = stock_data,
    fitted_model = fitted_model,
    importance_summary = importance_summary
  )
})

fit_summary <- map_dfr(model_results, function(result) {
  tibble(
    stock = result$stock_name,
    nodesize = result$fitted_model$nodesize,
    mtry = result$fitted_model$mtry,
    apparent_r2 = result$fitted_model$apparent_r2
  )
})

fitted_series <- map_dfr(model_results, function(result) {
  tibble(
    stock = result$stock_name,
    year = result$stock_data$Year,
    observed_recruitment = result$stock_data$Recruitment,
    predicted_recruitment = result$fitted_model$forest$predicted
  )
})

importance_summary_long <- map_dfr(model_results, function(result) {
  result$importance_summary |>
    mutate(stock = result$stock_name, .before = 1)
})

importance_summary_wide <- importance_summary_long |>
  dplyr::select(stock, variable, relative_importance) |>
  pivot_wider(
    names_from = stock,
    values_from = relative_importance
  )

write.csv(fit_summary, "recruitment_model_fit_summary.csv", row.names = FALSE)
write.csv(fitted_series, "recruitment_model_fitted_series.csv", row.names = FALSE)
write.csv(importance_summary_long, "recruitment_variable_importance_long.csv", row.names = FALSE)
write.csv(importance_summary_wide, "recruitment_variable_importance_wide.csv", row.names = FALSE)

message("Saved model summary to: recruitment_model_fit_summary.csv")
message("Saved fitted time series to: recruitment_model_fitted_series.csv")
message("Saved variable importance summary to: recruitment_variable_importance_long.csv")
message("Saved relative importance table to: recruitment_variable_importance_wide.csv")