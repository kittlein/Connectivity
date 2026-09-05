# Recruitment Random Forest for her.27.20-24

# Packages

library(readxl)
library(dplyr)
library(randomForestSRC)
library(MLmetrics)
library(Boruta)
library(ggplot2)
library(patchwork)
library(openxlsx)


# Stock and input file

stock <- "her.27.1-24a514a" 
input_file <- "herring_recruitment_ocean_1990_2026.xlsx"

# Recruitment age and environmental lag

recruitment_age <- 2L
ocean_lag <- recruitment_age + 1L

# Random Forest importance function for Boruta

getImp_rfsrc <- function(x, y, ntree = 1000, ...) {
  
  dd <- data.frame(
    response = y,
    x,
    check.names = FALSE
  )
  
  fit <- randomForestSRC::rfsrc(
    response ~ .,
    data = dd,
    ntree = ntree,
    importance = "permute",
    ...
  )
  
  imp <- fit$importance
  imp <- imp[colnames(x)]
  
  as.numeric(imp)
}


# Read stock data

stock_data <- readxl::read_excel(
  input_file,
  sheet = stock
) |>
  dplyr::mutate(
    year = as.integer(year),
    Recruitment = as.numeric(Recruitment)
  ) |>
  dplyr::arrange(year)


# Identify oceanographic predictors

ocean_cols <- setdiff(
  names(stock_data),
  c("year", "Recruitment")
)

# Recruitment time series

recruitment_data <- stock_data |>
  dplyr::select(
    year,
    Recruitment
  )


# Oceanographic time series

ocean_data <- stock_data |>
  dplyr::select(
    year,
    dplyr::all_of(ocean_cols)
  ) |>
  dplyr::rename(
    ocean_year = year
  ) |>
  dplyr::mutate(
    year = ocean_year + ocean_lag
  )


# Match recruitment with lagged oceanographic conditions

analysis_all <- recruitment_data |>
  dplyr::left_join(
    ocean_data,
    by = "year"
  ) |>
  dplyr::arrange(year)


# Keep rows with recruitment and complete oceanographic data

analysis_data <- analysis_all |>
  dplyr::filter(
    !is.na(Recruitment)
  ) |>
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(ocean_cols),
      ~ is.finite(.x)
    )
  )


# Prepare Random Forest data

dat <- analysis_data |>
  dplyr::select(
    -year,
    -ocean_year
  ) |>
  data.frame(
    check.names = FALSE
  )


# Boruta variable selection

set.seed(123)

bor <- Boruta::Boruta(
  Recruitment ~ .,
  data = dat,
  getImp = getImp_rfsrc,
  ntree = 1000,
  maxRuns = 200,
  doTrace = 0
)


# Boruta diagnostics

bor_stats <- Boruta::attStats(
  bor
)

bor_stats <- bor_stats |>
  tibble::rownames_to_column(
    "variable"
  ) |>
  dplyr::arrange(
    dplyr::desc(meanImp)
  )

cat("\nBoruta decisions:\n")
print(table(bor_stats$decision))

print(bor_stats)


# Confirmed variables

selected_vars <- Boruta::getSelectedAttributes(
  bor,
  withTentative = FALSE
)

cat(
  "\nConfirmed variables:\n",
  paste(selected_vars, collapse = "\n"),
  "\n"
)


# Stop if Boruta finds no confirmed predictor

if (length(selected_vars) == 0L) {
  
  tentative_vars <- Boruta::getSelectedAttributes(
    bor,
    withTentative = TRUE
  )
  
  cat(
    "\nConfirmed + tentative variables:\n",
    paste(tentative_vars, collapse = "\n"),
    "\n"
  )
  
  stop(
    "Boruta did not select any confirmed variables. ",
    "Inspect bor_stats before proceeding."
  )
}


# Final Random Forest data

dat_rf <- dat |>
  dplyr::select(
    Recruitment,
    dplyr::all_of(selected_vars)
  ) |>
  data.frame(
    check.names = FALSE
  )

rf_formula <- reformulate(
  selected_vars,
  response = "Recruitment"
)


# Fit Random Forest

set.seed(123)

rf <- randomForestSRC::rfsrc(
  formula = rf_formula,
  data = dat_rf,
  ntree = 1000,
  importance = "permute"
)


# Model performance

r2_fit <- MLmetrics::R2_Score(
  rf$predicted,
  rf$yvar
)

r2_oob <- MLmetrics::R2_Score(
  rf$predicted.oob,
  rf$yvar
)

cat(
  "\nR2 =", round(r2_fit, 3),
  "\nR2 OOB =", round(r2_oob, 3),
  "\n"
)


# Variable labels

var_labels <- c(
  "chla_min_minFst" = "Min chlorophyll-a min FST",
  "chla_max_minFst" = "Max chlorophyll-a min FST",
  "sst_min_minFst" = "Min temperature min FST",
  "sst_max_minFst" = "Max temperature min FST",
  "sss_min_minFst" = "Min salinity min FST",
  "sss_max_minFst" = "Max salinity min FST",
  "uo_min_minFst" = "Min zonal velocity min FST",
  "uo_max_minFst" = "Max zonal velocity min FST",
  "vo_min_minFst" = "Min meridional velocity min FST",
  "vo_max_minFst" = "Max meridional velocity min FST",
  "chla_min_maxFst" = "Min chlorophyll-a max FST",
  "chla_max_maxFst" = "Max chlorophyll-a max FST",
  "sst_min_maxFst" = "Min temperature max FST",
  "sst_max_maxFst" = "Max temperature max FST",
  "sss_min_maxFst" = "Min salinity max FST",
  "sss_max_maxFst" = "Max salinity max FST",
  "uo_min_maxFst" = "Min zonal velocity max FST",
  "uo_max_maxFst" = "Max zonal velocity max FST",
  "vo_min_maxFst" = "Min meridional velocity max FST",
  "vo_max_maxFst" = "Max meridional velocity max FST"
)


# Permutation importance

vi <- randomForestSRC::vimp(
  rf,
  importance = "permute"
)$importance

importance_table <- data.frame(
  variable = names(vi),
  importance = as.numeric(vi),
  check.names = FALSE
)

max_importance <- max(
  importance_table$importance,
  na.rm = TRUE
)

importance_table$relative_importance <-
  importance_table$importance / max_importance

importance_table$label <- unname(
  var_labels[importance_table$variable]
)

missing_labels <- is.na(
  importance_table$label
)

importance_table$label[missing_labels] <-
  importance_table$variable[missing_labels]

importance_table <- importance_table |>
  dplyr::arrange(
    dplyr::desc(relative_importance)
  )


# Select four most important predictors

top_vars <- head(
  importance_table$variable,
  min(4L, nrow(importance_table))
)

cat(
  "\nTop variables:\n",
  paste(top_vars, collapse = "\n"),
  "\n"
)


# Partial dependence

partial_list <- lapply(
  top_vars,
  function(v) {
    
    x_values <- sort(
      unique(
        dat_rf[[v]]
      )
    )
    
    po <- randomForestSRC::partial(
      rf,
      partial.xvar = v,
      partial.values = x_values,
      oob = TRUE
    )
    
    pd <- randomForestSRC::get.partial.plot.data(
      po
    )
    
    data.frame(
      variable = v,
      x = pd$x,
      partial_dependence = pd$yhat
    )
  }
)

partial_table <- dplyr::bind_rows(
  partial_list
)


# Recruitment scale for partial plots

ymax <- max(
  analysis_data$Recruitment,
  na.rm = TRUE
)

if (ymax >= 5e6) {
  
  div_y <- 1e6
  exp_y <- 6
  
} else {
  
  div_y <- 1e5
  exp_y <- 5
}


# Prediction table

prediction_table <- data.frame(
  year = analysis_data$year,
  ocean_year = analysis_data$ocean_year,
  observed = rf$yvar,
  fitted = rf$predicted,
  oob = rf$predicted.oob
)


# Recruitment fit plot

fit_plot <- ggplot2::ggplot(
  prediction_table,
  ggplot2::aes(x = year)
) +
  ggplot2::geom_line(
    ggplot2::aes(
      y = observed,
      colour = "Observed"
    ),
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      y = observed,
      colour = "Observed"
    ),
    size = 1.8
  ) +
  ggplot2::geom_line(
    ggplot2::aes(
      y = fitted,
      colour = "Fitted"
    ),
    linewidth = 0.9
  ) +
  ggplot2::geom_line(
    ggplot2::aes(
      y = oob,
      colour = "OOB"
    ),
    linewidth = 0.9
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Observed" = "black",
      "Fitted" = "red",
      "OOB" = "darkgreen"
    )
  ) +
  ggplot2::labs(
    title = stock,
    subtitle = bquote(
      R^2 == .(sprintf("%.2f", r2_fit)) ~
        "   " ~
        R[OOB]^2 == .(sprintf("%.2f", r2_oob)) ~
        "   lag =" ~ .(ocean_lag)
    ),
    x = "Year",
    y = "Recruitment",
    colour = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 12
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    legend.position = "top"
  )


# Variable importance plot

importance_plot_data <- importance_table |>
  dplyr::mutate(
    label = factor(
      label,
      levels = rev(label)
    )
  )

importance_plot <- ggplot2::ggplot(
  importance_plot_data,
  ggplot2::aes(
    x = label,
    y = relative_importance
  )
) +
  ggplot2::geom_col(
    width = 0.7
  ) +
  ggplot2::coord_flip() +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  ggplot2::labs(
    title = stock,
    subtitle = bquote(
      R^2 == .(sprintf("%.2f", r2_fit)) ~
        "   " ~
        R[OOB]^2 == .(sprintf("%.2f", r2_oob))
    ),
    x = NULL,
    y = "Relative variable importance"
  ) +
  ggplot2::theme_classic(
    base_size = 11
  )


# Partial-dependence plots

partial_plots <- lapply(
  top_vars,
  function(v) {
    
    pd <- partial_table |>
      dplyr::filter(
        variable == v
      )
    
    variable_label <- unname(
      var_labels[v]
    )
    
    if (is.na(variable_label)) {
      variable_label <- v
    }
    
    ggplot2::ggplot(
      pd,
      ggplot2::aes(
        x = x,
        y = partial_dependence / div_y
      )
    ) +
      ggplot2::geom_line(
        linewidth = 0.9
      ) +
      ggplot2::labs(
        title = variable_label,
        x = NULL,
        y = bquote(PD ~ (x ~ 10^.(exp_y)))
      ) +
      ggplot2::theme_classic(
        base_size = 11
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 10,
          hjust = 0
        )
      )
  }
)


# Combined effects figure

effects_figure <- importance_plot /
  patchwork::wrap_plots(
    partial_plots,
    ncol = 2
  )

print(fit_plot)
print(effects_figure)


# Model summary

model_summary <- data.frame(
  stock = stock,
  recruitment_age = recruitment_age,
  ocean_lag = ocean_lag,
  n = nrow(analysis_data),
  recruitment_year_min = min(analysis_data$year),
  recruitment_year_max = max(analysis_data$year),
  ocean_year_min = min(analysis_data$ocean_year),
  ocean_year_max = max(analysis_data$ocean_year),
  R2 = r2_fit,
  R2_OOB = r2_oob
)
