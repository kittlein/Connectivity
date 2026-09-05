################################################################################
# ENVIRONMENTAL WINDOWS AT EITHER MIN-FST OR MAX-FST LOCATION
# CHLA DELAY + BORUTA SELECTION + RF TUNING + RANDOM FOREST
################################################################################

library(terra); library(dplyr); library(lubridate);
library(readxl); library(randomForestSRC); library(Boruta)
library(ale)

################################################################################
# USER SETTINGS
################################################################################

fst.use <- "minFST"       # OPTIONS: "minFST" or "maxFST"

s.age <- c(2, 0, 0, 0, 0, 0, 2, 1, 1)
lags <- 0:36; 
windows <- c(1, 2, 3, 4); 
chla.delays <- 0:3; 
reference.month <- 1

use.log.recruitment <- TRUE; ntree.rf <- 5000; ntree.tune <- 5000; maxRuns.boruta <- 2000; seed.rf <- 123

if (!fst.use %in% c("minFST", "maxFST")) stop("fst.use must be either 'minFST' or 'maxFST'.")


################################################################################
# READ MONTHLY OCEANOGRAPHIC VARIABLES
################################################################################

chla <- terra::rast("monthly_ocean_layers/chla_monthly_1993-2026-herring.nc")
uo <- terra::rast("monthly_ocean_layers/uo_monthly_1993-2026-herring.nc")
vo <- terra::rast("monthly_ocean_layers/vo_monthly_1993-2026-herring.nc")
sss <- terra::rast("monthly_ocean_layers/sss_monthly_1993-2026-herring.nc")
sst <- terra::rast("monthly_ocean_layers/sst_monthly_1993-2026-herring.nc")


################################################################################
# READ MINIMUM-FST AND MAXIMUM-FST LOCATIONS
################################################################################

coords_file <- "herring_stock_FST_extreme_locations.xlsx"

min_fst_points <- readxl::read_excel(coords_file, sheet = "min_FST") |> dplyr::transmute(stock = as.character(stock), fst_location = "minFST", predicted_fst = as.numeric(fst), longitude = as.numeric(longitude), latitude = as.numeric(latitude))

max_fst_points <- readxl::read_excel(coords_file, sheet = "max_FST") |> dplyr::transmute(stock = as.character(stock), fst_location = "maxFST", predicted_fst = as.numeric(fst), longitude = as.numeric(longitude), latitude = as.numeric(latitude))

for(fila in 1:9){             # Stock is defined according to the row in min_fst_points
################################################################################
# SELECT STOCK AND FST LOCATION
################################################################################

stock.name <- min_fst_points$stock[fila]

points.use <- if (fst.use == "minFST") min_fst_points else max_fst_points

fila.use <- match(stock.name, points.use$stock)

if (is.na(fila.use)) stop(paste("Stock", stock.name, "was not found in the selected", fst.use, "table."))

prefix <- if (fst.use == "minFST") "minFst" else "maxFst"

cat("\nStock:", stock.name, "\n"); cat("FST location used:", fst.use, "\n")
cat("Predicted FST:", points.use$predicted_fst[fila.use], "\n")
cat("Coordinates:", points.use$longitude[fila.use], points.use$latitude[fila.use], "\n\n")


################################################################################
# READ RECRUITMENT
################################################################################

Reclu <- readxl::read_excel("herring_recruitment_1990_2026.xlsx", sheet = stock.name, skip = 8)

Reclu <- data.frame(Year = as.integer(Reclu$Year), Recruitment = as.numeric(Reclu$Recruitment))

if (use.log.recruitment && any(Reclu$Recruitment <= 0, na.rm = TRUE)) stop("Recruitment contains values <= 0 and therefore cannot be log-transformed.")


################################################################################
# CHECK MONTHLY TIME DIMENSION
################################################################################

fechas <- as.Date(terra::time(chla))

stopifnot(length(fechas) == terra::nlyr(chla), length(fechas) == terra::nlyr(uo), length(fechas) == terra::nlyr(vo), length(fechas) == terra::nlyr(sss), length(fechas) == terra::nlyr(sst))


################################################################################
# EXTRACT OCEANOGRAPHIC SERIES AT THE SELECTED FST LOCATION
################################################################################

chla_vals <- as.numeric(unlist(terra::extract(chla, points.use[fila.use, c("longitude", "latitude")], method = "simple", ID = FALSE), use.names = FALSE))

uo_vals <- as.numeric(unlist(terra::extract(uo, points.use[fila.use, c("longitude", "latitude")], method = "simple", ID = FALSE), use.names = FALSE))

vo_vals <- as.numeric(unlist(terra::extract(vo, points.use[fila.use, c("longitude", "latitude")], method = "simple", ID = FALSE), use.names = FALSE))

sss_vals <- as.numeric(unlist(terra::extract(sss, points.use[fila.use, c("longitude", "latitude")], method = "simple", ID = FALSE), use.names = FALSE))

sst_vals <- as.numeric(unlist(terra::extract(sst, points.use[fila.use, c("longitude", "latitude")], method = "simple", ID = FALSE), use.names = FALSE))

stopifnot(length(fechas) == length(chla_vals), length(fechas) == length(uo_vals), length(fechas) == length(vo_vals), length(fechas) == length(sss_vals), length(fechas) == length(sst_vals))


################################################################################
# CREATE MONTHLY OCEANOGRAPHIC DATA FRAME
################################################################################

ocean <- data.frame(fecha = fechas, chla = chla_vals, uo = uo_vals, vo = vo_vals, sss = sss_vals, sst = sst_vals)

names(ocean)[-1] <- paste0(prefix, ".", c("chla", "uo", "vo", "sss", "sst"))

env.vars <- paste0(prefix, ".", c("chla", "uo", "vo", "sss", "sst"))

chla.vars <- paste0(prefix, ".chla")

physical.vars <- paste0(prefix, ".", c("uo", "vo", "sss", "sst"))

ocean$fecha <- as.Date(ocean$fecha); ocean <- ocean[order(ocean$fecha), c("fecha", env.vars)]; rownames(ocean) <- NULL

cat("Oceanographic predictors:\n"); print(env.vars); cat("\n")


################################################################################
# PARAMETER GRID
#
# lag       = lag of physical variables relative to January of recruitment year
# window    = number of consecutive months averaged
# chla_delay = number of months Chl-a occurs later than the physical variables
#
# Effective Chl-a lag = lag - chla_delay
################################################################################

parameter.grid <- expand.grid(lag = lags+12*s.age[fila], window = windows, chla_delay = chla.delays, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

parameter.grid <- parameter.grid[parameter.grid$chla_delay <= parameter.grid$lag, , drop = FALSE]; rownames(parameter.grid) <- NULL

stopifnot(all(c("lag", "window", "chla_delay") %in% names(parameter.grid)))

cat("Number of lag x window x Chl-a delay combinations:", nrow(parameter.grid), "\n\n")


################################################################################
# CONSTRUCT ENVIRONMENTAL WINDOWS
################################################################################

make_window_data <- function(lag.months, window.months, chla.delay, reference.month, Reclu, ocean, env.vars) {
  
  if (length(lag.months) != 1L || is.na(lag.months)) stop("lag.months must contain one valid value.")
  if (length(window.months) != 1L || is.na(window.months)) stop("window.months must contain one valid value.")
  if (length(chla.delay) != 1L || is.na(chla.delay)) stop("chla.delay must contain one valid value.")
  
  chla.vars <- env.vars[grepl("\\.chla$", env.vars)]; physical.vars <- env.vars[!grepl("\\.chla$", env.vars)]
  
  ans <- lapply(seq_len(nrow(Reclu)), function(i) {
    
    reference.date <- as.Date(sprintf("%04d-%02d-01", as.integer(Reclu$Year[i]), reference.month))
    
    physical.end <- lubridate::`%m-%`(reference.date, lubridate::period(month = lag.months))
    
    physical.start <- lubridate::`%m-%`(physical.end, lubridate::period(month = window.months - 1))
    
    chla.lag <- lag.months - chla.delay
    
    chla.end <- lubridate::`%m-%`(reference.date, lubridate::period(month = chla.lag))
    
    chla.start <- lubridate::`%m-%`(chla.end, lubridate::period(month = window.months - 1))
    
    z.chla <- ocean[ocean$fecha >= chla.start & ocean$fecha <= chla.end, c("fecha", chla.vars), drop = FALSE]
    
    complete.chla <- nrow(z.chla) == window.months && length(unique(format(z.chla$fecha, "%Y-%m"))) == window.months && all(complete.cases(z.chla[, chla.vars, drop = FALSE]))
    
    chla.means <- if (complete.chla) vapply(z.chla[, chla.vars, drop = FALSE], function(x) mean(x, na.rm = TRUE), numeric(1)) else setNames(rep(NA_real_, length(chla.vars)), chla.vars)
    
    z.physical <- ocean[ocean$fecha >= physical.start & ocean$fecha <= physical.end, c("fecha", physical.vars), drop = FALSE]
    
    complete.physical <- nrow(z.physical) == window.months && length(unique(format(z.physical$fecha, "%Y-%m"))) == window.months && all(complete.cases(z.physical[, physical.vars, drop = FALSE]))
    
    physical.means <- if (complete.physical) vapply(z.physical[, physical.vars, drop = FALSE], function(x) mean(x, na.rm = TRUE), numeric(1)) else setNames(rep(NA_real_, length(physical.vars)), physical.vars)
    
    env.means <- c(chla.means, physical.means)
    
    data.frame(Year = Reclu$Year[i], Recruitment = Reclu$Recruitment[i], lag = lag.months, window = window.months, chla_delay = chla.delay, physical_start = physical.start, physical_end = physical.end, chla_start = chla.start, chla_end = chla.end, t(env.means), check.names = FALSE)
  })
  
  dplyr::bind_rows(ans)
}


################################################################################
# GENERATE ALL LAG x WINDOW x CHLA-DELAY DATASETS
################################################################################

lag.window.data <- lapply(seq_len(nrow(parameter.grid)), function(i) make_window_data(lag.months = parameter.grid$lag[i], window.months = parameter.grid$window[i], chla.delay = parameter.grid$chla_delay[i], reference.month = reference.month, Reclu = Reclu, ocean = ocean, env.vars = env.vars))


################################################################################
# OPTIONAL CHECK OF TEMPORAL WINDOWS
################################################################################

test.index <- which(parameter.grid$lag == 12 & parameter.grid$window == 3 & parameter.grid$chla_delay == 2)

if (length(test.index) == 1L) print(head(lag.window.data[[test.index]][, c("Year", "physical_start", "physical_end", "chla_start", "chla_end", env.vars)], 10))


################################################################################
# USE EXACTLY THE SAME RECRUITMENT YEARS IN ALL MODELS
################################################################################

valid.years <- Reduce(intersect, lapply(lag.window.data, function(z) z$Year[complete.cases(z[, c("Recruitment", env.vars), drop = FALSE])]))

if (length(valid.years) == 0L) stop("No common recruitment years remain across all parameter combinations.")

cat("\nCommon years used in all models:", min(valid.years), "-", max(valid.years), "\n")
cat("Number of recruitment years:", length(valid.years), "\n")
cat("Number of environmental predictors entering Boruta:", length(env.vars), "\n\n")


################################################################################
# BORUTA + TUNING OF MTRY AND NODESIZE + RANDOM FOREST
################################################################################

evaluate_model_rf <- function(i, lag.window.data, parameter.grid, valid.years, env.vars, use.log.recruitment = TRUE, ntree = 2000, ntree.tune = 500, maxRuns = 200, seed = 123) {
  
  dat <- lag.window.data[[i]]; dat <- dat[dat$Year %in% valid.years, , drop = FALSE]
  
  dat$Y <- if (use.log.recruitment) log(dat$Recruitment) else dat$Recruitment
  
  dat.mod <- dat[, c("Y", env.vars), drop = FALSE]; dat.mod <- dat.mod[complete.cases(dat.mod), , drop = FALSE]
  
  if (nrow(dat.mod) < 10L) return(data.frame(lag = parameter.grid$lag[i], window = parameter.grid$window[i], chla_delay = parameter.grid$chla_delay[i], n = nrow(dat.mod), nvar = NA_integer_, selected = NA_character_, mtry = NA_integer_, nodesize = NA_integer_, R2fit = NA_real_, R2OOB = NA_real_, RMSEfit = NA_real_, RMSEOOB = NA_real_))
  
  set.seed(seed + i); boruta.form <- stats::reformulate(env.vars, response = "Y")
  
  boruta.fit <- Boruta::Boruta(formula = boruta.form, data = dat.mod, maxRuns = maxRuns, doTrace = 0)
  
  boruta.fixed <- Boruta::TentativeRoughFix(boruta.fit)
  
  selected.vars <- Boruta::getSelectedAttributes(boruta.fixed, withTentative = FALSE)
  
  if (length(selected.vars) == 0L) return(data.frame(lag = parameter.grid$lag[i], window = parameter.grid$window[i], chla_delay = parameter.grid$chla_delay[i], n = nrow(dat.mod), nvar = 0L, selected = "None", mtry = NA_integer_, nodesize = NA_integer_, R2fit = NA_real_, R2OOB = NA_real_, RMSEfit = NA_real_, RMSEOOB = NA_real_))
  
  rf.form <- stats::reformulate(selected.vars, response = "Y"); rf.dat <- dat.mod[, c("Y", selected.vars), drop = FALSE]
  
  set.seed(seed + i)
  
  tune.fit <- randomForestSRC::tune(formula = rf.form, data = rf.dat, nodesize.try = 1:10, ntree.try = ntree.tune, nsplit = 10, method = "grid", do.best = FALSE, trace = FALSE, seed = seed + i)
  
  best.nodesize <- as.integer(tune.fit$optimal["nodesize"]); best.mtry <- as.integer(tune.fit$optimal["mtry"])
  
  if (!is.finite(best.nodesize) || !is.finite(best.mtry)) stop(paste("Tuning failed for parameter-grid row", i))
  
  set.seed(seed + i)
  
  rf.fit <- randomForestSRC::rfsrc(formula = rf.form, data = rf.dat, ntree = ntree, mtry = best.mtry, nodesize = best.nodesize, nsplit = 10, importance = "permute", forest = TRUE)
  
  obs <- rf.dat$Y; pred.fit <- as.numeric(rf.fit$predicted); pred.oob <- as.numeric(rf.fit$predicted.oob)
  
  SST <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  
  R2fit <- 1 - sum((obs - pred.fit)^2, na.rm = TRUE) / SST; R2OOB <- 1 - sum((obs - pred.oob)^2, na.rm = TRUE) / SST
  
  RMSEfit <- sqrt(mean((obs - pred.fit)^2, na.rm = TRUE)); RMSEOOB <- sqrt(mean((obs - pred.oob)^2, na.rm = TRUE))
  
  data.frame(lag = parameter.grid$lag[i], window = parameter.grid$window[i], chla_delay = parameter.grid$chla_delay[i], n = nrow(rf.dat), nvar = length(selected.vars), selected = paste(selected.vars, collapse = ", "), mtry = best.mtry, nodesize = best.nodesize, R2fit = R2fit, R2OOB = R2OOB, RMSEfit = RMSEfit, RMSEOOB = RMSEOOB)
}


################################################################################
# RUN ALL RANDOM FOREST MODELS
################################################################################

rf.results <- dplyr::bind_rows(lapply(seq_len(nrow(parameter.grid)), function(i) evaluate_model_rf(i = i, lag.window.data = lag.window.data, parameter.grid = parameter.grid, valid.years = valid.years, env.vars = env.vars, use.log.recruitment = use.log.recruitment, ntree = ntree.rf, ntree.tune = ntree.tune, maxRuns = maxRuns.boruta, seed = seed.rf)))

rf.results <- rf.results[order(-rf.results$R2OOB), ]; rownames(rf.results) <- NULL

print(head(rf.results, 30))


################################################################################
# SELECT BEST MODEL ACCORDING TO OOB R2
################################################################################

if (!any(is.finite(rf.results$R2OOB))) stop("No model produced a finite OOB R2.")

best.rf <- rf.results[which.max(rf.results$R2OOB), ]

best.index <- which(parameter.grid$lag == best.rf$lag & parameter.grid$window == best.rf$window & parameter.grid$chla_delay == best.rf$chla_delay)

stopifnot(length(best.index) == 1L)

dat.best <- lag.window.data[[best.index]]; dat.best <- dat.best[dat.best$Year %in% valid.years, , drop = FALSE]

selected.vars <- trimws(strsplit(best.rf$selected, ",")[[1]])

dat.best$Y <- if (use.log.recruitment) log(dat.best$Recruitment) else dat.best$Recruitment

rf.dat <- dat.best[, c("Y", selected.vars), drop = FALSE]; ok <- complete.cases(rf.dat)

rf.dat <- rf.dat[ok, , drop = FALSE]; dat.best <- dat.best[ok, , drop = FALSE]

rf.form <- stats::reformulate(selected.vars, response = "Y")

best.mtry <- as.integer(best.rf$mtry); best.nodesize <- as.integer(best.rf$nodesize)


################################################################################
# RECONSTRUCT BEST MODEL WITH EXACTLY THE TUNED PARAMETERS
################################################################################

set.seed(seed.rf + best.index)

best.fit <- randomForestSRC::rfsrc(formula = rf.form, data = rf.dat, ntree = ntree.rf, mtry = best.mtry, nodesize = best.nodesize, nsplit = 10, importance = "permute", forest = TRUE)

pred.fit <- as.numeric(best.fit$predicted); pred.oob <- as.numeric(best.fit$predicted.oob)

SST <- sum((rf.dat$Y - mean(rf.dat$Y, na.rm = TRUE))^2, na.rm = TRUE)

R2.fit <- 1 - sum((rf.dat$Y - pred.fit)^2, na.rm = TRUE) / SST; R2.oob <- 1 - sum((rf.dat$Y - pred.oob)^2, na.rm = TRUE) / SST

RMSE.fit <- sqrt(mean((rf.dat$Y - pred.fit)^2, na.rm = TRUE)); RMSE.oob <- sqrt(mean((rf.dat$Y - pred.oob)^2, na.rm = TRUE))


################################################################################
# REPORT BEST MODEL
################################################################################

cat("\nBest model\n")
cat("Stock:", stock.name, "\n")
cat("FST location:", fst.use, "\n")
cat("Physical lag:", best.rf$lag, "months\n")
cat("Window:", best.rf$window, "months\n")
cat("Chl-a delay:", best.rf$chla_delay, "months\n")
cat("Effective Chl-a lag:", best.rf$lag - best.rf$chla_delay, "months\n")
cat("Selected variables:", paste(selected.vars, collapse = ", "), "\n")
cat("mtry:", best.mtry, "\n")
cat("nodesize:", best.nodesize, "\n")
cat("R2 fit:", round(R2.fit, 4), "\n")
cat("R2 OOB from search:", round(best.rf$R2OOB, 4), "\n")
cat("R2 OOB reconstructed:", round(R2.oob, 4), "\n")
cat("RMSE fit:", round(RMSE.fit, 4), "\n")
cat("RMSE OOB:", round(RMSE.oob, 4), "\n\n")


################################################################################
# BACK-TRANSFORM PREDICTIONS FOR PLOTTING
################################################################################

if (use.log.recruitment) {
  
  observed <- exp(rf.dat$Y); predicted <- exp(pred.fit); predicted.oob <- exp(pred.oob)
  
} else {
  
  observed <- rf.dat$Y; predicted <- pred.fit; predicted.oob <- pred.oob
}


################################################################################
# OBSERVED, FITTED AND OOB RECRUITMENT SERIES
################################################################################

ylim <- range(c(observed, predicted, predicted.oob), na.rm = TRUE)

plot.name <- paste0(stock.name, "  ", fst.use)

main.txt <- bquote(atop(bold(.(plot.name)), R[fit]^2 == .(sprintf("%.2f", R2.fit)) ~~ R[OOB]^2 == .(sprintf("%.2f", R2.oob)) ~~ lag == .(best.rf$lag) ~~ window == .(best.rf$window) ~~ chla.delay == .(best.rf$chla_delay)))

par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1))

plot(dat.best$Year, observed, type = "b", pch = 19, lwd = 2, ylim = ylim, las = 1, xlab = "Year", ylab = "Recruitment", main = main.txt, cex.main = 1.1)

lines(dat.best$Year, predicted, type = "l", lwd = 2, col = "blue")

lines(dat.best$Year, predicted.oob, type = "l", lwd = 2, col = "red")

legend("topright", legend = c("Observed", "Predicted", "OOB predicted"), col = c("black", "blue", "red"), pch = c(19, NA, NA), lty = 1, lwd = 2, bty = "n")


################################################################################
# PARTIAL RESPONSE PLOTS FOR VARIABLES SELECTED BY BORUTA
################################################################################

nvar <- length(selected.vars); ncol.plot <- min(2, nvar); nrow.plot <- ceiling(nvar / ncol.plot)

par(mfrow = c(nrow.plot, ncol.plot), mar = c(4.5, 4.5, 2, 1))

for (xvar in selected.vars) {
  
  x.grid <- seq(min(rf.dat[[xvar]], na.rm = TRUE), max(rf.dat[[xvar]], na.rm = TRUE), length.out = 50)
  
  pd.obj <- randomForestSRC::partial(best.fit, partial.xvar = xvar, partial.values = x.grid)
  
  pd <- randomForestSRC::get.partial.plot.data(pd.obj)
  
  partial.response <- if (use.log.recruitment) exp(pd$yhat) else pd$yhat
  
  plot(pd$x, partial.response, type = "l", lwd = 2, las = 1, xlab = xvar, ylab = "Partial recruitment response")
  
  points(pd$x, partial.response, pch = 19, cex = 0.5)
}

par(mfrow = c(1, 1))


################################################################################
# OPTIONAL: SAVE COMPLETE MODEL SEARCH
################################################################################

# write.csv(rf.results, paste0("RF_", fst.use, "_lag_window_chlaDelay_", stock.name, ".csv"), row.names = FALSE)

################################################################################
# SPATIAL RANDOMIZATION TEST
#
# Compare the OOB R2 obtained at the selected minFST/maxFST location against
# OOB R2 values obtained from random marine cells within the stock polygon.
#
# IMPORTANT:
#   Boruta selection is NOT repeated.
#   Tuning is NOT repeated.
#   lag, window, chla_delay, selected variables, mtry and nodesize are fixed.
#   Recruitment years are also kept fixed.
################################################################################


################################################################################
# SETTINGS
################################################################################

n.random <- 1000

seed.spatial <- 987

focal.R2OOB <- R2.oob

best.lag <- as.integer(best.rf$lag)

best.window <- as.integer(best.rf$window)

best.chla.delay <- as.integer(best.rf$chla_delay)

best.mtry <- as.integer(best.rf$mtry)

best.nodesize <- as.integer(best.rf$nodesize)

selected.vars.fixed <- selected.vars


################################################################################
# STOCK POLYGON
#
# Replace the following line with your actual polygon object.
# stock.poly can be either an sf object or a terra SpatVector.
################################################################################

# EXAMPLE:
stock.poly <- terra::vect(paste0("stocks/", stock.name, ".kml"))


################################################################################
# CONVERT POLYGON TO SPATVECTOR AND MATCH THE OCEANOGRAPHIC CRS
################################################################################

if (inherits(stock.poly, "sf") || inherits(stock.poly, "sfc")) stock.poly <- terra::vect(stock.poly)

if (!inherits(stock.poly, "SpatVector")) stop("stock.poly must be an sf/sfc object or a terra SpatVector.")

stock.poly <- terra::project(stock.poly, terra::crs(sst))


################################################################################
# BUILD A MARINE-CELL MASK INSIDE THE STOCK POLYGON
#
# SST is used here only to identify valid ocean cells. Candidate locations are
# subsequently checked against ALL oceanographic variables and required dates.
################################################################################

sea.mask <- terra::ifel(!is.na(sst[[1]]), 1, NA)

sea.stock <- terra::crop(sea.mask, stock.poly)

sea.stock <- terra::mask(sea.stock, stock.poly)

n.valid.cells <- terra::global(!is.na(sea.stock), "sum", na.rm = TRUE)[1, 1]

cat("\nMarine raster cells available inside stock polygon:", n.valid.cells, "\n")

if (!is.finite(n.valid.cells) ) stop("Too few valid marine cells inside the stock polygon.")

n.random <- ifelse(n.valid.cells < n.random, n.valid.cells, n.random) 
################################################################################
# RANDOMLY SELECT MARINE CELLS
#
# Sampling raster cells rather than arbitrary coordinates avoids selecting
# terrestrial locations and prevents several random points from representing
# different positions within exactly the same raster cell.
################################################################################

set.seed(seed.spatial)

random.points <- terra::spatSample(sea.stock, size = n.random, method = "random",
                                   na.rm = TRUE, as.points = TRUE, values = FALSE)

random.xy <- terra::crds(random.points)

random.locations <- data.frame(random_id = seq_len(nrow(random.xy)),
                               longitude = random.xy[, 1],
                               latitude = random.xy[, 2])


################################################################################
# FUNCTION TO EXTRACT ALL MONTHLY OCEANOGRAPHIC SERIES AT ONE RANDOM LOCATION
################################################################################

extract_random_ocean <- function(longitude, latitude, chla, uo, vo, sss, sst, fechas, prefix) {
  
  xy <- data.frame(longitude = longitude, latitude = latitude)
  
  x.chla <- as.numeric(unlist(terra::extract(chla, xy, method = "simple", ID = FALSE), use.names = FALSE))
  
  x.uo <- as.numeric(unlist(terra::extract(uo, xy, method = "simple", ID = FALSE), use.names = FALSE))
  
  x.vo <- as.numeric(unlist(terra::extract(vo, xy, method = "simple", ID = FALSE), use.names = FALSE))
  
  x.sss <- as.numeric(unlist(terra::extract(sss, xy, method = "simple", ID = FALSE), use.names = FALSE))
  
  x.sst <- as.numeric(unlist(terra::extract(sst, xy, method = "simple", ID = FALSE), use.names = FALSE))
  
  z <- data.frame(fecha = fechas, chla = x.chla, uo = x.uo, vo = x.vo, sss = x.sss, sst = x.sst)
  
  names(z)[-1] <- paste0(prefix, ".", c("chla", "uo", "vo", "sss", "sst"))
  
  z
}


################################################################################
# FUNCTION TO FIT EXACTLY THE SELECTED MODEL AT ONE RANDOM LOCATION
################################################################################

fit_random_location <- function(j) {
  
  random.ocean <- extract_random_ocean(longitude = random.locations$longitude[j],
                                       latitude = random.locations$latitude[j],
                                       chla = chla, uo = uo, vo = vo, sss = sss, sst = sst,
                                       fechas = fechas, prefix = prefix)
  
  random.dat <- make_window_data(lag.months = best.lag,
                                 window.months = best.window,
                                 chla.delay = best.chla.delay,
                                 reference.month = reference.month,
                                 Reclu = Reclu,
                                 ocean = random.ocean,
                                 env.vars = env.vars)
  
  random.dat <- random.dat[random.dat$Year %in% dat.best$Year, , drop = FALSE]
  
  random.dat$Y <- if (use.log.recruitment) log(random.dat$Recruitment) else random.dat$Recruitment
  
  required.columns <- c("Y", selected.vars.fixed)
  
  if (!all(required.columns %in% names(random.dat))) {
    
    return(data.frame(random_id = j, longitude = random.locations$longitude[j],
                      latitude = random.locations$latitude[j], n = NA_integer_,
                      R2OOB = NA_real_, RMSEOOB = NA_real_))
  }
  
  random.rf.dat <- random.dat[, required.columns, drop = FALSE]
  
  if (!all(complete.cases(random.rf.dat))) {
    
    return(data.frame(random_id = j, longitude = random.locations$longitude[j],
                      latitude = random.locations$latitude[j],
                      n = sum(complete.cases(random.rf.dat)),
                      R2OOB = NA_real_, RMSEOOB = NA_real_))
  }
  
  if (nrow(random.rf.dat) != nrow(rf.dat)) {
    
    return(data.frame(random_id = j, longitude = random.locations$longitude[j],
                      latitude = random.locations$latitude[j], n = nrow(random.rf.dat),
                      R2OOB = NA_real_, RMSEOOB = NA_real_))
  }
  
  random.form <- stats::reformulate(selected.vars.fixed, response = "Y")
  
  set.seed(seed.rf + best.index)
  
  random.fit <- randomForestSRC::rfsrc(formula = random.form, data = random.rf.dat,
                                       ntree = ntree.rf, mtry = best.mtry,
                                       nodesize = best.nodesize, nsplit = 10,
                                       importance = "none", forest = FALSE)
  
  random.pred.oob <- as.numeric(random.fit$predicted.oob)
  
  random.obs <- random.rf.dat$Y
  
  random.SST <- sum((random.obs - mean(random.obs, na.rm = TRUE))^2, na.rm = TRUE)
  
  random.R2OOB <- 1 - sum((random.obs - random.pred.oob)^2, na.rm = TRUE) / random.SST
  
  random.RMSEOOB <- sqrt(mean((random.obs - random.pred.oob)^2, na.rm = TRUE))
  
  data.frame(random_id = j, longitude = random.locations$longitude[j],
             latitude = random.locations$latitude[j], n = nrow(random.rf.dat),
             R2OOB = random.R2OOB, RMSEOOB = random.RMSEOOB)
}


################################################################################
# RUN THE 100 RANDOM-LOCATION MODELS
################################################################################

random.results <- dplyr::bind_rows(lapply(seq_len(n.random), fit_random_location))

cat("\nNumber of random locations successfully fitted:",
    sum(is.finite(random.results$R2OOB)), "of", n.random, "\n")


################################################################################
# REMOVE LOCATIONS THAT COULD NOT BE EVALUATED
################################################################################

random.results.valid <- random.results[is.finite(random.results$R2OOB), , drop = FALSE]

if (nrow(random.results.valid) == 0L) stop("None of the random locations produced a valid OOB R2.")


################################################################################
# EMPIRICAL COMPARISON WITH THE FOCAL MINFST/MAXFST LOCATION
################################################################################

random.mean <- mean(random.results.valid$R2OOB)

random.sd <- sd(random.results.valid$R2OOB)

random.quantiles <- quantile(random.results.valid$R2OOB,
                             probs = c(0.025, 0.05, 0.50, 0.95, 0.975),
                             na.rm = TRUE)

random.percentile <- mean(random.results.valid$R2OOB <= focal.R2OOB) * 100

empirical.p <- (1 + sum(random.results.valid$R2OOB >= focal.R2OOB)) /
  (1 + nrow(random.results.valid))


################################################################################
# REPORT RESULTS
################################################################################

cat("\nSpatial randomization test\n")
cat("Stock:", stock.name, "\n")
cat("Focal FST location:", fst.use, "\n")
cat("Focal R2 OOB:", round(focal.R2OOB, 4), "\n")
cat("Mean random R2 OOB:", round(random.mean, 4), "\n")
cat("SD random R2 OOB:", round(random.sd, 4), "\n")
cat("Median random R2 OOB:", round(random.quantiles["50%"], 4), "\n")
cat("95% random interval:", round(random.quantiles["2.5%"], 4), "-",
    round(random.quantiles["97.5%"], 4), "\n")
cat("Percentile of focal location:", round(random.percentile, 1), "%\n")
cat("Empirical upper-tail P:", round(empirical.p, 4), "\n\n")


################################################################################
# HISTOGRAM OF RANDOM-LOCATION OOB R2
################################################################################

hist(random.results.valid$R2OOB, breaks = 30, las = 1,
     xlab = expression(R[OOB]^2),
     main = paste(stock.name, "-", fst.use),
     xlim = range(c(random.results.valid$R2OOB, focal.R2OOB), na.rm = TRUE))

abline(v = focal.R2OOB, col = "red", lwd = 3)

abline(v = random.mean, lty = 2, lwd = 2)

legend("topright",
       legend = c(paste0(fst.use, " R2 OOB = ", sprintf("%.2f", focal.R2OOB)),
                  paste0("Random mean = ", sprintf("%.2f", random.mean)),
                  paste0("Empirical P = ", sprintf("%.3f", empirical.p))),
       col = c("red", "black", NA), lty = c(1, 2, NA), lwd = c(3, 2, NA),
       bty = "n")


################################################################################
# OPTIONAL: BOXPLOT WITH FOCAL VALUE
################################################################################

boxplot(random.results.valid$R2OOB, horizontal = TRUE, las = 1,
        xlab = expression(R[OOB]^2), ylim=c(-1,1),
        main = paste(stock.name, "-", fst.use, "vs random marine locations"))

points(focal.R2OOB, 1, pch = 19, cex = 1.5, col = "red")

################################################################################
# Accumulated local effects
# Accumulated local effects

ale.bins <- 6L

pred_rfsrc_ale <- function(
    object,
    newdata,
    type = "response"
) {
  
  pred <- stats::predict(
    object,
    newdata = as.data.frame(newdata)
  )
  
  as.numeric(pred$predicted)
}

dat.ale <- rf.dat[
  ,
  c("Y", selected.vars),
  drop = FALSE
]

dat.ale <- tibble::as_tibble(dat.ale)

if (!all(c("Y", selected.vars) %in% names(dat.ale))) {
  stop("ALE variables do not match the fitted-model data.")
}

ale.unique <- vapply(
  dat.ale[, selected.vars, drop = FALSE],
  function(x) {
    length(unique(x[!is.na(x)]))
  },
  integer(1)
)

ale.vars <- selected.vars[
  ale.unique >= 2L
]

if (length(ale.vars) == 0L) {
  
  warning("No selected predictor has enough unique values for ALE.")
  
  ale.fit <- NULL
  ale.plots <- NULL
  
} else {
  
  ale.fit <- suppressWarnings(
    ale::ALE(
      model = best.fit,
      x_cols = ale.vars,
      data = dat.ale,
      y_col = "Y",
      pred_fun = pred_rfsrc_ale,
      pred_type = "response",
      max_num_bins = ale.bins,
      output_stats = TRUE,
      p_values = NULL,
      boot_it = 0,
      parallel = 0,
      silent = TRUE
    )
  )
  
  ale.plots <- suppressWarnings(
    plot(ale.fit)
  )
  
  print(
    ale.plots,
    ncol = min(2L, length(ale.vars))
  )
}

# Save complete analysis output

run.info <- list(
  date = Sys.time(),
  R.version = R.version.string,
  platform = R.version$platform,
  session.info = utils::sessionInfo()
)

settings <- list(
  stock = stock.name,
  fila = fila,
  fst.use = fst.use,
  recruitment.age = s.age[fila],
  lags = lags,
  windows = windows,
  chla.delays = chla.delays,
  reference.month = reference.month,
  use.log.recruitment = use.log.recruitment,
  ntree.rf = ntree.rf,
  ntree.tune = ntree.tune,
  maxRuns.boruta = maxRuns.boruta,
  seed.rf = seed.rf,
  seed.spatial = seed.spatial,
  n.random.requested = 1000L,
  n.random.used = n.random
)

focal.location <- list(
  stock = stock.name,
  fst.location = fst.use,
  predicted.fst = points.use$predicted_fst[fila.use],
  longitude = points.use$longitude[fila.use],
  latitude = points.use$latitude[fila.use]
)

best.model.info <- list(
  best.rf = best.rf,
  best.index = best.index,
  selected.vars = selected.vars,
  mtry = best.mtry,
  nodesize = best.nodesize,
  lag = best.rf$lag,
  window = best.rf$window,
  chla.delay = best.rf$chla_delay,
  effective.chla.lag = best.rf$lag - best.rf$chla_delay,
  R2.fit = R2.fit,
  R2.oob = R2.oob,
  RMSE.fit = RMSE.fit,
  RMSE.oob = RMSE.oob
)

recruitment.plot.data <- data.frame(
  Year = dat.best$Year,
  observed = observed,
  predicted = predicted,
  predicted.oob = predicted.oob
)

randomization.info <- list(
  focal.R2OOB = focal.R2OOB,
  random.results = random.results,
  random.results.valid = random.results.valid,
  random.mean = random.mean,
  random.sd = random.sd,
  random.quantiles = random.quantiles,
  random.percentile = random.percentile,
  empirical.p = empirical.p,
  random.locations = random.locations
)

ale.info <- list(
  variables = ale.vars,
  unique.values = ale.unique,
  data = dat.ale,
  fit = ale.fit,
  plots = ale.plots,
  bins = ale.bins
)

salida <- list(
  run.info = run.info,
  settings = settings,
  focal.location = focal.location,
  recruitment = Reclu,
  ocean = ocean,
  parameter.grid = parameter.grid,
  valid.years = valid.years,
  lag.window.data = lag.window.data,
  rf.results = rf.results,
  best.model.info = best.model.info,
  best.fit = best.fit,
  rf.dat = rf.dat,
  dat.best = dat.best,
  recruitment.plot.data = recruitment.plot.data,
  randomization = randomization.info,
  ale = ale.info
)

output.file <- paste0(
  "salida.",
  stock.name,
  ".rds"
)

saveRDS(
  salida,
  file = output.file,
  compress = "xz"
)

cat("\nRDS output saved as:", output.file, "\n")
cat(
  "File size:",
  round(file.info(output.file)$size / 1024^2, 2),
  "MB\n"
)

################################################################################
# OPTIONAL: SAVE RANDOMIZATION RESULTS
################################################################################

# write.csv(random.results,
#           paste0("RF_spatial_randomization_", fst.use, "_", stock.name, ".csv"),
#           row.names = FALSE)
}