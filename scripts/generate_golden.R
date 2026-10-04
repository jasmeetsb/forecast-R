#!/usr/bin/env Rscript
# Generate golden test fixtures for forecast-jax validation
# Run from forecast-r/: Rscript scripts/generate_golden.R [OUTDIR]
# Output goes to OUTDIR (default: ../forecast-jax/tests/golden/)
#
# Regenerates all 71 fixtures in forecast-jax/tests/golden/. The sections from
# "Naive / SNaive / RWF / Mean" through "ETS" are the original Phase 1-2 script
# (forecast-r commit f139f935) and are unchanged. The sections after them
# rebuild the generators for the fixtures added later in forecast-jax:
#   df44710  ets_engine
#   5c27a41  auto_arima, ndiffs
#   b10c3b1  findfrequency, fourier, ma, mstl, seasonaldummy, stlf
#   57c867c  tslm, splinef, arfima, dshw
# Requires this repository's forecast package installed in R_LIBS_USER (its
# R/C++ source is identical to the CRAN 9.0.2 release) and jsonlite >= 1.8.5.
# With R 4.3.3, forecast 9.0.2.9000 built from commit 435bf15c and jsonlite
# 1.8.8, every file this script writes is byte-identical to the committed
# fixture (verified 2026-10-04).

.libPaths(Sys.getenv("R_LIBS_USER"))
library(forecast)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) {
  args[1]
} else {
  file.path(dirname(getwd()), "forecast-jax", "tests", "golden")
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- Provenance (printed to stdout only; nothing is added to the fixtures) ---
git_head <- suppressWarnings(tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) character(0)
))
if (length(git_head) != 1 || !is.null(attr(git_head, "status"))) {
  git_head <- "(working directory is not a git checkout)"
}
dep_pkgs <- c("jsonlite", "fracdiff", "urca", "Rcpp", "RcppArmadillo")
cat("R version:         ", R.version.string, "\n")
cat("forecast version:  ", format(packageVersion("forecast")), "\n")
cat("forecast path:     ", find.package("forecast"), "\n")
cat("forecast built:    ", packageDescription("forecast")$Built, "\n")
cat("dependencies:      ", paste(dep_pkgs, vapply(dep_pkgs, function(p) {
  format(packageVersion(p))
}, ""), collapse = ", "), "\n")
cat("working directory: ", getwd(), "\n")
cat("git HEAD:          ", git_head, "\n")

cat("Generating golden test fixtures to:", outdir, "\n")

# --- Datasets ---
datasets <- list(
  airpass = AirPassengers,
  usdeaths = USAccDeaths,
  www = WWWusage,
  lynx = lynx,
  nile = Nile
)

# --- Helper: serialize forecast result ---
serialize_forecast <- function(fc, h_val) {
  list(
    mean = as.numeric(fc$mean),
    lower_80 = as.numeric(fc$lower[, "80%"]),
    upper_80 = as.numeric(fc$upper[, "80%"]),
    lower_95 = as.numeric(fc$lower[, "95%"]),
    upper_95 = as.numeric(fc$upper[, "95%"])
  )
}

# ============================================================
# Naive / SNaive / RWF / Mean
# ============================================================
for (name in names(datasets)) {
  y <- datasets[[name]]
  freq <- frequency(y)
  h <- if (freq > 1) 2 * freq else 10

  # --- naive ---
  fc <- naive(y, h = h)
  result <- list(
    function_name = "naive",
    dataset = name,
    input = list(y = as.numeric(y), frequency = freq, h = h),
    output = list(
      method = fc$method,
      fitted = as.numeric(fc$fitted),
      residuals = as.numeric(fc$residuals),
      forecast = serialize_forecast(fc, h)
    )
  )
  write_json(result, file.path(outdir, paste0("naive_", name, ".json")),
             auto_unbox = TRUE, digits = 10)

  # --- snaive (only for seasonal data) ---
  if (freq > 1) {
    fc <- snaive(y, h = h)
    result <- list(
      function_name = "snaive",
      dataset = name,
      input = list(y = as.numeric(y), frequency = freq, h = h),
      output = list(
        method = fc$method,
        fitted = as.numeric(fc$fitted),
        residuals = as.numeric(fc$residuals),
        forecast = serialize_forecast(fc, h)
      )
    )
    write_json(result, file.path(outdir, paste0("snaive_", name, ".json")),
               auto_unbox = TRUE, digits = 10)
  }

  # --- rwf with drift ---
  fc <- rwf(y, h = h, drift = TRUE)
  result <- list(
    function_name = "rwf_drift",
    dataset = name,
    input = list(y = as.numeric(y), frequency = freq, h = h, drift = TRUE),
    output = list(
      method = fc$method,
      sigma2 = fc$model$sigma2,
      drift = fc$model$par$drift,
      drift_se = fc$model$par$drift.se,
      fitted = as.numeric(fc$fitted),
      residuals = as.numeric(fc$residuals),
      forecast = serialize_forecast(fc, h)
    )
  )
  write_json(result, file.path(outdir, paste0("rwf_drift_", name, ".json")),
             auto_unbox = TRUE, digits = 10)

  # --- meanf ---
  fc <- meanf(y, h = h)
  result <- list(
    function_name = "meanf",
    dataset = name,
    input = list(y = as.numeric(y), frequency = freq, h = h),
    output = list(
      method = fc$method,
      fitted = as.numeric(fc$fitted),
      residuals = as.numeric(fc$residuals),
      forecast = serialize_forecast(fc, h)
    )
  )
  write_json(result, file.path(outdir, paste0("meanf_", name, ".json")),
             auto_unbox = TRUE, digits = 10)
}

# ============================================================
# BoxCox.lambda
# ============================================================
for (name in names(datasets)) {
  y <- datasets[[name]]
  if (min(y, na.rm = TRUE) > 0) {
    lam_g <- BoxCox.lambda(y, method = "guerrero")
    lam_l <- BoxCox.lambda(y, method = "loglik")
    result <- list(
      function_name = "BoxCox.lambda",
      dataset = name,
      input = list(y = as.numeric(y)),
      output = list(
        lambda_guerrero = lam_g,
        lambda_loglik = lam_l
      )
    )
    write_json(result, file.path(outdir, paste0("boxcox_lambda_", name, ".json")),
               auto_unbox = TRUE, digits = 10)
  }
}

# ============================================================
# Accuracy (training set) for naive methods
# ============================================================
for (name in names(datasets)) {
  y <- datasets[[name]]
  fc <- naive(y, h = 10)
  acc <- accuracy(fc)
  result <- list(
    function_name = "accuracy_naive",
    dataset = name,
    output = list(
      ME = acc[1, "ME"],
      RMSE = acc[1, "RMSE"],
      MAE = acc[1, "MAE"],
      MPE = acc[1, "MPE"],
      MAPE = acc[1, "MAPE"],
      MASE = acc[1, "MASE"],
      ACF1 = acc[1, "ACF1"]
    )
  )
  write_json(result, file.path(outdir, paste0("accuracy_naive_", name, ".json")),
             auto_unbox = TRUE, digits = 10)
}

# ============================================================
# ETS (for Phase 2 — generate now so they're ready)
# ============================================================
for (name in names(datasets)) {
  y <- datasets[[name]]
  fit <- ets(y)
  fc <- forecast(fit, h = 24)
  result <- list(
    function_name = "ets",
    dataset = name,
    input = list(y = as.numeric(y), frequency = frequency(y), model = "ZZZ"),
    output = list(
      method = fit$method,
      par = as.list(fit$par),
      initstate = as.numeric(fit$initstate),
      fitted = as.numeric(fitted(fit)),
      residuals = as.numeric(residuals(fit)),
      aic = fit$aic,
      aicc = fit$aicc,
      bic = fit$bic,
      sigma2 = fit$sigma2,
      loglik = fit$loglik,
      components = as.character(fit$components),
      states_dim = dim(fit$states),
      forecast = serialize_forecast(fc, 24)
    )
  )
  write_json(result, file.path(outdir, paste0("ets_", name, ".json")),
             auto_unbox = TRUE, digits = 10)
}

# ============================================================
# ETS engine replay: ets() (forecast-jax df44710)
# Records the parameters and full initial state of ets(y) as input, and the
# fitted values, residuals, state matrix, lik = -2 * loglik and sigma2 as
# output. The integer codes match pegelsresid.C(): error A=1, M=2;
# trend/season N=0, A=1, M=2. Unused parameters are written the way
# pegelsresid.C() sets them (beta = 0, gamma = 0, phi = 1).
# digits = 15 gives 16-17 significant digits (jsonlite >= 1.8.5), enough to
# replay the recursion exactly. digits = NA would not: since jsonlite 1.8.7 it
# caps output at 15 digits.
# ============================================================
for (name in c("airpass", "nile", "www")) {
  y <- datasets[[name]]
  fit <- ets(y)
  comp <- fit$components
  error_type <- switch(comp[1], A = 1L, M = 2L)
  trend_type <- switch(comp[2], N = 0L, A = 1L, M = 2L)
  season_type <- switch(comp[3], N = 0L, A = 1L, M = 2L)
  damped <- as.logical(comp[4])
  result <- list(
    function_name = "ets_engine",
    dataset = name,
    input = list(
      y = as.numeric(y),
      frequency = frequency(y),
      error_type = error_type,
      trend_type = trend_type,
      season_type = season_type,
      damped = damped,
      alpha = unname(fit$par["alpha"]),
      beta = if (trend_type > 0) unname(fit$par["beta"]) else 0,
      gamma = if (season_type > 0) unname(fit$par["gamma"]) else 0,
      phi = if (damped) unname(fit$par["phi"]) else 1,
      init_state = as.numeric(fit$states[1, ])
    ),
    output = list(
      method = fit$method,
      fitted = as.numeric(fitted(fit)),
      residuals = as.numeric(residuals(fit)),
      # jsonlite writes a one-column state matrix (ts) as a flat array and a
      # multi-column one as an array of rows
      states = fit$states,
      lik = -2 * fit$loglik,
      sigma2 = fit$sigma2
    )
  )
  write_json(result, file.path(outdir, paste0("ets_engine_", name, ".json")),
             auto_unbox = TRUE, digits = 15)
}

# ============================================================
# ARIMA: auto.arima(), forecast(), ndiffs(), nsdiffs() (forecast-jax 5c27a41)
# ============================================================
for (name in names(datasets)) {
  y <- datasets[[name]]
  freq <- frequency(y)
  h <- if (freq > 1) 2 * freq else 10

  # --- auto.arima + forecast ---
  # The ndiffs/nsdiffs fields come from the raw series, so they can differ from
  # the d chosen inside auto.arima. For usdeaths auto.arima applies ndiffs
  # after the seasonal difference and picks d = 1, while ndiffs(y) is 0.
  fit <- auto.arima(y)
  fc <- forecast(fit, h = h)
  result <- list(
    function_name = "auto_arima",
    dataset = name,
    input = list(y = as.numeric(y), frequency = freq),
    output = list(
      method = fc$method,
      order = fit$arma[c(1, 6, 2)],              # p, d, q
      seasonal_order = fit$arma[c(3, 7, 4, 5)],  # P, D, Q, period
      par = as.list(coef(fit)),
      aic = fit$aic,
      aicc = fit$aicc,
      bic = fit$bic,
      sigma2 = fit$sigma2,
      loglik = fit$loglik,
      ndiffs = ndiffs(y),
      nsdiffs = if (freq > 1) nsdiffs(y) else 0,
      fitted = as.numeric(fitted(fit)),
      residuals = as.numeric(residuals(fit)),
      forecast = c(list(h = h), serialize_forecast(fc, h))
    )
  )
  write_json(result, file.path(outdir, paste0("auto_arima_", name, ".json")),
             auto_unbox = TRUE, digits = 10)

  # --- ndiffs (KPSS default and ADF) / nsdiffs (seasonal-strength default) ---
  result <- list(
    function_name = "ndiffs",
    dataset = name,
    input = list(y = as.numeric(y), frequency = freq),
    output = list(
      ndiffs_kpss = ndiffs(y, test = "kpss"),
      ndiffs_adf = ndiffs(y, test = "adf"),
      nsdiffs = if (freq > 1) nsdiffs(y) else 0
    )
  )
  write_json(result, file.path(outdir, paste0("ndiffs_", name, ".json")),
             auto_unbox = TRUE, digits = 10)
}

# ============================================================
# Decomposition and seasonal regressors (forecast-jax b10c3b1)
# ============================================================

# --- findfrequency() ---
for (name in c("airpass", "usdeaths", "www", "lynx")) {
  y <- datasets[[name]]
  result <- list(
    function_name = "findfrequency",
    dataset = name,
    input = list(y = as.numeric(y)),
    output = list(frequency = findfrequency(y))
  )
  write_json(result, file.path(outdir, paste0("findfrequency_", name, ".json")),
             auto_unbox = TRUE, digits = 10)
}

# --- fourier(): in-sample matrix (row-major) plus in-sample and h = 24 future
#     shapes, K = 1, 2, 3 ---
y <- AirPassengers
for (K in 1:3) {
  ft <- fourier(y, K = K)
  ft_future <- fourier(y, K = K, h = 24)
  result <- list(
    function_name = "fourier",
    dataset = "airpass",
    input = list(y = as.numeric(y), frequency = frequency(y), K = K, h = 24),
    output = list(
      insample = unname(ft),
      insample_nrow = nrow(ft),
      insample_ncol = ncol(ft),
      future_nrow = nrow(ft_future),
      future_ncol = ncol(ft_future)
    )
  )
  write_json(result, file.path(outdir, paste0("fourier_airpass_K", K, ".json")),
             auto_unbox = TRUE, digits = 10)
}

# --- seasonaldummy() ---
y <- AirPassengers
sd_mat <- seasonaldummy(y)
result <- list(
  function_name = "seasonaldummy",
  dataset = "airpass",
  input = list(y = as.numeric(y), frequency = frequency(y)),
  output = list(matrix = unname(sd_mat), nrow = nrow(sd_mat), ncol = ncol(sd_mat))
)
write_json(result, file.path(outdir, "seasonaldummy_airpass.json"),
           auto_unbox = TRUE, digits = 10)

# --- ma(): even order 12 (centred 2x12) and odd order 7; NA ends -> "NA" ---
y <- AirPassengers
result <- list(
  function_name = "ma",
  dataset = "airpass",
  input = list(y = as.numeric(y), frequency = frequency(y)),
  output = list(
    ma12 = as.numeric(ma(y, order = 12)),
    ma7 = as.numeric(ma(y, order = 7))
  )
)
write_json(result, file.path(outdir, "ma_airpass.json"),
           auto_unbox = TRUE, digits = 10)

# --- mstl() ---
# The committed fixtures label the dataset with the snake_case form of the R
# dataset name ("_air_passengers", "_u_s_acc_deaths"). The file names use the
# short names. Seasonal columns follow the fixed fields in the order mstl()
# returns them.
mstl_sets <- c(airpass = "AirPassengers", usdeaths = "USAccDeaths")
for (name in names(mstl_sets)) {
  y <- get(mstl_sets[[name]])
  fit <- mstl(y)
  out <- list(
    trend = as.numeric(fit[, "Trend"]),
    remainder = as.numeric(fit[, "Remainder"]),
    ncol = ncol(fit),
    colnames = colnames(fit)
  )
  for (col in grep("^Seasonal", colnames(fit), value = TRUE)) {
    out[[col]] <- as.numeric(fit[, col])
  }
  result <- list(
    function_name = "mstl",
    dataset = tolower(gsub("([A-Z])", "_\\1", mstl_sets[[name]])),
    input = list(y = as.numeric(y), frequency = frequency(y)),
    output = out
  )
  write_json(result, file.path(outdir, paste0("mstl_", name, ".json")),
             auto_unbox = TRUE, digits = 10)
}

# --- stlf(): forecast fields are flat in output, not nested ---
for (name in c("airpass", "usdeaths")) {
  y <- datasets[[name]]
  fc <- stlf(y, h = 24)
  result <- list(
    function_name = "stlf",
    dataset = name,
    input = list(y = as.numeric(y), frequency = frequency(y), h = 24),
    output = c(
      list(method = fc$method),
      serialize_forecast(fc, 24),
      list(
        fitted = as.numeric(fc$fitted),
        residuals = as.numeric(fc$residuals)
      )
    )
  )
  write_json(result, file.path(outdir, paste0("stlf_", name, ".json")),
             auto_unbox = TRUE, digits = 10)
}

# ============================================================
# tslm / splinef / arfima / dshw (forecast-jax 57c867c)
# forecast.lm() returns interval matrices without column names, and dshw()
# returns no intervals. This serializer therefore takes interval columns by
# position (default level = c(80, 95)) and writes [] when they are missing.
# ============================================================
serialize_forecast_d1 <- function(fc) {
  interval <- function(m, j) if (is.null(m)) numeric(0) else as.numeric(m[, j])
  list(
    mean = as.numeric(fc$mean),
    lower_80 = interval(fc$lower, 1),
    upper_80 = interval(fc$upper, 1),
    lower_95 = interval(fc$lower, 2),
    upper_95 = interval(fc$upper, 2)
  )
}

# --- tslm() + forecast.lm(): trend + season if seasonal, else trend ---
for (name in c("airpass", "usdeaths", "www", "nile")) {
  y <- datasets[[name]]
  freq <- frequency(y)
  h <- if (freq > 1) 2 * freq else 10
  formula_rhs <- if (freq > 1) "trend + season" else "trend"
  fit <- tslm(as.formula(paste("y ~", formula_rhs)))
  fc <- forecast(fit, h = h)
  fit_summary <- summary(fit)
  result <- list(
    function_name = "tslm",
    dataset = name,
    input = list(y = as.numeric(y), frequency = freq, formula = formula_rhs, h = h),
    output = list(
      method = fc$method,
      coefficients = as.list(coef(fit)),
      sigma = fit_summary$sigma,
      r_squared = fit_summary$r.squared,
      adj_r_squared = fit_summary$adj.r.squared,
      fitted = as.numeric(fitted(fit)),
      residuals = as.numeric(residuals(fit)),
      forecast = serialize_forecast_d1(fc)
    )
  )
  write_json(result, file.path(outdir, paste0("tslm_", name, ".json")),
             auto_unbox = TRUE, digits = 10)
}

# --- splinef(): AirPassengers with lambda = 0 (log), the others untransformed ---
splinef_sets <- list(
  airpass_boxcox = list(y = AirPassengers, lambda = 0),
  nile = list(y = Nile),
  usdeaths = list(y = USAccDeaths),
  www = list(y = WWWusage)
)
for (name in names(splinef_sets)) {
  y <- splinef_sets[[name]]$y
  lambda <- splinef_sets[[name]]$lambda
  freq <- frequency(y)
  h <- if (freq > 1) 2 * freq else 10
  fc <- splinef(y, h = h, lambda = lambda)
  input <- list(y = as.numeric(y), frequency = freq, h = h)
  if (!is.null(lambda)) {
    input$lambda <- lambda
  }
  result <- list(
    function_name = "splinef",
    dataset = name,
    input = input,
    output = list(
      method = fc$method,
      fitted = as.numeric(fc$fitted),
      residuals = as.numeric(fc$residuals),
      forecast = serialize_forecast_d1(fc)
    )
  )
  write_json(result, file.path(outdir, paste0("splinef_", name, ".json")),
             auto_unbox = TRUE, digits = 10)
}

# --- arfima() + forecast.fracdiff() ---
# ar and ma are kept in fracdiff's sign convention. Empty ar/ma vectors are
# written as NULL, which jsonlite encodes as {} (the committed form). sigma2 is
# also {} because arfima objects have no sigma2 element.
for (name in c("www", "lynx", "nile")) {
  y <- datasets[[name]]
  h <- 10
  fit <- arfima(y)
  fc <- forecast(fit, h = h)
  result <- list(
    function_name = "arfima",
    dataset = name,
    input = list(y = as.numeric(y), frequency = frequency(y), h = h),
    output = list(
      d = fit$d,
      ar = if (length(fit$ar) > 0) fit$ar else NULL,
      ma = if (length(fit$ma) > 0) fit$ma else NULL,
      sigma2 = fit$sigma2,
      fitted = as.numeric(fitted(fit)),
      residuals = as.numeric(residuals(fit)),
      forecast = serialize_forecast_d1(fc)
    )
  )
  write_json(result, file.path(outdir, paste0("arfima_", name, ".json")),
             auto_unbox = TRUE, digits = 10)
}

# --- dshw(): synthetic double-seasonal series and the first 3 weeks of taylor ---
# Synthetic: linear level x multiplicative period-4 and period-12 patterns x
# log-normal noise (seed 42, t = 1..120). This is the R form of
# _make_double_seasonal() in forecast-jax tests/unit/test_dshw_unit.py.
set.seed(42)
t <- 1:120
dshw_synthetic <- (10 + 0.02 * t) * (1 + 0.3 * sin(2 * pi * t / 4)) *
  (1 + 0.2 * cos(2 * pi * t / 12)) * exp(rnorm(120, 0, 0.05))
dshw_sets <- list(
  synthetic = list(label = "synthetic", y = dshw_synthetic,
                   period1 = 4, period2 = 12, h = 12),
  taylor = list(label = "taylor_3wk", y = taylor[1:1008],  # 3 x 336 half-hours
                period1 = 48, period2 = 336, h = 48)
)
for (name in names(dshw_sets)) {
  ds <- dshw_sets[[name]]
  fc <- dshw(ds$y, period1 = ds$period1, period2 = ds$period2, h = ds$h)
  result <- list(
    function_name = "dshw",
    dataset = ds$label,
    input = list(y = as.numeric(ds$y), period1 = ds$period1,
                 period2 = ds$period2, h = ds$h),
    output = list(
      method = fc$method,
      fitted = as.numeric(fc$fitted),
      residuals = as.numeric(fc$residuals),
      forecast = serialize_forecast_d1(fc)
    )
  )
  write_json(result, file.path(outdir, paste0("dshw_", name, ".json")),
             auto_unbox = TRUE, digits = 10)
}

cat("Done. Generated", length(list.files(outdir, pattern = "\\.json$")), "golden test files.\n")
