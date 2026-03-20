#!/usr/bin/env Rscript
# Generate golden test fixtures for forecast-jax validation
# Run from forecast-r/: Rscript scripts/generate_golden.R
# Output goes to: ../forecast-jax/tests/golden/

.libPaths(Sys.getenv("R_LIBS_USER"))
library(forecast)
library(jsonlite)

outdir <- file.path(dirname(getwd()), "forecast-jax", "tests", "golden")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

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

cat("Done. Generated", length(list.files(outdir, pattern = "\\.json$")), "golden test files.\n")
