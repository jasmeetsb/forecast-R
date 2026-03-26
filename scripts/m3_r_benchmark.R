#!/usr/bin/env Rscript
# Run R forecast methods on all M3 series and export results as JSONL.
# Usage: Rscript m3_r_benchmark.R [method] [period]
# method: ets, arima, theta, naive (default: all)
# period: YEARLY, QUARTERLY, MONTHLY, OTHER (default: all)

library(Mcomp)
library(forecast)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
method_filter <- if (length(args) >= 1) args[1] else "all"
period_filter <- if (length(args) >= 2) toupper(args[2]) else "all"

outdir <- file.path(Sys.getenv("HOME"), "github/rforecast-conversion/forecast-jax/tests/benchmark/r_results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# sMAPE helper
smape <- function(actual, forecast) {
  denom <- abs(actual) + abs(forecast)
  ratio <- ifelse(denom > 0, abs(actual - forecast) / denom, 0)
  200 * mean(ratio)
}

run_method <- function(method_name, series_list) {
  outfile <- file.path(outdir, paste0(method_name, ".jsonl"))
  # Write header
  header <- list(
    method = method_name,
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    forecast_version = as.character(packageVersion("forecast")),
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    n_series = length(series_list)
  )
  writeLines(toJSON(header, auto_unbox = TRUE), outfile)

  for (i in seq_along(series_list)) {
    s <- series_list[[i]]
    sn <- s$sn
    period <- s$period
    h <- s$h
    freq <- frequency(s$x)
    train <- s$x
    test <- s$xx

    result <- tryCatch({
      t0 <- proc.time()
      if (method_name == "ets") {
        fit <- ets(train)
        fc <- forecast(fit, h = h)
        model_str <- fit$method
        aic_val <- fit$aic
      } else if (method_name == "arima") {
        fit <- auto.arima(train)
        fc <- forecast(fit, h = h)
        model_str <- paste0("ARIMA", paste0("(", paste(arimaorder(fit), collapse=","), ")"))
        aic_val <- fit$aic
      } else if (method_name == "theta") {
        fc <- thetaf(train, h = h)
        model_str <- "Theta"
        aic_val <- NA
      } else if (method_name == "naive") {
        if (freq > 1) {
          fc <- snaive(train, h = h)
          model_str <- "Seasonal naive"
        } else {
          fc <- naive(train, h = h)
          model_str <- "Naive"
        }
        aic_val <- NA
      }
      elapsed <- (proc.time() - t0)[3]
      fc_vals <- as.numeric(fc$mean)
      smape_val <- smape(as.numeric(test), fc_vals)
      # MASE: scale from training set
      scale_err <- mean(abs(diff(train, lag = max(1, freq))))
      mase_val <- if (scale_err > 1e-10) mean(abs(as.numeric(test) - fc_vals)) / scale_err else NA

      list(sn = sn, period = period, status = "OK",
           model_string = model_str, aic = aic_val,
           forecast = fc_vals, smape = smape_val, mase = mase_val,
           time_s = as.numeric(elapsed), error = NULL)
    }, error = function(e) {
      list(sn = sn, period = period, status = "FAILED",
           model_string = NA, aic = NA,
           forecast = rep(NA, h), smape = NA, mase = NA,
           time_s = NA, error = conditionMessage(e))
    })

    line <- toJSON(result, auto_unbox = TRUE, digits = 10)
    write(line, file = outfile, append = TRUE)

    if (i %% 100 == 0) {
      cat(sprintf("[%s] %d/%d series done\n", method_name, i, length(series_list)))
    }
  }
  cat(sprintf("[%s] Complete: %d series → %s\n", method_name, length(series_list), outfile))
}

# Filter series
series <- M3
if (period_filter != "all") {
  series <- series[sapply(series, function(s) toupper(s$period) == period_filter)]
  cat("Filtered to", length(series), "series for period", period_filter, "\n")
}

methods <- if (method_filter == "all") c("naive", "theta", "arima", "ets") else method_filter

for (m in methods) {
  cat(sprintf("Starting %s on %d series...\n", m, length(series)))
  run_method(m, series)
}
cat("All done.\n")
