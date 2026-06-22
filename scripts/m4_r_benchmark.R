#!/usr/bin/env Rscript
# M4 competition benchmark -- R `forecast` baseline for forecast-jax.
#
# Runs ONE (method, period) per invocation on the normalized M4 JSON
# ({sn, period, frequency, h, x, xx}) and writes JSONL matching the M3
# r_results schema so R and forecast-jax are scored on identical inputs
# with identical accuracy measures.
#
# Usage:
#   Rscript m4_r_benchmark.R <method> <PERIOD> --data-file <path> \
#       [--out-dir <dir>] [--limit <N>]
#
#   method : naive theta arima ets tbats bats nnetar
#   PERIOD : YEARLY QUARTERLY MONTHLY WEEKLY DAILY HOURLY
#   --limit: cap series count (smoke testing)
#   --cores: fit series in parallel across N fork workers (mclapply); default 1.
#            Each series is fit on its own core (the standard "forecasting at
#            scale" pattern); a single fit is itself single-threaded.
#
# Output : <out-dir>/<method>_<period>.jsonl
#   default out-dir = sibling forecast-jax/tests/benchmark/py_results/m4_r,
#   else ./py_results/m4_r
#
# Method defs + metrics mirror forecast-jax EXACTLY (see
# m3_benchmark.run_single_inprocess / _run_single_extra and
# forecast_jax.metrics.accuracy):
#   naive  -> snaive() if frequency>1 else naive()
#   theta  -> thetaf()
#   arima  -> auto.arima()   (R defaults: stepwise search)
#   ets    -> ets()          (R defaults)
#   tbats  -> tbats(); bats -> bats()   (ts frequency = seasonal period)
#   nnetar -> set.seed(42); nnetar()
#   sMAPE  = mean(200*|a-f| / (|a|+|f|)), 0 when |a|+|f| == 0
#   MASE   = mean(|a-f|) / mean(|train[m:] - train[:-m]|), m = max(1, frequency)
#
# Requires: forecast, jsonlite

suppressPackageStartupMessages({
  library(forecast)
  library(jsonlite)
  library(parallel)
})

## ---- arg parsing -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript m4_r_benchmark.R <method> <PERIOD> --data-file <path> ",
       "[--out-dir <dir>] [--limit N]")
}
method <- tolower(args[1])
period <- toupper(args[2])
rest <- args[-(1:2)]

get_flag <- function(flag, default = NULL) {
  i <- which(rest == flag)
  if (length(i) == 0) return(default)
  rest[i[1] + 1]
}
data_file <- get_flag("--data-file")
out_dir   <- get_flag("--out-dir")
limit_arg <- get_flag("--limit")
cores_arg <- get_flag("--cores")
if (is.null(data_file)) stop("--data-file is required")
limit <- if (is.null(limit_arg)) Inf else as.integer(limit_arg)
cores <- if (is.null(cores_arg)) 1L else max(1L, as.integer(cores_arg))

valid_methods <- c("naive", "theta", "arima", "ets", "tbats", "bats", "nnetar")
if (!(method %in% valid_methods)) {
  stop(sprintf("Unknown method '%s'; choose one of: %s",
               method, paste(valid_methods, collapse = " ")))
}

## ---- resolve output dir ----------------------------------------------------
all_args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", all_args[grep("^--file=", all_args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
if (is.null(out_dir)) {
  # Prefer sibling forecast-jax/tests/benchmark/py_results/m4_r. The
  # py_results/m4_r leg may not exist yet, so anchor the existence check on
  # the stable benchmark dir and create the rest below.
  bench_dir <- file.path(script_dir, "..", "..", "forecast-jax",
                         "tests", "benchmark")
  out_dir <- if (dir.exists(bench_dir)) file.path(bench_dir, "py_results", "m4_r")
             else file.path(getwd(), "py_results", "m4_r")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(out_dir, sprintf("%s_%s.jsonl", method, tolower(period)))

## ---- load data -------------------------------------------------------------
data <- fromJSON(data_file, simplifyVector = TRUE,
                 simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
series <- unname(data)
if (is.finite(limit)) series <- series[seq_len(min(limit, length(series)))]
n <- length(series)
cat(sprintf("Loaded %d series from %s\n", n, data_file))
cat(sprintf("Running %s on %s -> %s\n", method, period, outfile))

## ---- accuracy measures (mirror forecast_jax.metrics.accuracy) --------------
smape_fn <- function(actual, fc) {
  denom <- abs(actual) + abs(fc)
  ratio <- ifelse(denom > 0, abs(actual - fc) / denom, 0)
  200 * mean(ratio)
}
mase_fn <- function(train, actual, fc, frequency) {
  m <- max(1L, frequency)
  if (length(train) <= m) return(NA_real_)
  scale <- mean(abs(train[(m + 1):length(train)] - train[1:(length(train) - m)]))
  if (!is.finite(scale) || scale < 1e-10) return(NA_real_)
  mean(abs(actual - fc)) / scale
}

## ---- per-series fit + forecast ---------------------------------------------
fit_forecast <- function(method, y, h, frequency) {
  aic <- NA_real_
  if (method == "naive") {
    fc <- if (frequency > 1) snaive(y, h = h) else naive(y, h = h)
    model <- if (frequency > 1) "Seasonal naive" else "Naive"
  } else if (method == "theta") {
    fc <- thetaf(y, h = h)
    model <- "Theta"
  } else if (method == "arima") {
    fit <- auto.arima(y)
    fc <- forecast(fit, h = h)
    ord <- arimaorder(fit)
    model <- sprintf("ARIMA(%d,%d,%d)", ord["p"], ord["d"], ord["q"])
    if (!is.na(ord["P"]) && (ord["P"] + ord["D"] + ord["Q"]) > 0) {
      model <- sprintf("%s(%d,%d,%d)[%d]", model,
                       ord["P"], ord["D"], ord["Q"], ord["Frequency"])
    }
    aic <- tryCatch(as.numeric(fit$aic), error = function(e) NA_real_)
  } else if (method == "ets") {
    fit <- ets(y)
    fc <- forecast(fit, h = h)
    model <- fit$method
    aic <- tryCatch(as.numeric(fit$aic), error = function(e) NA_real_)
  } else if (method == "tbats") {
    # use.parallel=FALSE: tbats' internal cluster collides with the outer
    # mclapply workers (socket-port exhaustion). The model SELECTION is identical
    # either way (parallel only distributes the search), so accuracy is unchanged.
    fit <- tbats(y, use.parallel = FALSE)
    fc <- forecast(fit, h = h)
    model <- tryCatch(trimws(capture.output(print(fit))[1]),
                      error = function(e) "TBATS")
  } else if (method == "bats") {
    fit <- bats(y, use.parallel = FALSE)
    fc <- forecast(fit, h = h)
    model <- tryCatch(trimws(capture.output(print(fit))[1]),
                      error = function(e) "BATS")
  } else if (method == "nnetar") {
    set.seed(42)
    fit <- nnetar(y)
    fc <- forecast(fit, h = h)
    model <- tryCatch(as.character(fit$method), error = function(e) NA_character_)
    if (length(model) == 0 || is.na(model[1])) model <- "NNAR"
  }
  list(mean = as.numeric(fc$mean)[seq_len(h)], model = model, aic = aic)
}

## ---- per-series worker (pure: returns one JSONL row, no shared state) -------
empty_obj <- setNames(list(), character(0))   # serializes to {}
# Reference M3 schema includes `aic` only for these methods; tbats/bats/nnetar
# omit the key entirely.
emit_aic <- method %in% c("naive", "theta", "arima", "ets")

process_one <- function(s) {
  sn <- s$sn
  per <- if (!is.null(s$period)) s$period else period
  freq <- as.integer(s$frequency)
  h <- as.integer(s$h)
  train <- as.numeric(s$x)
  test <- as.numeric(s$xx)
  tryCatch({
    y <- ts(train, frequency = freq)
    t0 <- Sys.time()
    r <- fit_forecast(method, y, h, freq)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    fc <- r$mean
    ma <- mase_fn(train, test, fc, freq)
    out <- list(sn = sn, period = per, status = "OK", model_string = r$model)
    if (emit_aic) out$aic <- if (is.finite(r$aic)) round(r$aic, 6) else NA
    out$forecast <- round(fc, 6)
    out$smape    <- smape_fn(test, fc)
    out$mase     <- if (is.finite(ma)) ma else NA
    out$time_s   <- round(elapsed, 3)   # isolated per-series fit time
    out$error    <- empty_obj
    out
  }, error = function(e) {
    out <- list(sn = sn, period = per, status = "FAILED", model_string = NA)
    if (emit_aic) out$aic <- NA
    out$forecast <- NA
    out$smape    <- NA
    out$mase     <- NA
    out$time_s   <- NA
    out$error    <- conditionMessage(e)
    out
  })
}

## ---- run (series fit in parallel across cores when --cores > 1) -------------
cat(sprintf("Fitting %d series with %d core(s)...\n", n, cores))
t_start <- Sys.time()
if (cores > 1L) {
  # mc.preschedule=FALSE => dynamic dispatch, balances heavy-tailed fit times.
  rows <- mclapply(series, process_one, mc.cores = cores, mc.preschedule = FALSE)
  for (i in seq_along(rows)) {        # rescue any worker that died (segfault/OOM)
    if (!is.list(rows[[i]]) || is.null(rows[[i]]$status)) {
      s <- series[[i]]
      fr <- list(sn = s$sn,
                 period = if (!is.null(s$period)) s$period else period,
                 status = "FAILED", model_string = NA)
      if (emit_aic) fr$aic <- NA
      fr$forecast <- NA; fr$smape <- NA; fr$mase <- NA; fr$time_s <- NA
      fr$error <- paste("worker died:", paste(as.character(rows[[i]]), collapse = " "))
      rows[[i]] <- fr
    }
  }
} else {
  rows <- lapply(series, process_one)
}
elapsed_all <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

## ---- write JSONL (header + one row per series, in input order) -------------
header <- list(
  method = method,
  r_version = as.character(getRversion()),
  forecast_version = as.character(packageVersion("forecast")),
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  n_series = n
)
con <- file(outfile, "w")
writeLines(toJSON(header, auto_unbox = TRUE), con)
for (row in rows) {
  writeLines(toJSON(row, auto_unbox = TRUE, na = "null", null = "null", digits = 10), con)
}
close(con)

## ---- summary ---------------------------------------------------------------
statuses <- vapply(rows, function(r) as.character(r$status), character(1))
n_ok <- sum(statuses == "OK"); n_fail <- n - n_ok
smapes <- vapply(rows, function(r) {
  v <- r$smape
  if (is.numeric(v) && length(v) == 1 && is.finite(v)) v else NA_real_
}, numeric(1))
cat(sprintf("\n  [%s/%s] DONE: %d OK, %d failed | avg sMAPE=%.2f | %.1fs wall (%d core)\n",
            method, tolower(period), n_ok, n_fail,
            mean(smapes, na.rm = TRUE), elapsed_all, cores))
