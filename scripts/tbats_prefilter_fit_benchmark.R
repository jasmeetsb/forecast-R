#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(forecast)
  library(jsonlite)
})

parse_args <- function(args) {
  out <- list(
    series_json = NULL,
    seasonal_periods = NULL,
    k_vector = NULL,
    repeat_count = 1L,
    top_k = 2L
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--series-json") {
      i <- i + 1L
      out$series_json <- args[[i]]
    } else if (key == "--seasonal-periods") {
      i <- i + 1L
      out$seasonal_periods <- args[[i]]
    } else if (key == "--k-vector") {
      i <- i + 1L
      out$k_vector <- args[[i]]
    } else if (key == "--repeat") {
      i <- i + 1L
      out$repeat_count <- as.integer(args[[i]])
    } else if (key == "--top-k") {
      i <- i + 1L
      out$top_k <- as.integer(args[[i]])
    } else {
      stop(sprintf("Unknown argument: %s", key))
    }
    i <- i + 1L
  }
  out
}

parse_numeric_list <- function(raw) {
  if (is.null(raw) || !nzchar(raw)) {
    return(numeric())
  }
  as.numeric(strsplit(raw, ",", fixed = TRUE)[[1]])
}

parse_integer_list <- function(raw) {
  if (is.null(raw) || !nzchar(raw)) {
    return(integer())
  }
  as.integer(strsplit(raw, ",", fixed = TRUE)[[1]])
}

load_series_json <- function(path) {
  raw <- jsonlite::read_json(path, simplifyVector = FALSE)
  if (is.null(raw$series) || length(raw$series) == 0) {
    stop(sprintf("No series found in %s", path))
  }
  list(
    series = lapply(raw$series, function(row) as.numeric(unlist(row))),
    frequency = if (!is.null(raw$frequency)) as.integer(raw$frequency) else 1L
  )
}

build_structures <- function(seasonal_periods, k_vector) {
  list(
    list(use_box_cox = FALSE, use_beta = FALSE, use_damped = FALSE),
    list(use_box_cox = TRUE, use_beta = FALSE, use_damped = FALSE),
    list(use_box_cox = FALSE, use_beta = TRUE, use_damped = FALSE),
    list(use_box_cox = TRUE, use_beta = TRUE, use_damped = FALSE),
    list(use_box_cox = FALSE, use_beta = TRUE, use_damped = TRUE),
    list(use_box_cox = TRUE, use_beta = TRUE, use_damped = TRUE)
  )
}

score_initial_structure <- function(y, structure, seasonal_periods, k_vector) {
  if (structure$use_box_cox && any(y <= 0)) {
    return(Inf)
  }

  makeTBATSWMatrix <- getFromNamespace("makeTBATSWMatrix", "forecast")
  updateTBATSGammaBold <- getFromNamespace("updateTBATSGammaBold", "forecast")
  updateTBATSGMatrix <- getFromNamespace("updateTBATSGMatrix", "forecast")
  makeTBATSFMatrix <- getFromNamespace("makeTBATSFMatrix", "forecast")
  calcTBATSFaster <- getFromNamespace("calcTBATSFaster", "forecast")
  calcWTilda <- getFromNamespace("calcWTilda", "forecast")
  makeXMatrix <- getFromNamespace("makeXMatrix", "forecast")

  alpha <- 0.09
  beta.v <- if (structure$use_beta) 0.05 else NULL
  small.phi <- if (structure$use_beta) {
    if (structure$use_damped) 0.999 else 1
  } else {
    NULL
  }
  adj_beta <- if (structure$use_beta) 1 else 0
  gamma.one.v <- rep(0, length(k_vector))
  gamma.two.v <- rep(0, length(k_vector))
  tau <- as.integer(2 * sum(k_vector))

  y_work <- y
  if (structure$use_box_cox) {
    lambda <- BoxCox.lambda(stats::ts(y, frequency = 1), method = "guerrero", lower = -1, upper = 2)
    lambda <- min(max(lambda, 0.01), 0.99)
    y_work <- as.numeric(BoxCox(y, lambda))
  }

  s.vector <- if (tau > 0) numeric(tau) else NULL
  x0_zero <- makeXMatrix(
    l = 0,
    b = if (structure$use_beta) 0 else NULL,
    s.vector = s.vector,
    d.vector = NULL,
    epsilon.vector = NULL
  )$x

  w <- makeTBATSWMatrix(
    smallPhi = small.phi,
    kVector = as.integer(k_vector),
    arCoefs = NULL,
    maCoefs = NULL,
    tau = tau
  )

  gamma.bold <- if (tau > 0) matrix(0, nrow = 1, ncol = tau) else NULL
  if (!is.null(gamma.bold)) {
    updateTBATSGammaBold(
      gammaBold = gamma.bold,
      kVector = as.integer(k_vector),
      gammaOne = gamma.one.v,
      gammaTwo = gamma.two.v
    )
  }

  g <- matrix(
    0,
    nrow = 1 + adj_beta + tau,
    ncol = 1
  )
  updateTBATSGMatrix(
    g = g,
    gammaBold = gamma.bold,
    alpha = alpha,
    beta = beta.v
  )

  F <- makeTBATSFMatrix(
    alpha = alpha,
    beta = beta.v,
    small.phi = small.phi,
    seasonal.periods = seasonal_periods,
    k.vector = as.integer(k_vector),
    gamma.bold.matrix = gamma.bold,
    ar.coefs = NULL,
    ma.coefs = NULL
  )

  D <- F - g %*% w$w.transpose
  n <- length(y_work)
  y_mat <- matrix(y_work, nrow = 1, ncol = n)
  y_hat <- matrix(0, nrow = 1, ncol = n)
  e <- matrix(0, nrow = 1, ncol = n)
  x <- matrix(0, nrow = length(x0_zero), ncol = n)

  calcTBATSFaster(
    y = y_mat,
    yHat = y_hat,
    wTranspose = w$w.transpose,
    F = F,
    x = x,
    g = g,
    e = e,
    xNought = x0_zero
  )

  y.tilda <- e
  w.tilda.transpose <- matrix(0, nrow = n, ncol = ncol(w$w.transpose))
  w.tilda.transpose[1, ] <- w$w.transpose
  w.tilda.transpose <- calcWTilda(
    wTildaTranspose = w.tilda.transpose,
    D = D
  )

  fit <- tryCatch(
    lm.fit(x = w.tilda.transpose, y = as.numeric(t(y.tilda))),
    error = function(e) NULL
  )

  if (is.null(fit) || is.null(fit$coefficients)) {
    x0_est <- matrix(0, nrow = ncol(w.tilda.transpose), ncol = 1)
  } else {
    coefs <- fit$coefficients
    coefs[is.na(coefs)] <- 0
    x0_est <- matrix(coefs, nrow = length(coefs), ncol = 1)
  }

  y_hat2 <- matrix(0, nrow = 1, ncol = n)
  e2 <- matrix(0, nrow = 1, ncol = n)
  x2 <- matrix(0, nrow = length(x0_est), ncol = n)
  calcTBATSFaster(
    y = y_mat,
    yHat = y_hat2,
    wTranspose = w$w.transpose,
    F = F,
    x = x2,
    g = g,
    e = e2,
    xNought = x0_est
  )

  sse <- sum(e2 * e2)
  if (!is.finite(sse) || sse <= 0) {
    return(Inf)
  }
  as.numeric(n * log(sse))
}

fit_structure <- function(y, structure, seasonal_periods, k_vector) {
  fitSpecificTBATS <- getFromNamespace("fitSpecificTBATS", "forecast")
  freq <- max(1L, as.integer(max(seasonal_periods)))
  y_ts <- stats::ts(y, frequency = freq)
  tryCatch(
    fitSpecificTBATS(
      y = y_ts,
      use.box.cox = structure$use_box_cox,
      use.beta = structure$use_beta,
      use.damping = structure$use_damped,
      seasonal.periods = seasonal_periods,
      k.vector = as.integer(k_vector),
      ar.coefs = NULL,
      ma.coefs = NULL,
      bc.lower = 0,
      bc.upper = 1,
      biasadj = FALSE
    ),
    error = function(e) NULL
  )
}

run_benchmark <- function(series, seasonal_periods, k_vector, repeat_count, top_k) {
  structures <- build_structures(seasonal_periods, k_vector)
  best_aic <- NULL
  best_structure <- NULL
  fit_counts <- NULL

  t0 <- proc.time()[["elapsed"]]
  for (rep_idx in seq_len(repeat_count)) {
    best_aic <- rep(Inf, length(series))
    best_structure <- rep(NA_integer_, length(series))
    fit_counts <- integer(length(series))
    for (series_idx in seq_along(series)) {
      y <- series[[series_idx]]
      scores <- vapply(
        structures,
        function(structure) score_initial_structure(y, structure, seasonal_periods, k_vector),
        numeric(1)
      )
      finite_idx <- which(is.finite(scores))
      if (length(finite_idx) == 0) {
        next
      }
      ordered <- finite_idx[order(scores[finite_idx])]
      primary <- ordered[seq_len(min(top_k, length(ordered)))]

      best_fit <- NULL
      best_fit_aic <- Inf
      best_fit_idx <- NA_integer_

      for (structure_idx in primary) {
        fit <- fit_structure(y, structures[[structure_idx]], seasonal_periods, k_vector)
        fit_counts[series_idx] <- fit_counts[series_idx] + 1L
        if (!is.null(fit) && is.finite(fit$AIC) && fit$AIC < best_fit_aic) {
          best_fit <- fit
          best_fit_aic <- fit$AIC
          best_fit_idx <- structure_idx - 1L
        }
      }

      if (is.null(best_fit) && length(primary) < length(ordered)) {
        for (structure_idx in ordered[-seq_along(primary)]) {
          fit <- fit_structure(y, structures[[structure_idx]], seasonal_periods, k_vector)
          fit_counts[series_idx] <- fit_counts[series_idx] + 1L
          if (!is.null(fit) && is.finite(fit$AIC) && fit$AIC < best_fit_aic) {
            best_fit <- fit
            best_fit_aic <- fit$AIC
            best_fit_idx <- structure_idx - 1L
          }
        }
      }

      best_aic[series_idx] <- best_fit_aic
      best_structure[series_idx] <- best_fit_idx
    }
  }
  elapsed <- proc.time()[["elapsed"]] - t0

  list(
    backend = "r",
    method = "tbats-prefit",
    n_series = length(series),
    n_structures = length(structures),
    repeat_count = repeat_count,
    seasonal_periods = seasonal_periods,
    k_vector = k_vector,
    top_k = top_k,
    wall_time_s = as.numeric(elapsed),
    series_per_s = length(series) * repeat_count / elapsed,
    aic_head = as.numeric(best_aic[seq_len(min(3, length(best_aic)))]),
    best_structure_head = as.integer(best_structure[seq_len(min(5, length(best_structure)))]),
    fit_count_mean = mean(fit_counts),
    fit_count_max = max(fit_counts)
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$series_json)) {
  stop("--series-json is required")
}

seasonal_periods <- parse_numeric_list(args$seasonal_periods)
k_vector <- parse_integer_list(args$k_vector)
if (length(seasonal_periods) == 0 || length(k_vector) == 0) {
  stop("--seasonal-periods and --k-vector are required")
}
if (length(seasonal_periods) != length(k_vector)) {
  stop("seasonal_periods and k_vector must have the same length")
}

loaded <- load_series_json(args$series_json)
result <- run_benchmark(
  series = loaded$series,
  seasonal_periods = seasonal_periods,
  k_vector = k_vector,
  repeat_count = args$repeat_count,
  top_k = args$top_k
)
cat(toJSON(result, auto_unbox = TRUE, digits = 10, pretty = TRUE))
