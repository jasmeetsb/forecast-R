#!/usr/bin/env Rscript
# Export M3 competition dataset to JSON for Python benchmark.
library(Mcomp)
library(jsonlite)

outfile <- file.path(Sys.getenv("HOME"), "github/rforecast-conversion/forecast-jax/tests/benchmark/m3_data.json")

series_list <- lapply(seq_along(M3), function(i) {
  s <- M3[[i]]
  list(
    sn = s$sn,
    type = s$type,
    period = s$period,
    frequency = frequency(s$x),
    h = s$h,
    x = as.numeric(s$x),
    xx = as.numeric(s$xx)
  )
})
names(series_list) <- sapply(M3, function(s) s$sn)

cat("Exporting", length(series_list), "M3 series to", outfile, "\n")
write_json(series_list, outfile, auto_unbox = TRUE, digits = 10, pretty = FALSE)
cat("Done. File size:", file.size(outfile), "bytes\n")
