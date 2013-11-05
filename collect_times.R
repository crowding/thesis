#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(plyr)
  library(vadr)
  library(rstan)
  source("library.R")
})

files <- dir("models", "*.fit.RData", full.names=TRUE)

main <- function(out, ...) {
  files <- c(...)

  timings <- ldply(files, get_times_and_overshoots, .parallel=TRUE)
  save(file=out, timings)
}

modelfile <- files[[10]]

get_times_and_overshoots <- function(modelfile) {
  cat(modelfile, "\n")
  e <- load2env(modelfile)
  adply(e$fits, 1, function(row) {
    bind[fit=bind[fit], optimized=optimized, ...=group] <- row
    time <- (if ("time" %in% names(attributes(fit)))
             sum(attr(fit, "time")[c("user.self", "user.child")]) else 0)
    overshoot <- (if ("overshoot" %in% names(attributes(fit)))
                  attr(fit, "overshoot") else
                  data.frame(model_name=e$model_name, group))
    cbind(fit=NA, optimized=NA, overshoot, time)
  })
}

run_as_command()
