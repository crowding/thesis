#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(plyr)
  library(reshape2)
  library(vadr)
  library(rstan)
  library(lubridate)
  source("library.R")
})

#files <- dir("models", "*.fit.RData", full.names=TRUE)

main <- function(out, ..., cutoff = as.POSIXct(0)) {
  files <- c(...)
  timings <- collect_times(files)
  save(file=out, timings)
}

collect_times <- function(files) {
  ldply(files, get_times_and_overshoots, .parallel=TRUE)
}

run_partial <- function(out="models/timings_partial.RData",
                        files=dir("models", "*.fit.RData", full.names=TRUE),
                        cutoff,
                        cutoff.interval=as.difftime(20, units="hours")) {
  info <- file.info(files)
  if (missing(cutoff)) {
    files <- chain(
      info, name_rows,
      arrange(desc(mtime)),
      mutate(interval=c(as.difftime(0, units="secs"), -diff(mtime)),
             minterval = cummax(interval)),
      subset(as.difftime(minterval, units="secs") < cutoff.interval),
      .$.rownames)
  } else {
    files <- rownames(subset(info, mtime > cutoff))
  }
  timings <- collect_times(files)
  save(file=out, timings)
}

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
    version <-
    cbind(fit=NA, optimized=NA, overshoot, time)
  })
}

both_versions <- function(
  versions=data.frame(version=c("1.3.0", "2.0.1"),
                      file=c("models/timings.1.3.0.RData",
                             "models/timings_partial.RData"))) {
  mdply(versions, function(version, file) {
    e <- load2env(as.character(file))
    e$timings
  })
}

bad <- function(n=5) {
  data <- both_versions()
  chain(data,
        dcast(model_name + subject~version, value.var="time"),
        mutate(slowdown=`2.0.1`/`1.3.0`),
        subset(is.finite(slowdown)),
        arrange(desc(slowdown)),
        head(20))
}

run_as_command()


