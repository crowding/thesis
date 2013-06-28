library(rstan)
library(ptools)
library(plyr)
library(reshape2)
source("library.R")

modelfiles <- dir(pattern=".fit.RData")

plotStanFits <- function(modelfiles) {
  fitdata <- llply(modelfiles, get_samples)
  bind[samples=samples, optimized=optimized] <- reduce_fits(fitdata)
}

get_samples <- function(modelfile) {
  x <- load2env(modelfile)
  mutate(x$fits,
         fit=lapply(fit, as.data.frame),
         model_name=x$model@model_name)
}

unfactor <- colwise(function(x) if (is.factor(x)) as.character(x) else x)

reduce_fits <- function(fits) {
  samples <- ldply(fits, function(x) {
    n <- names(x) %-% c("fit", "optimized")
    frame <- rbind.fill %()% ((Map %<<% x)(
      function(fit, optimized, ...)
      cbind(fit, data.frame(..., stringsAsFactors=FALSE))))
    unfactor(melt(frame, id.vars=n))
  })
  optimized <- ldply(fits, function(x) {
    n <- names(x) %-% c("fit", "optimized")
    frame <- rbind.fill %()% ((Map %<<% x)(
      function(fit, optimized, ...)
      cbind(as.data.frame(optimized),
            data.frame(..., stringsAsFactors=FALSE))))
    unfactor(melt(frame, id.vars=n))
  })
  list(samples=samples, optimized=optimized)
}

main <- function(...) {
  plotStanFits(c(...))
}

violinPlot <- function(samples, optimized) {}

runAsCommand()
