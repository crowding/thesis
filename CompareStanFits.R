library(rstan)
library(ptools)
library(plyr)
library(reshape2)
library(ggplot2)
source("library.R")
theme_set(theme_bw())
theme_update(panel.grid.major=element_blank(),
             panel.grid.minor=element_blank())

modelfiles <- dir(pattern=".fit.RData")

plotStanFits <- function(modelfiles=dir(pattern=".fit.RData")) {
  fitdata <- llply(modelfiles, get_samples)
  bind[samples=samples, optimized=optimized] <- collect_model_data(fitdata)
  violinPlot(samples, optimized)
}

get_samples <- function(modelfile) {
  #load fit file, extract samples from each fit,
  #remember model name, return named list
  x <- load2env(modelfile)
  mutate(x$fits,
         fit=lapply(fit, as.data.frame),
         model_name=x$model@model_name)
}

unfactor <- colwise(function(x) if (is.factor(x)) as.character(x) else x)

#unpack a fitlist into two long-format data frames
collect_model_data <- function(model.list) {
  #fitlist is a list with one item per model
  list(samples = ldply(model.list, collect_fit_data, "fit"),
       optimized = ldply(model.list, collect_fit_data, "optimized"))
}

collect_fit_data <- function(model, extract,
                             ignore = c("fit", "optimized") %-% extract) {
  #"model" is a list of lists:
  #"fit": data.frame posterior samples, one per fit
  #"optimized": max-likelihood fit, one per fit
  #and other, arbitrary identifying variables, one per fit
  model[ignore] <- NULL
  identifiers <- names(model) %-% ignore %-% extract
  chain(
    model,
    (Map %<<% .)(list), #zip
    lapply(., mkchain(
      c(., .[[extract]]), put(.[[extract]], NULL),
      data.frame %()% .,
      mutate(., .n=1:nrow(.)), #remember the sample id
      melt(id.vars=c(".n", identifiers)),
      unfactor
      )),
    rbind.fill %()% .
    )
}

violinPlot <- function(samples, optimized) {
  #shows the posterior distributions over each aprameter for each subject
  (ggplot(samples)
   + aes(subject, value, fill=model_name, color=model_name)
   + facet_wrap("variable", scales="free")
   + geom_violin(size=0.1, alpha=0.5, position="identity")
   + geom_point(data=optimized, shape=4)
   )
}

crossPlot <- function(samples, optimized, filter, subsample=200) {
  #shows posterior distribusions two at a time
  if (!missing(filter)) {
    samples <- match_df(samples, filter)
    optimized <- match_df(samples, optimized)
  }

  matchnames <- names(samples) %-% "value"

  if (subsample < Inf) {
    samples <- ddply(
      samples, matchnames %-% n,
      function(x)
      if (nrow(x) > subsample) x[sample(nrow(x), subsample),] else x)
  }

  xvar <- yvar <- unique(samples$variable)

  plotdata <- chain(
    unique(samples$variable),
    expand.grid(xvar=xvar, yvar=yvar),
    subset(xvar < yvar),
    merge(samples, by.x="xvar", by.y="variable"),
    merge(samples,
          by.x=revalue(matchnames, c("variable"="yvar")), by.y=matchnames))

  (ggplot(plotdata)
   )
}

main <- function(...) {
  plotStanFits(c(...))
}

run_as_command()
