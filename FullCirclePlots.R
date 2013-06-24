suppressPackageStartupMessages({
  library(grid)
  library(ggplot2)
  library(rstan)
  library(ptools)
  library(reshape2)
  source("library.R")
})

prediction_dataset <-
    mkchain[
            .,
            fold=FALSE,
            ordinate = "displacement",
            ordinate.values = get_sampling(e$data, menergy, ordinate),
            menergy](
              e=.,
              .$data,
              .[e$splits %v% e$model_split %-% ordinate],
              unique,
              merge(quickdf(structure(list(ordinate.values),
                                      names=ordinate)),
                    by=c()),
              if(fold) refold(.) else .
              )

get_sampling = function(data, grid, ordinate) {
  unique(data[[ordinate]], grid[[ordinate]])
}

maxll <- function(stanenv, split) {
  #select best observed lieklihood from sample
  chain(
    stanenv$fits, merge(split), .$fit[[1]],
    d=as.data.frame, with(which.min(lp__)), d[.,]
    )
}

optimized <- function(stanenv, split, startpoint=maxll(stanenv, split)) {
  #obtain coefficients by gradient descent from the starting point.
  data <- merge(stanenv$data, split)
  model <- merge(stanenv$fits, split)$fit[[1]]
  print(startpoint)
  standata <- stanenv$stan_format(stanenv$format_data(data))
  l = optimizing(model@stanmodel, standata, init=startpoint)
  c(l$par, list(l$lp__))
}
optimized(e, data.frame(subject="pbm"))

column_quantiles <- mkchain(
  lapply(.,
    mkchain(
      quantile(c(0.05, 0.95), names=FALSE),
      put(names(.), c("min", "max")))),
  unlist(recursive=FALSE))

predict.stanenv <- function(stanenv, newdata=stanenv$data,
                            #function that selects some "canonical" coefs
                            selector=optimized,
                            #function that computes summary data along all coefs
                            summary=NULL,
                            samples=250
                            ) {
  ddply_along(newdata, e$model_split, .progress="text", function(split, chunk) {
    fit <- merge(stanenv$fits, split)$fit[[1]]
    coefs <- as.data.frame(fit)
    #usually we select the canonical model by optimization
    selected_coefs <- selector(stanenv, split)
    #note that we expect stan_predict to work on EITHER single coefs
    #and many data or single data and many coefs
    fit <- stanenv$stan_predict(coefs=selected_coefs, chunk)
    if (!is.null(summary)) {
      summary <- adply(chunk, 1, mkchain(e$stan_predict(coefs=coefs), column_quantiles)
                       , .progress="text")
      cbind(chunk, fit=fit, summary)
    } else {
      cbind(chunk, fit=fit)
    }
  })
}

infile <- "SlopeModel.fit.RData"
grid <- "motion_energy.csv"
plotfile <- "SlopeModel.plots.pdf"

main <- function(infile="SlopeModel.fit.RData", grid="motion_energy.csv",
                 plotfile="SlopeModel.plots.pdf") {
  e <- load2env(infile);
  class(e) <- c("stanenv", class(e));
  menergy <- read.csv(grid)
  pred_dataset <- prediction_dataset(e, menergy=menergy, fold=TRUE)
  cairo_pdf(plotfile, onefile=TRUE)
  crossPlots(e)
#  fullCirclePlots(e)
  on.exit(dev.off)
}

makeTitle <- function(...) {
  bind[fit=fit, ...=group] <- list(...)
  paste(fit@model_name,
        paste(names(group), toupper(as.character(group)),
              collapse=", "),
        sep=", ")
}

fullCirclePlots <- function(e, ...) {
  extraArgs <- dots(...)
  (Map %<<% e$fits)(function(...) {
    bind[fit=fit, ...=group] <- list(...)
    title <- makeTitle(...)
    chunk <- match_df(fits$data, group)
    fullCirclePlot(predictable, data=chunk)
  })
}

predictable <- function(stanenv, data) {
  put(class(list(stanenv=stanenv, data=data)), "predictable")
}

predict.predictable <- function (object, newdata=object$data, ...) {
  predict(object$stanenv, newdata, ...)
}

fullCirclePlot <- function(fits, data, pred.dataset,
                           style=c("bubble", "binned"), fold=FALSE) {
    style <- match.arg(style)
  subdata <- match_df(data,
                      data.frame(subject=subject, stringsAsFactors=FALSE),
                      on="subject")
  switch(style, bubble = {
    plotdata <- mkrates(refold(subdata, fold=fold))
  }, binned = {
    plotdata <- bin_along_resid(model, subdata,
                                "response", splits, "displacement", fold=fold)
  })
  print((ggplot(plotdata)
         + displacement_scale
         + proportion_scale
         + content_color_scale
         + facet_spacing_experiment
         + plotPredictions(model, data=data, fold=fold, ...)
         + geom_point()
         + (switch(style, bubble=balloon, binned=geom_point()))
         + labs(title = "Data and model fits for subject " %++% subject)
         ))
}

crossPlots <- function(e, ...) {
  extraArgs <- dots(...)
  (Map %<<% e$fits)(function(...) {
    bind[fit=fit, ...=group] <- list(...)
    title <- paste(e$model@model_name,
                   paste(names(group), toupper(as.character(group)),
                         collapse=", "),
                   sep=", ")
    print(title)
    print((crossPlot %<<% extraArgs)(fit=fit, list(labs(title=title))))
    NULL
  })
}

crossPlot <- function(fit, ..., subsample=1000) {
  #plot x-y things...
  plotdata <- chain(
    fit, short=as.matrix,
    .[sample(min(subsample, nrow(.))),],
    long=melt,
    .$parameters, unique,
    expand.grid(facet.x=., facet.y=.),
    subset(as.numeric(facet.x) > as.numeric(facet.y)),
    merge(long,
          by.x="facet.x", by.y="parameters"),
    rename(c(value="x")),
    merge(long,
          by.x=c("facet.y", "iterations"),
          by.y=c("parameters", "iterations")),
    rename(c(value="y")),
    mutate(lp__=short[iterations, "lp__"])
    )
  (ggplot(plotdata)
   + aes(x=x, y=y)
   + facet_grid(facet.y ~ facet.x, scales="free")
   + geom_point(size=0.5, alpha=0.2)
   + theme_bw()
   + theme(strip.text.x=element_text(angle=-90, hjust=1),
           strip.text.y=element_text(angle=0, hjust=0),
           strip.background=element_blank(),
           axis.text.x=element_text(angle=-90, vjust=0.5, hjust=0,
                                    size=rel(0.6), face=2),
           axis.text.y=element_text(angle=0, hjust=1, size=rel(0.5), face=2),
           panel.grid.major=element_line(color=hcl(120, 80, 95)),
           panel.grid.minor=element_blank(),
           panel.border=element_rect(color="gray90"),
           panel.margin=unit(1, "mm"),
           aspect.ratio=1)
   + scale_y_continuous(breaks = curr(pretty, n=3))
   + scale_x_continuous(breaks = curr(pretty, n=3))
   + labs(x="", y="")
   + list(...)
   )
}

run_as_command()
