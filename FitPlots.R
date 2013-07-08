suppressPackageStartupMessages({
  library(grid)
  library(plyr)
  library(ggplot2)
  library(fields)
  library(data.table)
  library(rstan)
  library(ptools)
  library(reshape2)
  source("library.R")
  source("scales.R")
  theme_set(theme_bw())
  #theme_update(panel.border=element_blank())
})

infile <- "OnlyMotionEnergy.fit.RData"
grid <- "motion_energy.csv"
plotfile <- "OnlyMotionEnergy.plots.pdf"

main <- function(infile="OnlyMotionEnergy.fit.RData",
                 grid="motion_energy.csv",
                 plotfile="OnlyMotionEnergy.plots.pdf") {
  e <- load2env(infile)
  class(e) <- c("stanenv", class(e))
  #inject a motion-energy interpolator if necessary. An interpolator modifies
  #a data frame. This one interpolates over "norm_diff"
  if (!exists("interpolator", e, inherits=FALSE)) {
    print("building interpolator for motion energy...")
    menergy <- read.csv(grid)
    e$interpolator <- interpolator(menergy)
  }
  cairo_pdf(plotfile, onefile=TRUE)
  fullCirclePlots(e, fold=TRUE)
  crossPlots(e)
  on.exit(dev.off)
}

prediction_dataset <-
    function(fit, data=fit$data,
             fold=FALSE,
             ordinate = "displacement",
             menergy=data,
             ordinate.values = get_sampling(fit$stanenv$data, menergy, ordinate)){
      chain(data, .[fit$stanenv$splits %v% fit$stanenv$model_split %-% ordinate],
            unique,
            merge(quickdf(structure(list(ordinate.values),
                                    names=ordinate)),
                  by=c()),
            if(fold) refold(.) else .)
    }

get_sampling = function(data, grid=data, ordinate) {
  unique(c(data[[ordinate]], grid[[ordinate]]))
}

maxll <- function(stanenv, split) {
  #select best observed lieklihood from sample
  chain(
    stanenv$fits, merge(split), .$fit[[1]],
    d=as.data.frame, with(which.min(lp__)), d[.,]
    )
}

colwise_se <- mkchain(colwise(sd)(.), put(names(.), paste0(names(.), ".sd")))

colwise_se_frame <- mkchain(colwise(sd)(.))

optimized <- function(stanenv, split, startpoint=maxll(stanenv, split)) {
  if ("optimized" %in% names(stanenv$fits)) {
    return(merge(stanenv$fits, split)$optimized[[1]])
  }
  #obtain coefficients by gradient descent from the starting point.
  data <- merge(stanenv$data, split)
  model <- merge(stanenv$fits, split)$fit[[1]]
  print(startpoint)
  standata <- stanenv$stan_format(stanenv$format_data(data))
  l = optimizing(model@stanmodel, standata, init=startpoint)
  c(l$par, list(l$lp__))
}

column_quantiles <- mkchain(
  lapply(.,
    mkchain(
      quantile(c(0.05, 0.95), names=FALSE),
      put(names(.), c("min", "max")))),
  unlist(recursive=FALSE))

predict.stanenv <- function(object, newdata=object$data,
                            #function that selects some "canonical" coefs
                            selector=optimized,
                            #function that computes summary data along all coefs
                            summary=NULL,
                            samples=Inf
                            ) {
  ddply_along(
    newdata, object$model_split, .progress=if(!is.null(summary)) "text" else "none",
    function(split, chunk) {
      fit <- merge(object$fits, split)$fit[[1]]
      coefs <- as.data.frame(fit)
      #usually we select the canonical model by optimization
      #note that we expect stan_predict to work on EITHER single coefs
      #and many data or single data and many coefs
      if (!is.null(summary)) {
        #we either summarize over all coefs, or select one set of coefs
        #selected_coefs <- selector(object, split)
        ix <- sample(nrow(coefs), min(samples, nrow(coefs)))
        subset_coefs <- coefs[ix,]
        summary <- adply(chunk, 1,
                         mkchain(object$stan_predict(coefs=subset_coefs), summary)
                         , .progress="text")
        cbind(chunk, summary)
      } else {
        selected_coefs <- selector(object, split)
        fit <- object$stan_predict(coefs=selected_coefs, chunk)
        cbind(chunk, fit)
      }
    })
}

interpolator <- function(
    menergy,
    interpolating=c("target_number_shown", "content", "displacement"),
    interpolated=c("norm_diff"),
    matched=c()
    ) {
  menergy <- chain(menergy, subset(grid==TRUE), add_energies,
                   .[c(interpolating, interpolated, matched)],
                   data.table)
  function(data) {
    #a hack -- can't find a 3-d interpolator at the moment.
    #Just ask what plane we're operating over.
    bind[...=outofplane, inplane, inplane[2]] <-
        chain(data, .[interpolating],
              vapply(mkchain(unique, length), 0), sort, names)
    setkeyv(menergy, outofplane)
    ddply(data, outofplane, function(chunk) {
      ldply(interpolated, function(interp.var) {
        interp.over <- menergy[unique(chunk[outofplane])]
        grid <- acast(interp.over,
                      lapply(inplane, mkchain(as.name, as.quoted)),
                      value.var=interp.var)
        interp.obj <- list(x=as.numeric(dimnames(grid)[[1]]),
                           y=as.numeric(dimnames(grid)[[2]]),
                           z=grid)
        loc <- do.call("cbind", chunk[inplane])
        interp <- interp.surface(interp.obj, loc)
        put(chunk[[interp.var]], interp)
      })
    })
  }
}

makeTitle <- function(...) {
  bind[fit=fit, optimized=, ...=group] <- list(...)
  paste(fit@model_name,
        paste(names(group), toupper(as.character(group)),
              collapse=", "),
        sep=", ")
}

fullCirclePlots <- function(e, fold=FALSE, ...) {
  extraArgs <- dots(...)
  (Map %<<% e$fits)(function(...) {
    bind[fit=fit, optimized=optimized, ...=group] <- list(...)
    title <- makeTitle(...)
    print(title)
    chunk <- merge(e$data, group)
    fullCirclePlot(predictable(e), data=chunk, group=group,
                   optimized=optimized, splits=e$splits, fold=fold)
    NULL
  })
  NULL
}

fullCirclePlot <- function(fits, data, group, optimized, splits, predictions,
                           style=c("bubble", "binned"), fold=FALSE) {
  style <- match.arg(style)

  switch(style, bubble = {
    plotdata <- mkrates(refold(data, fold=fold), splits=splits)
  }, binned = {
    plotdata <- bin_along_resid(model, data,
                                "response", splits, "displacement", fold=fold)
  })

  preddata <- prediction_dataset(plotdata, fit=fits)
  predictions <- chain(preddata,
                       predict(fits, ., type="response", se.fit=TRUE))
  if (fold) {
    predictions2 <- chain(preddata, fold_trials(fold=TRUE),
                          predict(fits, ., type="response", se.fit=TRUE))
    predictions$fit <- (predictions$fit + (1-predictions2$fit))/2
    predictions$se.fit <- (predictions$se.fit + predictions2$se.fit)/2
  }
  predictions <- cbind(preddata, predictions)

  #if folding need to average folded and unfolded predictions
  print((ggplot(plotdata)
         + displacement_scale
         + proportion_scale
         + content_color_scale
         + facet_spacing_rows
         + prediction_layer(predictions)
         + switch(style, bubble=balloon, binned=geom_point())
         + labs(title = "Data and model fits for observer "
                %++% toupper(group$subject))
         + theme(aspect.ratio=0.25)
         ))
}

predictable <- function(stanenv, data=stanenv$data) {
  structure(list(stanenv=stanenv, data=data), class="predictable")
}

#pretend to act something like "glm.predict"
predict.predictable <- function (
  object, newdata=object$data,
  se.fit=FALSE, type=c("response", "terms", "link"), ...) {
  type <- match.arg(type)
  newdata$.order <- seq_len(nrow(newdata))
  if (exists("interpolator", object$stanenv, inherits=FALSE)) {
    newdata <- object$stanenv$interpolator(newdata)
  }
  newdata <- newdata[order(newdata$.order),]
  if (se.fit) {
    df <- predict(object$stanenv, newdata, select=optimized)
    dfse <- predict(object$stanenv, newdata, summary=colwise_se_frame)
    switch(type,
           response=list(fit=df$response, se.fit=dfse$response),
           link=list(fit=df$link, se.fit=dfse$link),
           terms=list(fit=df, se.fit=dfse))
  } else {
    df = predict(object$stanenv, newdata)
    switch(type, response=df$response, link=df$link, terms=df)
  }
}

interpolate <- function(...) UseMethod("interpolate")

interpolate.predictable <- function(
    object, newdata=object$data) {
  object$stanenv$interpolator(newdata)
}

crossPlots <- function(e, ...) {
  extraArgs <- dots(...)
  (Map %<<% e$fits)(function(...) {
    bind[fit=fit, ...=group, optimized=] <- list(...)
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
