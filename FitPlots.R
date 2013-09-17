 suppressPackageStartupMessages({
  library(grid)
  library(plyr)
  library(ggplot2)
  library(data.table)
  library(rstan)
  library(vadr)
  library(reshape2)
  source("library.R")
  source("scales.R")
  source("icons.R")
  source("density_library.R")
  source("ndinterp.R")
  theme_set(theme_bw())
  #theme_update(panel.border=element_blank())
})

infile <- "Hemifield.fit.RData"
grid <- "motion_energy.csv"
plotfile <- "Hemifield.plots.pdf"

main <- function(infile="Hemifield.fit.RData",
                 grid="motion_energy.csv",
                 plotfile="Hemifield.plots.pdf",
                 plots=c("numdensity", "pfun", "sampling")) {
  e <- load2env(infile)
  class(e) <- c("stanenv", class(e))
  #inject a motion-energy interpolator if necessary. An interpolator modifies
  #a data frame. This one interpolates over "norm_diff"
  if (!exists("interpolator", e, inherits=FALSE)) {
    message("building interpolator for motion energy...")
    menergy <- read.csv(grid)
    e$interpolator <- do.call(interpolator, c(list(menergy), e$interpolator.args))
  }
  if (!missing(plotfile)) {
    cairo_pdf(plotfile, onefile=TRUE)
    on.exit(dev.off, add=TRUE)
  }
  for (i in plots) {
    switch(i, numdensity={
      if (any(with(e$data, target_number_shown < target_number_all))) {
        message("plotting number/density data...")
        numdensity_plot(e)
      }
    }, pfun={
      message("plotting psychometric functions...")
      fullCirclePlots(e, fold=TRUE)
    }, contour={
      message("plotting contours...")
      contourPlots(e, fold=TRUE, menergy=menergy)
    }, sampling={
      message("plotting sampling...")
      crossPlots(e)
    }, stop("unknown plot type"))
  }
}

slice <- function(object, selector) UseMethod("slice")

slice.stanenv <- function(object, selector) {
  structure(list(model = object, selector = selector,
                 data = match_df(object$data, selector)),
            class="stanslice")
}

predict.stanslice <- function(object,
                              newdata=object$data,
                              type="response") {
  newdata <- cbind(newdata, unrowname(object$selector))
  predict(predictable(object$model), newdata=newdata)
}

contourPlots <- function(object, fold=fold, menergy) {
  contour_env <- new.env(parent=environment())

  #actually source locally
  source <- function(..., local=ignored) {do.call(base::source, envir=parent.frame(),
                                   list(..., local=TRUE))}
  with(contour_env, source("contours.R"))
  a_ply(object$fits, 1, function(row) {
    bind[fit=, optimized=, ...=group] <- row
    print(group)
    group$full_circle <- 1
    slice <- slice(object, group)
    contour_env$plot_contours(model=slice, subject=group$subject,
                              motion.energy=menergy, fold=TRUE)
  })
}

numdensity_plot <- function(object) {
  segment.data <- extract_segment(object$data, subjects=object$fits$subject,
                                  splits=object$splits)
  folded.segment.data <- extract_segment(object$data, subjects=object$fits$subject,
                                         splits=object$splits,
                                         fold=TRUE, spindle=TRUE, collapse=TRUE)
  x <- predictable(object)
  predictions <- cbind(segment.data, predict(x, segment.data, se.fit=TRUE))

  folded.predictions <- mutilate.predictions(
      predictions, fold=TRUE, spindle=TRUE, collapse=TRUE)

  print(plot.spacing %+% segment.data + facet_wrap(~ subject + displacement + content)
        + density_prediction_layers(predictions, connect="number")
        + errorbars(facet="label", segment.data))
  print(plot.spacing %+% folded.segment.data + facet_wrap(~ label)
        + density_prediction_layers(folded.predictions, connect="number")
        + errorbars(facet="label", folded.segment.data))
  NULL
}

prediction_dataset <-
    function(fit, data=fit$data,
             fold=FALSE,
             ordinate = "displacement",
             menergy = data,
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
      newdata, object$model_split,
      .progress=if(!is.null(summary)) "text" else "none",
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

invoke <- function(data, f, ...) f %()% data

interpolator <- function(
    menergy,
    interpolating=c("displacement", "content", "spacing", "extent", "fullcircle"),
    interpolated=c("norm_diff", "energy_diff"),
    matched=c()
    ) {
  menergy <- chain(menergy, add_energies, subset(as.logical(grid)),
                   mutate(extent = 2*pi*(target_number_shown/target_number_all),
                          spacing = 2*pi*eccentricity/target_number_all,
                          fullcircle = target_number_all == target_number_shown),
                   data.table)
  x <- function(data) {
    data <- chain(data,
                  mutate(extent=2*pi*(target_number_shown/target_number_all),
                         spacing = 2*pi*eccentricity/target_number_all,
                         fullcircle = target_number_all == target_number_shown))

    #we have to interpolate any columns that don't have a match in the grid.
    count_unmatched_values <-
        mkchain( #which column value are not matched in the grid
                list(data[[.]], menergy[[.]]),
                lapply(unique),
                invoke(match), is.na, sum)
    interpolate.by <-
        chain(interpolating,
              . %in% names(data),
              interpolating[.], #column names
              vapply(count_unmatched_values, 0),
              .[.>0],
              sort, rev, names)
    match.by <- interpolating %-% interpolate.by

        #assert that all
    #i.e.
    setkeyv(menergy, match.by)
    chunker <- function(chunk) {
      interp <- sapply(
          interpolated, USE.NAMES=TRUE, simplify=FALSE,
          function(interp.var) {
            if (length(interpolate.by) >= 1) {
              interp.over <- menergy[unique(chunk[match.by])]
              grid <- acast(interp.over,
                            lapply(interpolate.by, mkchain(as.name, as.quoted)),
                            value.var=interp.var, fun.aggregate=mean)
              names(dimnames(grid)) <- interpolate.by
              if (any(dim(grid) == 1)) {
                stop("Columns ",
                     paste(names(dimnames(grid))[dim(grid)==1], collapse=", "),
                     " have only one entry in grid")
              }
              ref <- lapply(dimnames(grid), as.numeric)
              interp <- interp.nd(chunk[names(dimnames(grid))], grid, ref, rule=2)
            } else {
              interp <- menergy[chunk[match.by]][[interp.var]]
            }
          })
      cbind(chunk, quickdf(interp))
    }
    ddply(data, match.by, chunker)
  }
  x
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
  data <- subset(data, target_number_shown == target_number_all)

  switch(style, bubble = {
    plotdata <- mkrates(refold(data, fold=fold), splits=splits)
  }, binned = {
    plotdata <- bin_along_resid(model, data,
                                "response", splits, "displacement", fold=fold)
  })

  preddata <- prediction_dataset(plotdata, fit=fits)
  predictions <- chain(preddata,
                       predict(fits, ., type="response"))
  ## predictions <- chain(preddata,
                       ## predict(fits, ., type="response", se.fit=TRUE))
  if (fold) {
    ## predictions2 <- chain(preddata, fold_trials(fold=TRUE),
    ##                       predict(fits, ., type="response", se.fit=TRUE))
    predictions2 <- chain(preddata, fold_trials(fold=TRUE),
                          predict(fits, ., type="response"))
    #predictions$fit <- (predictions$fit + (1-predictions2$fit))/2
    predictions <- (predictions + (1-predictions2))/2
    #predictions$se.fit <- sqrt((predictions$se.fit^2 + predictions2$se.fit^2)/2)
  }
  #predictions <- cbind(preddata, predictions)
  predictions <- cbind(preddata, fit=predictions)
  #error bars don't add much here and are horrible to compute
  predictions$se.fit <- 0

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
    df <- df[order(df$.order),]
    dfse <- predict(object$stanenv, newdata, summary=colwise_se_frame)
    dfse <- dfse[order(df$.order),]
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

condition_warn <- function(samples, columns = c("value"),
                           include = function(data) data$parameters != "lp__",
                           flag = function(x) NULL) {
  cond <- lapply(samples[columns], function(value)
                 ((!is.finite(value)) | (abs(value) > 1E3)))
  includes <- include(samples)
  cond <- Reduce(`|`, cond) & includes
  if (any(cond & includes)) {
    message("Infinite or nan or large values!")
    message(deparse(substitute(samples)), ":", "\n")
    print(unique(samples[cond,"parameters", drop=FALSE]))
  }
  if(sum(!cond) == 0) stop("no no no no no")
  samples[!cond, , drop=FALSE]
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
  flag <- FALSE;
  setFlag <- function(x) flag <<- TRUE;
  plotdata <- chain(
    fit,
    short=as.matrix,
    .[sample(min(subsample, nrow(.))),],
    melt,
    long=condition_warn(flag=setFlag),
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
  if(flag) browser();
  (ggplot(plotdata)
   + aes(x=x, y=y)
   + facet_grid(facet.y ~ facet.x, scales="free")
   + geom_point(size=0.15)
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

try_command <- mkchain(strsplit(" +"), unlist, .[-1:-2], main %()% .)

run_as_command()
