library(data.table)
source("ndinterp.R")
source("density_library.R")

load_stanfit <- function(file, menergy=read.csv("motion_energy.csv")) {
  e <- load2env(file)
  class(e) <- c("stanenv", class(e))
  #inject a motion-energy interpolator if necessary. An interpolator modifies
  #a data frame. This one interpolates over "norm_diff"
  if (!exists("interpolator", e, inherits=FALSE)) {
    message("building interpolator for motion energy...")
    e$interpolator <- do.call(interpolator, c(list(menergy), e$interpolator.args))
  }
  e
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
  predict(predictable(object$model), newdata=newdata,)
}

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
                           mkchain(object$stan_predict(coefs=subset_coefs),
                                   summary),
                           .progress="text")
          cbind(chunk, summary)
        } else {
          selected_coefs <- selector(object, split)
          fit <- object$stan_predict(coefs=selected_coefs, chunk)
          if (any(is.na(fit))) {
            browser()
          }
          cbind(chunk, fit)
        }
      })
}

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

maxll <- mkchain(as.data.frame, .[which.max(.$lp__)[[1]],])

invoke <- function(data, f, ...) f %()% data

interpolator <- function(
  menergy,
  interpolating=c("displacement", "content", "spacing",
                  "extent", "fullcircle", "contrast"),
  interpolated=c("norm_diff", "energy_diff"),
  matched=c()
  ) {
  menergy <- chain(menergy, add_energies,
                   mutate(extent = 2*pi*(target_number_shown/target_number_all),
                          spacing = 2*pi*eccentricity/target_number_all,
                          fullcircle = target_number_all == target_number_shown,
                          contrast = round_any(content_cw + content_ccw, 0.01)),
#                   subset(content_ccw + content_cw == 0.5),
                   data.table)
  x <- function(data) {
    data <- chain(data,
                  mutate_when_missing(
                      content_cw = (1+content)/4,
                      content_ccw = (1-content)/4),
                  mutate(extent=2*pi*(target_number_shown/target_number_all),
                         spacing = 2*pi*eccentricity/target_number_all,
                         fullcircle = target_number_all == target_number_shown,
                         contrast = round_any(content_cw + content_ccw, 0.01)))

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

    setkeyv(menergy, match.by)
    if (length(interpolate.by) == 0) {
      interp <- cbind(
        data,
        as.data.frame(menergy[data[match.by], mult="first"])[interpolated])
      if (any %()% lapply(interp[interpolated], is.na)) {
        warning("Trying setting contrast to 0.5")
        if ("contrast" %in% match.by) {
          interp <- cbind(data,
                          as.data.frame(menergy[data[match.by %-% "contrast"],
                                                mult="first"])[interpolated])
        }
        if (any %()% lapply(interp[interpolated], is.na)) {
          stop("can't find matches")
        } else {
          return(interp)
        }
      }
      return(interp)
    }
    chunker <- function(chunk) {
      interp <- sapply(
          interpolated, USE.NAMES=TRUE, simplify=FALSE,
          function(interp.var) {
            interp.over <- menergy[unique(chunk[match.by])]
            grid <- acast(interp.over,
                          lapply(interpolate.by, mkchain(as.name, as.quoted)),
                          value.var=interp.var, fun.aggregate=mean, drop=TRUE)
            names(dimnames(grid)) <- interpolate.by
            if (any(dim(grid) == 1)) {
              stop("Columns ",
                   paste(names(dimnames(grid))[dim(grid)==1], collapse=", "),
                   " have only one entry in grid")
            }
            ref <- lapply(dimnames(grid), as.numeric)
            interp <- interp.nd(chunk[names(dimnames(grid))], grid, ref, rule=2)

            if (length(interp) != nrow(chunk)) {
              cat("oops!\n")
              browser()
            }
            if (any(is.na(interp))) {
              cat("wat!\n")
              browser()
            }
            interp
          })
      cbind(chunk, quickdf(interp))
    }
    ddply(data, match.by, chunker)
  }
  x
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
    result = switch(type,
           response=list(fit=df$response, se.fit=dfse$response),
           link=list(fit=df$link, se.fit=dfse$link),
           terms=list(fit=df, se.fit=dfse))
  } else {
    df = predict(object$stanenv, newdata)
    result = switch(type, response=df$response, link=df$link, terms=df)
  }
  result
}

interpolate <- function(...) UseMethod("interpolate")

interpolate.predictable <- function(
    object, newdata=object$data) {
  object$stanenv$interpolator(newdata)
}
