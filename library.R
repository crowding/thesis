suppressPackageStartupMessages({
  #library(R.devices)
  library(stringr)
  library(psyphy)
  library(ptools)
})

## function to build a nonlinear term as used by the gnm package.
## called like;
##
## nonlinearTerm(parameterTerms...)(predictorTerms...)(expression)(predictor...,
## start) where "expression" is written in terms of "parameterTerms"
## and "predictorTerms", and "predictor" nominates data terms
## corresponding to predictor terms. "start" is a named vector of
## starting values corresponding to the parameters.
##
##   displacementTerm <- (nonlinearTerm(cs, beta_dx)(spacing, displacement)
##                     ((2 - 2/(1+exp(-cs/spacing))) * beta_dx * displacement))
## nonlinearTerm()
nonlinearTerm <- function(..., start=NULL) {
  predictors <- quote_args(...)
  predictors[is.missing(predictors)] <- list(1)
  function(...) {
    variables <- quote_args(...)
    #
    function(expr) {
      expr <- substitute(expr)
      eval(template(structure(
        class="nonlin",
        function(.=...(variables), start, trace=FALSE) {
          ...( lapply(names(variables),
                      function(x) template( .(as.name(x))
                                           <- substitute(.(as.name(x))))))
          c(list( predictors = alist(...(predictors))
                 , variables = list(...(lapply(names(variables), as.name)))
                 , term = function(predLabels, varLabels) {
                   t = as.list(parse(text=c(predLabels, varLabels)))
                   names(t) <- .(c(names(predictors), names(variables)))
                   x <- paste(deparse(substitute(.(expr), t)), collapse="")
                   if(trace) {
                     print(formula)
                     print(x)
                   }
                   x
                 }),
            if (missing(start)) NULL else list(start=function(x) {
              if(trace) {
                print(start[names(x)])
              }
              start[names(x)]
            }))
        })))
    }
  }
}

binom.fam <- binomial(link=logit.2asym(g=0.025, lam=0.025))

drop_recursive <- function(df) df[!vapply(df, is.recursive, FALSE)]

predict_from_model_frame <- function(
  models, newdata, fold=TRUE, spindle=TRUE, collapse=FALSE) {
  ##take a data frame with a list of models, and the variables to
  ##match by, produce predictions for the folding data.
  newdata_missing <- missing(newdata)
  chain(models,
        adply(1, function(row) {
          bind[model=bind[model], ...=group] <- as.list(row)
          if (newdata_missing) {
            newdata <- predict_from_model(model)
          } else {
            predict_from_model(model,
                               match_df(newdata, quickdf(group),
                                        on = names(newdata) %^% names(group)))
          }
        }),
        drop_recursive,
        if (any(fold, spindle, collapse)) {
          mutilate.predictions(., fold=fold, spindle=spindle, collapse=collapse)
        } else .)
}

predict_from_model <- function(
  model, newdata=model$data, se.fit=TRUE) {
  #chunk prediction because predict.gnm does something odd
  newdata$.chunk <- floor(seq_len(nrow(newdata))/1000)
  ddply(newdata, ".chunk",
        function(chunk) {
          if (se.fit && tryCatch(error=function(x) FALSE, {
                pred <- predict(model, newdata=chunk,
                                type="response", se.fit=TRUE)
                TRUE
              })) {
            cbind(chunk, pred[1:2], model=NA)
          } else {
            cbind(chunk
                  , fit = predict(model, newdata=chunk, type="response"))
          }
        })
}

##We'll be modeling raw data, but plotting folded/spindled. Here's a
##function that "re-folds" the predictions so that they can be plotted
##on a folded plot.
mutilate.predictions <-
  function(pred,
           fold=abs(diff(range(sign(pred$content)))) > 1,
           spindle=length(unique(pred$side)) > 1,
           collapse=FALSE) {
    columns <- c(as.quoted(splits),
                 if (spindle) NULL else as.quoted("side"))
    chain(pred,
          refold(fold),
          ddply_keeping_unique_cols(
            columns, summarize,
            fit = mean(fit), se.fit = sqrt(sum(se.fit^2))),
          if(collapse) collapse(.) else .,
          labeler)
  }

folding_predict <- function(model, newdata=model$data, type="link", fold=FALSE, ...) {
  pred <- predict(model, newdata=newdata, type=type, ...)
  if(fold) {
    pred2 <- predict(model, newdata=fold_trials(newdata, TRUE), type=type, ...)
    switch(class(pred),
           list = Map(a=pred, b=pred2, n = names(pred), f=function(a, b, n) {
             switch(n,
                    se.fit = (a + b) / 2 / sqrt(2),
                    fit=switch(type,
                      response=(a + (1-b)) / 2,
                      terms = (a - b) / 2),
                    a+b/2)
           }),
           switch(type,
             response=(pred+(1-pred2))/2,
             (pred - pred2) / 2))
  } else {
    pred
  }
}

collapse <- function(data) {
  #collapses different sides and direction contents together (as for
  #most subjects in this experiment these don't matter.)
  data <- chain(data, subset(abs(content) > 0 & displacement/sign(content) < 0.45))
  args <- dots(
    data,
    segment.config.vars %v% segment.experiment.vars %-% c("displacement", "content"),
    summarize)

  if ("n_obs" %in% names(data)) {
    if ("fit" %in% names(data)) {
      args <- args %__% dots(fit = mean(fit*n_obs)/sum(n_obs),
                             se.fit = sqrt(sum((se.fit^2)*n_obs)/sum(n_obs)))
    }
    args <- args %__% dots(n_ccw = sum(n_ccw), n_cw = sum(n_cw),
                           n_obs = sum(n_obs), p = n_cw/n_obs)
  } else {
    if ("fit" %in% names(data)) {
      args <- args %__% dots(fit = mean(fit), se.fit = sqrt(mean(se.fit^2)))
    } else {
      stop("collapse has to have count data or prediction data")
    }
  }
  ddply_keeping_unique_cols %()% args
}

match_to <- function(x, from) {
  vapply(x, function(x) from[which.min(abs(log(x) - log(from)))], 0)
}

#label function for each facet
labeler <- function(data) {
  l <- if ("exp_type" %in% names(data)) {
    unlist(use.names=FALSE, dlply(data, "exp_type", with, {
      monodisp <- exists("displacement") && length(unique(displacement))==1
      monocon <- exists("content") && length(unique(content))==1
      switch(
        as.character(exp_type[[1]]),
        numdensity = {
          paste0(
            if (!monodisp && !monocon) "Observer " else "",
            sprintf("%s", toupper(subject) ),
            if (monodisp) paste0(" d=", format(displacement, digits=2)) else "",
            if (monocon) paste0(" C=", format(content, digits=2)) else ""
            )
        },
        content = sprintf("Content, observer %s", toupper(subject)),
        spacing = sprintf("Spacing, observer %s", toupper(subject)))
    }))
  } else {
    with(data,
         if (exists("full_circle") && length(unique(full_circle)) > 1) {
           sprintf("Observer %s", toupper(subject))
         } else {
           sprintf("Observer %s", toupper(subject))
         })
  }
  cbind(data, label=l)
}

zip <- function(l, collate=c) {
  do.call("mapply", c(list(FUN=collate, SIMPLIFY=FALSE), l))
}

#' Try to bin values coming from staircase data into fewer
#' values.
#'
#' Uses individually-counted data, not trial data. (because everyone
#' hates my bubble scatter plots.)  Doing this honestly, i.e. most
#' faithfully representing the difference between model and data,
#' requires doing this with reference to the model (so that31
#' differences between prediction and model leave the same
#' residuals). Adds new columns "n" and "p".
#' @title
#' @param data Data on the individual trial level
#' @param responsevar
#' @param predictvar
#' @param split Column names that define each separate psychometric function.
#' @param along Column names that define the abscisssa of each psychometric function.
#' @param n_bins The number of bins to split evenly into
#' @param binsize The bin size to use. Either this or n_bins but be defined.
#' @return A data frame with the "split" and "along" columns plus "n" and "p"
#' @author Peter Meilstrup
bin_along <- function(data, responsevar, predictvar,
                      split, along, bins=6) {
  ddply(data, split, function(chunk) {
    a <- chunk[[along]]
    chunk$.bin <- floor((order(a) - 1)/length(a) * bins)
    ddply(chunk, .(.bin), function(x) {
      names <- c("p", along, "nobs")
      quickdf(structure(list(mean(x[[responsevar]]), mean(x[[along]]),
                             nrow(x)),
                        names=names))
    })
  })
}

dlply_along <- function(.data, .variables, .fun=NULL, ...
                        , .progress = "none", .inform = FALSE, .drop = TRUE
                        , .parallel = FALSE, .paropts = NULL) {
  split <- plyr:::splitter_d(.data, as.quoted(.variables), drop=.drop)
  rows <- plyr:::splitter_a(attr(split, "split_labels"), 1)
  llply(  .data=seq_len(length(split))
        , .fun=function(i) .fun(rows[[i]], split[[i]], ...)
        , .progress = .progress, .inform = .inform
        , .parallel = .parallel, .paropts = .paropts)
}

motion_energy_model <- function(model, energy) {
  model$motion.energy <- energy
  class(model) <- union("motion_energy_model", class(model))
  model
}

predict.motion_energy_model <- function(object, newdata=model$data, ...) {
  newdata <- attach_motion_energy(newdata, object$motion.energy)
  NextMethod("predict", object=object, newdata=newdata, ...)
}

bin_along_resid <-
  function(model, data, responsevar, split, along,
           bins=6, restrict, fold=FALSE)
  UseMethod("bin_along_resid")

bin_along_resid.motion_energy_model <-
  function(model, data, responsevar, split, along,
           bins=6, restrict, fold=FALSE) {
    if (is.missing(restrict) || is.null(restrict)) {
      restrict <- unique(model$motion.energy[[along]])
    }
    NextMethod("bin_along_resid", model, restrict=restrict)
  }

#bin observations, using an "average" that retains the Pearson
#residual with respect to the model. The intention is to bin towards
#values that can be plotted on top of the model fit and can be used to
#visually assess model fit.
bin_along_resid.default <- function(model, data, responsevar, split, along,
                                    bins=6, restrict, fold=FALSE) {
  # if we are binning "folded"
  missing.restrict <- missing(restrict)
  data <- unmkrates(data)
  data$fit <- predict(model, newdata=data, type="response")
  chain(
    data,
    mutate(., .pred=predict(model, newdata=., type="response")),
    refold(fold),
    mutate(.pred=ifelse(folded, 1-.pred, .pred)),
    ddply_keeping_unique_cols(split, function(chunk) chain(
      chunk,
      mutate(.binvar=eval.quoted(as.quoted(along))[[1]],
             .bin=floor((order(.binvar) - 1)/length(.binvar) * bins)),
      ddply(".bin", function(x) {
        mean_along <- mean(x$.binvar)
        if (!missing.restrict) {
          #actually we want to restrict to values that are observed...
          mean_along <- take_nearest(mean_along, data[[along]])
        }
        l <- structure(list(mean_along), names=along)
        total_obs <- sum(x[[responsevar]])
        total_pred <- sum(x$.pred)
        total_var <- sum(x$.pred * (1-(x$.pred)))
        pearson_resid <- (total_obs - total_pred) / sqrt(total_var)
        quickdf(c(l, list(n_obs = nrow(x), total_obs=total_obs,
                          total_pred=total_pred, total_var=total_var,
                          pearson_resid=pearson_resid)))
      })))) -> binned
  # if we are folding then we have to make a "folded" prediction
  binned$pred <- folding_predict(model, newdata=binned, type="response")
  # this can produce valuce slightly outside [0,1]
  # this should produce identical Pearson residuals.
  binned <- mutate(binned,
                   p = (n_obs*pred + pearson_resid * sqrt(n_obs*(pred)*(1-pred)))/n_obs,
                   new_resid = n_obs*(p - pred)/sqrt(n_obs*(pred)*(1-pred)) )
  #assert that the residual is the same.
  if (any(with(binned, abs(new_resid - pearson_resid) > 0.01)))
    stop("not binning correctly")
  binned
}

##Given some high dimensional data and a grid to reduce it to, do the
##reduction by preserving the residual.
bin_grid_resid <- function(model, grid, data=model$data, coords, fold=FALSE) {
  #Before forcing data onto the grid, calculate residuals from the model
  data$.bin.pred <- folding_predict(fold=fold, model, newdata=data, type="response")
  data <- mutate(data,
                 .bin.total_n = if (exists("n_obs")) n_obs else 1,
                 .bin.total_yes = if (exists("n_obs")) n_cw else response,
                 .bin.total_var = (.bin.pred * (1-.bin.pred)) * .bin.total_n,
                 .bin.total_pred = .bin.pred * .bin.total_n,
                 .bin.total_resid = .bin.total_yes - .bin.total_pred)
  grid$.bin.pred <- folding_predict(fold=fold, model, newdata=grid, type="response")

  #at this point, we can refold if necessary.
  grid <- refold(grid, fold=fold)

  #select grid bins for each data point. Assumes grid is rectilinear(enough).
  coordcols <- paste0(coords, ".coord.", seq_along(coords))
  for (i in seq_along(coords)) {
    data[coordcols[[i]]] <-
      take_nearest(data[[coords[[i]]]], grid[[coords[[i]]]])
  }

  #Map data onto the grid, ensuring that everything matches...
  data$.bin.right_check <- seq_len(nrow(data))
  gridded <- merge(grid, suffix(data, ".data"),
                   by.x = coords,
                   by.y = paste0(coordcols, ".data"))
  stopifnot(data$.bin.right_check.data %in% gridded)
  sumcols <- grep("\\.bin", colnames(gridded), value=TRUE) %-% ".bin.pred"
  gridded <- ddply(gridded, coords, function(chunk) {
    data.frame(pred=chunk$.bin.pred[[1]], lapply(chunk[sumcols], sum))
  })

  gridded <- mutate(
    gridded,
    pearson_resid = (.bin.total_resid.data)/sqrt(.bin.total_var.data),
    n_obs = .bin.total_n.data,
    p = (n_obs*pred
         + pearson_resid * sqrt(n_obs*pred*(1-pred)))/n_obs,
    .new_resid = n_obs*(p - pred)/sqrt(n_obs*(pred)*(1-pred)),
    )
  #assert that the residual is the same.
  stopifnot(with(gridded, abs(.new_resid - pearson_resid) <= 0.01))
  drop_columns(gridded, grep("^\\.", names(gridded), value=TRUE))
}

suffix <- function(df, suffix) {
  names(df) <- paste0(names(df), suffix)
  df
}

shuffle <- function(df) df[sample(seq_len(nrow(df))),]

sample.data.frame <- function(df, howmany=10, replace=FALSE, prob=NULL)
  df[sample(seq_len(nrow(df)), howmany, replace, prob),]

wrap <- function(x, spacing) x - spacing*round(x/spacing)

load2env <- function(file, env=new.env()) {
  load(file, envir=env)
  env
}

unattr <- function(x) `attributes<-`(x, NULL)

idf <- function(df) {
  out <- idf(df);
  class(out) <- union(class(out), "data.frame");
  out
}

do.rename <- function(data, folding=TRUE) {
  replacements <- if (folding) {
    c(folded_direction_content="content",
      folded_displacement="displacement",
      folded_response_with_carrier="response",
      target_spacing="spacing"
      )
  } else {
    c(abs_direction_content="content",
      abs_displacement="displacement",
      abs_response_cw="response",
      target_spacing="spacing"
      )
  }
  chain(rename(data, replacements, warn_missing=FALSE),
        mutate(bias = if (folding) 0 else 1))
}

str_match_matching <- function(...) {
  x <- str_match(...)
  x[!is.na(x[,1]), , drop=FALSE]
}

mutate_when_has <- function(.data, ...) {
  stopifnot(is.data.frame(.data) || is.list(.data) || is.environment(.data))
  cols <- list_quote(...)
  for (col in names(cols)) {
    if (col %in% names(.data))
    .data[[col]] <- eval(cols[[col]], .data, parent.frame())
  }
  .data
}

mutate_when_missing <- function(.data, ...) {
  stopifnot(is.data.frame(.data) || is.list(.data) || is.environment(.data))
  cols <- list_quote(...)
  for(col in names(cols)) {
    if (!col %in% names(.data))
      .data[[col]] <- eval(cols[[col]], .data, parent.frame())
  }
  .data
}

ziprbind <- function(l, collector=rbind.fill) Map %<<% dots(f=collector) %()% l

bound_prob <- function(x) pmax(pmin(x, 1), 0)

library(binom)
binom_se <- function(n, p) sqrt(p*(1-p)/n)
binom_se_upper <- function(n, p)
  binom.confint(bound_prob(p)*n, n, methods="wilson", p=0.682)$upper
binom_se_lower <- function(n, p)
  binom.confint(bound_prob(p)*n, n, methods="wilson", p=0.682)$lower

fold_trials <- function(data, fold.trial) {
  #fold.trial a vector of booleans to say which trials to fold.
  refold_me <- function(trials, fold) {
    #how to "fold" the motion energy calculation and other paired fields
    cw_cols <- str_match_matching(sort(colnames(trials)), "(.*)_cw(.*)")
    ccw_cols <- str_match_matching(sort(colnames(trials)), "(.*)_ccw(.*)")
    diff_cols <- str_match_matching(sort(colnames(trials)), "(.*)_diff(.*)")
    content_cols <- str_match_matching(sort(colnames(trials)), "(.*)content(.*)")
    if (any(cw_cols[,c(2,3)] != ccw_cols[,c(2,3)])) stop("hmm")
    Map(a=cw_cols[,1], b=ccw_cols[,1],
        f=function(a,b) {trials[fold, c(a,b)] <<- trials[fold, c(b,a)]; NULL})
    lapply(diff_cols[,1], function(c) {trials[fold, c] <<- -trials[fold, c]; NULL})
    lapply(content_cols[,1], function(c) {trials[fold, c] <<- -trials[fold, c]; NULL})
    trials
  }
  #p <- NA
  chain(data,
        mutate(folded=(if (exists("folded"))
                       xor(folded, fold.trial) else fold.trial)),
        #n_ccw and n_cw were already taken care of...
        mutate_when_has(
          displacement = (ifelse(folded, -displacement, displacement)),
          response = ifelse(folded, !response, response),
          fit = ifelse(folded, 1-fit, fit),
          p = ifelse(folded, 1-p, p),
          .bin.total_yes = ifelse(
            folded, .bin.total_n - .bin.total_yes, .bin.total_yes),
          .bin.total_pred = ifelse(
            folded, .bin.total_n - .bin.total_pred, .bin.total_pred),
          .bin.total_resid = ifelse(
            folded, -.bin.total_resid, .bin.total_resid),
          .bin.pred = ifelse(
            folded, 1-.bin.pred, .bin.pred)),
        refold_me(fold.trial))
}

refold <- function(data, fold=TRUE) {
  fold.trial <- with(data, fold & ((content < 0)
                                    | (content == 0 & displacement < 0)))
  fold_trials(data, fold.trial)
}

#this is kind of ridiculously slow. Its main purpose it to carry
#around motion energy values which are computed for each unique trial.
ddply_keeping_unique_cols <- function(.data, .columns, .fun, ...) {
  kept.columns <- structure(names(.data), names=names(.data))
  d_ply(.data, .columns, function(df) {
    u <- vapply(kept.columns, function(i)  length(unique(df[[i]])) <= 1, FALSE)
    if (any(!u)) kept.columns <<- kept.columns[names(u[u])]
  })
  .columns <- c(as.quoted(.columns), as.quoted(kept.columns))
  ddply(.data, .columns, .fun, ...)
}

is_rates <- function(data,
                     splits=c("displacement", "content",
                       "spacing", "subject", "exp_type")) {
  if ("n_obs" %in% names(data)) {
    if ((all(data$n_obs == 1))) {
      if (all(count(data, splits)$freq == 1)) {
        return(NA) #stop("can't tell if rate-formatted or not....")
      } else {
        return(FALSE)
      }
    } else {
      if (all(count(data, splits)$freq == 1)) {
        return(TRUE)
      } else {
       return(NA) #stop("can't tell if rate-formatted or not....")
      }
    }
  } else {
    return(FALSE)
  }
}

#convert binary data into counted yes/no data
mkrates <- function(data,
                    splits=c("displacement", "content",
                      "spacing", "subject", "exp_type"), keep_unique=TRUE) {
  is <- is_rates(data, splits)
  if (is.na(is)) {
    data <- unmkrates(data)
  } else if (is) {
    return(data)
  }
  counter <- function(s) summarize(s,
    n_obs = length(response), p = mean(response),
    n_cw = sum(response), n_ccw = sum(!response))
  nullcounter <- function(s) mutate(s, n=I(c()), p=I(c()),
                                    n_cw=I(c()), n_ccw=I(c()))
  if (nrow(data) == 0) {
    nullcounter(data)
  } else {
    chain(data,
          (if (keep_unique) ddply_keeping_unique_cols else ddply)(splits, counter),
          arrange(desc(n_obs)))
  }
}

#undo mkrates, convert counted yes/no data into binary data.
unmkrates <- function(data, keep.count.cols=FALSE, columns=names(data)) {
  is <- is_rates(data, splits)
  if (!is.na(is)) {
    if (!is) {
      return(data)
    }
  }
  rows = inverse.rle(list(
    lengths = if ("n_obs" %in% names(data)) data$n_obs else data$n,
    values = seq_len(nrow(data))))
  responses = inverse.rle(list(
      lengths = rbind(data$n_cw, data$n_ccw),
      values = replicate(nrow(data), c(TRUE,FALSE))))
  data <- data[rows,columns, drop=FALSE]
  data$response <- responses
  if (keep.count.cols) {
    mutate(data, n_cw = ifelse(response, 1, 0), n_ccw = ifelse(response, 0, 1),
           n_obs=1, p=ifelse(response, 1, 0) )
  } else {
    data[c("n_cw", "n_ccw", "n", "p")] <- list()
  }
  data
}

take_nearest <- function(data, candidates) {
  #coerce values to the "nearest" of a set of candidates.
  restrict <- sort(unique(candidates))
  breaks <- c(-Inf, restrict[1:(length(restrict)-1)] + diff(restrict)/2, Inf)
  restrict[findInterval(data, breaks, rightmost.closed=TRUE)]
}

which.pmax <- function(..., na.rm=FALSE) {
  args <- list(...)
  Reduce(x=seq_along(args), init=list(-Inf, 0), function(aia, ib) {
    a <- aia[[1]]; ia <- aia[[2]]
    b <- args[[ib]]
    keep <- a >= b & !(na.rm & is.na(b))
    list(ifelse(keep, a, b), ifelse(keep, ia, ib))
  })[[2]]
}

which.pmin <- function(..., na.rm=FALSE) {
  args <- list(...)
  Reduce(x=seq_along(args), init=list(Inf, 0), function(aia, ib) {
    a <- aia[[1]]; ia <- aia[[2]]
    b <- args[[ib]]
    keep <- a <= b & !(na.rm & is.na(b))
    list(ifelse(keep, a, b), ifelse(keep, ia, ib))
  })[[2]]
}

abs.pmin <- function(..., na.rm=FALSE) {
  args <- list(...)
  select <- which.pmin %<<% c(na.rm=na.rm) %()% lapply(args, abs)
  arr <- cbind %()% lapply(args, unattr)
  arr[cbind(seq_len(nrow(arr)), select)]
}

abs.pmax <- function(..., na.rm=FALSE) {
  args <- list(...)
  select <- which.pmax %<<% c(na.rm=na.rm) %()% lapply(args, abs)
  arr <- cbind %()% lapply(args, unattr)
  arr[cbind(seq_len(nrow(arr)), select)]
}

#compute columns for total "local" motion energy and total "global" motion energy
recast_data <- function(
  data,
  number.factor=if ("number.factor" %in% names(data)) unique(data$number.factor) else 2,
  envelope.factor = if ("envelope.factor" %in% names(data))
    unique(data$number.factor) else number.factor,
  carrier.factor = if ("carrier.factor" %in% names(data))
    unique(data$carrier.factor) else number.factor) {
  #all of content_local, content_global and
  #content/spacing, should be similarly scaled _for full circle
  #stimuli_. this helps keep model fitting on the right starting
  #points.
  #similar goes for spacing and number_shown_a_spacing
  chain(data,
        drop_columns(c("envelope.factor", "carrier.factor")),
        mutate(carrier.factor=carrier.factor,
               envelope.factor=envelope.factor,
               check="foo"),
        mutate_when_missing(
          eccentricity = 20/3,
          target_number_shown = round(2*pi*eccentricity/spacing),
          target_number_all = round(2*pi*eccentricity/spacing),
          content_cw = (content+1/4),
          content_ccw = (1 - content)/4,
          side = factor("all", levels=c("all", "bottom", "left", "right", "top"))),
        mutate(
          content_local = content / spacing,
          content_global = content / 2/pi/eccentricity
          * abs.pmin(
            target_number_all / carrier.factor,
            target_number_shown),
          full_circle = target_number_shown == target_number_all,
          extent = spacing * target_number_shown,
          number_shown_as_spacing =
            2*pi*eccentricity / envelope.factor / abs.pmin(
              target_number_shown,
              target_number_all / envelope.factor)))
}

seq_range <- function(range, ...) seq(from=range[[1]], to=range[[2]], ...)
`%++%` <- function(x, y) paste(x, y, sep="")
`%-%` <- setdiff
`%v%` <- union
`%^%` <- intersect
`%call%` <- function(x, y) do.call(x, as.list(y), envir=parent.frame())
mask.na <- function(x, f) `[<-`(x, !f(x), value=NA)
invoke <- function(x, y) y %()% x

ddply_along <-
  function(.data, .variables, .fun=NULL, ...
           , .progress = "none", .inform = FALSE, .drop = TRUE
           , .parallel = FALSE, .paropts = NULL) {
    split <- plyr:::splitter_d(.data, as.quoted(.variables), drop=.drop)
    rows <- plyr:::splitter_a(attr(split, "split_labels"), 1)
    s <- as.list(seq_len(length(split)))
    attributes(s) <- attributes(split)
    class(s) <- c("split", "list")
    ldply(  .data=s
          , .fun=function(i) .fun(rows[[i]], split[[i]], ...)
          , .progress = .progress, .inform = .inform
          , .parallel = .parallel, .paropts = .paropts)
}

figure <- function(label, ...) {
  #can't import R.devices because of name conflicts
  if (R.devices:::devIsOpen(label)){
    R.devices:::devSet(label)
  } else {
    R.devices:::devNew(...)
    R.devices:::devSetLabel(label=label)
  }
}

strip_extension <- function(filename) {
  sub(  "((.)\\.[^.]*|)$", "\\2", filename)
}

replace_extension <- function(filename, new_extension, append="") {
  sub(  "((.)\\.[^.]*|)$"
      , paste("\\2", append, ".", new_extension, sep="")
      , filename)
}

unique_by <- function(data, columns) {
  dups <- duplicated(data[columns])
  data[!dups,]
}

drop_columns <- function(data, drop) {
  data[colnames(data)[!colnames(data) %in% drop]]
}

assert_finite <- function(x) {
  if(any(!is.finite(x))) stop("Some NaN values found") else x
}

motion_energy_calcs <- dots(
  contrast_diff = assert_finite(
    sqrt(rowSums(data[c(cw_cols)] / target_number_shown))
    - sqrt(rowSums(data[c(ccw_cols)] / target_number_shown))),
  total_e = assert_finite(rowSums(data[c(cw_cols, ccw_cols)])),
  energy_cw = assert_finite(rowSums(data[cw_cols]) / max(total_e)),
  energy_ccw = assert_finite(rowSums(data[ccw_cols]) / max(total_e)),
  energy_diff = assert_finite(energy_cw - energy_ccw),
  energy_total = assert_finite(total_e / max(total_e)),
  norm_diff = assert_finite(energy_diff / energy_total)
  )

add_energies <- function(data, drop=TRUE) {
  cw_cols <- grep("(.*)_cw\\.\\d$", names(data), value=TRUE)
  ccw_cols <- grep("(.*)_ccw\\.\\d$", names(data), value=TRUE)
  channel_cols <- grep("\\.\\d*$", names(data), value=TRUE)
  chain(data,
        mutate(target_number_shown = (
          ifelse(is.na(target_number_shown),
                 target_number_all, target_number_shown))),
        (here(mutate) %<<% motion_energy_calcs)(),
        drop_columns(c(cw_cols, ccw_cols, "total_e")),
        rename(c(abs_displacement="displacement",
                 abs_direction_content="content"), warn_missing=FALSE))
}

attach_motion_energy <- function(trials, motion_energy_frame) {
  chain(motion_energy_frame,
        mutate(
          target_number_shown = ifelse(
            is.na(target_number_shown),
            target_number_all, target_number_shown))
        ) -> menergy

  trials$left.check <- 1:nrow(trials)
  menergy$right.check <- 1:nrow(menergy)

  trials <- drop_columns(trials, names(motion_energy_calcs))
  joined <- join(trials, menergy, type="left",
                  by = intersect(names(trials), names(menergy)))

  if (any(dups <- duplicated(joined$left.check))) {
    #getting a lot of duplicated matches, for
    #some reason?
    warning("Multiple motion-energy matches")
    joined <- joined[!dups,]
  }
  if (any(missed <- is.na(joined$right.check))) {
    stop("Motion energy information not found for all trials")
  }
  drop_columns(joined, c("left.check", "right.check"))
}

#helper function for when merge doesn't retrieve anything
merge_check <- function(x, y, ..., by=intersect(names(x), names(y))) {
  sapply(by, function(col) nrow(merge(x, y, by=by %-% col)))
}

key_check <- function(x, y, ..., by=intersect(names(x), names(y))) {
  sapply(by, function(col) length(unique(intersect(x[[col]], y[[col]]))))
}

show_dup <- function(joined, which_dup = 1) {
  dups <- which(duplicated(joined$left.check))
  subset(joined, left.check %in% left.check[dups[which_dup]])
}

hosmerlem = function(model, newdata=model$data, groups=10) {
  #Hosmer-Lemeshow statistic for goodness of fit, for a binomial
  #response. We could use this in other ways, it's a
  #deviance test essentially.
  newdata$fit <- predict(model, newdata=newdata, type="response")
  #make sure it is binary data and not grouped
  newdata <- unmkrates(newdata, columns=c("fit"))
  bind[y, yhat] <- newdata[c("response", "fit")]
  cutyhat <- floor((order(yhat) - 1)/length(yhat) * groups)
  obs = xtabs(cbind(1 - y, y) ~ cutyhat)
  expect = xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq = sum((obs - expect)^2/expect)
  P = pchisq(-chisq, groups - 1)
  return(list(stat=chisq, dof=groups-1, p=P))
}

subset_deviance <- function(model, subset=TRUE, .enclos=parent.frame()) {
  mask <- eval(substitute(subset), model$data, .enclos)
  y <- model$y[mask]
  mu <- model$fitted.values[mask]
  wts <- model$prior.weights[mask]
  sum(model$family$dev.resids(y, mu, wts))
}

normalize <- function(x) x / mean(x)
standardize <- function(x) (x - mean(x)) / sd(x)
shift_min <- function(x, ...) x - min(x, ...)
normalize_min <- function(x) (x / min(x))

colwise_mutate <- function(...) {
  f <- colwise(...)
  function(df) {
    ff <- f(df)
    df[names(ff)] <- ff
    df
  }
}

captureWarnings <- function(expr, print=TRUE, raise) {
  warnings <- NULL
  i <- 0
  raiseMissing <- missing(raise)
  result <- withCallingHandlers(
    simpleWarning=function(w) {
      i <<- i + 1
      warnings[[i]] <<- w
      if (!raiseMissing && str_detect(w$message, raise)) {
        stop(w)
      }
      if(print) message(w)
      invokeRestart("muffleWarning")
    },
    tryCatch(
      list(result=expr, warnings=warnings, error=NULL),
      error=function(e) {
        if (!raiseMissing && str_detect(e$message, raise)) {
          stop(e)
        }
        if(print) message(e)
        list(result=NULL, warnings=warnings, error=e)
    }))
}

functions.of <- function(..., .envir=parent.frame()) {
  args <- quote_args(...)
  function(...) {
    lapply(quote_args(...), function(x) {
      eval(template(function(.=...(args)) .(x)), .envir)
    })
  }
}

asisify <- function(df) {
  quickdf(lapply(df, function(x) if (is.recursive(x)) I(x) else x))
}

nonrec <- function(df) {
  #only nonrecursive columns
  df[!vapply(df, is.recursive, FALSE)]
}

declare_data <- function(...) {
  #make a data frame from a bunch of list arguments, automatically
  #applyin asIs to recursive columns
  chain(list(...), zip,
        lapply(function(x) if(is.recursive(x)) I(x) else x),
        quickdf)
}

qq <- template

qe <- function(x, envir=parent.frame()) eval(
                    template %()% list(substitute(x), parent.frame()),
                    parent.frame())

factor_in_order <- function(x, filter=function(x) is.character(x) || is.factor(x)) {
  if (filter(x)) factor(x, levels=unique(x)) else x
}

put <- macro(function(assignment, value) {
  assignment_target <- function(x) {
    switch(
      class(x),
      name=x,
      character=as.name(x),
      call = switch(head<- class(x[[1]]),
        name = assignment_target(x[[2]]),
        character = assignment_target(x[[2]]),
        stop("that doesn't look like an assignment target")),
      stop("that doesn't look like an assignment target"))
  }
  target <- assignment_target(assignment)
  template((function() {.(assignment) <- .(value); .(target)})())
})

print.if.nonempty <- function(d) {
  name <- substitute(d)
  if (!empty(d)) {
    cat(name,":", "\n")
    print(d)
  }
}

evalqq <- function(x, envir=parent.frame())
  eval(template(substitute(x), envir))

exploreFun <- function(f, start=rep(1, length(formals(f))),
                       min=start/4, max=start*4, mark=start,
                       envir=parent.frame()) {
  #make a plot of a function's parameters
  anames <- chain(f, match.call(qq(f(...(start)))), .[-1], names)
  fcall <- substitute(f)
  text.all=deparse(as.call(c(list(fcall), put(names(start), anames))))
  title(text.all)
  screens = split.screen(c(1, length(anames)), erase=FALSE)
  tryCatch(
    Map(anames, start, min, max, mark, screens,
        f=function(aname, st, min, max, mark, s) {
          screen(s)

          call <- template(f(...(ifelse(aname==anames,
                                        list(as.name(aname)), start))))
          textcall <- deparse(put(call[[1]], fcall))
          cat(textcall, "\n")
          qe(curve(.(call), .(min), .(max), xname=.(aname), ylab=""),
             envir)
          abline(v=mark, col="gray50", lty="11")
        }
        )
    , finally=close.screen(screens))
}
