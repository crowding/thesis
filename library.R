suppressPackageStartupMessages({
  #library(R.devices)
  library(stringr)
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
        function(.=...(variables), start) {
          ...( lapply(names(variables),
                      function(x) template( .(as.name(x))
                                           <- substitute(.(as.name(x))))))
          c(list( predictors = alist(...(predictors))
                 , variables = list(...(lapply(names(variables), as.name)))
                 , term = function(predLabels, varLabels) {
                   t = as.list(parse(text=c(predLabels, varLabels)))
                   names(t) <- .(c(names(predictors), names(variables)))
                   deparse(substitute(.(expr), t))
                 }),
            if (missing(start)) NULL else list(start=function(x) start[names(x)]))
        })))
    }
  }
}

predict_from_model_frame <- function(models, newdata,
                                     fold=TRUE, spindle=TRUE, collapse=FALSE) {
  ##take a data frame with a list of models, and the variables to
  ##match by, produce predictions for the folding data.
  newdata_missing <- missing(newdata)
  chain(models,
        adply(1, function(row) {
          bind[model=bind[model], ...=group] <- as.list(row)
          if (newdata_missing) {
            predict_from_model(model)
          } else {
            predict_from_model(model,
                               match_df(newdata, quickdf(group),
                                        on = names(newdata) %^% names(group)))
          }
        }),
        if (any(fold, spindle, collapse)) {
          mutilate.predictions(., fold=fold, spindle=spindle, collapse=collapse)
        } else .)
}

predict_from_model <- function(model, newdata=model$data) {
  #chunk prediction because predict.gnm does something odd
  newdata$.chunk <- floor(seq_len(nrow(newdata))/1000)
  ddply(newdata, ".chunk",
        function(chunk) {
          pred <- predict(model, newdata=chunk, type="response", se.fit=TRUE)
          cbind(chunk, pred[1:2], model=NA)
        })
}

match_to <- function(x, from) {
  vapply(x, function(x) from[which.min(abs(log(x) - log(from)))], 0)
}

collapse <- function(data) {
  #collapses different sides and direction contents together (as for
  #most subjects in this experiment these don't matter.)
  args <- dots(
    chain(data, subset(abs(content) > 0 & displacement/sign(content) < 0.45)),
    segment.config.vars %v% segment.experiment.vars %-% c("displacement", "content"),
    summarize)

  if ("n" %in% names(data)) {
    if ("fit" %in% names(data)) {
      args <- args %__% dots(fit = mean(fit*n)/sum(n),
                             se.fit = sqrt(sum((se.fit^2)*n)/sum(n)))
    }
    args <- args %__% dots(n_ccw = sum(n_ccw), n_cw = sum(n_cw),
                           n = sum(n), p = n_cw/n)
  } else {
    if ("fit" %in% names(data))
      args <- args %__% dots(fit = mean(fit), se.fit = sqrt(mean(se.fit^2)))
  }
  ddply_keeping_unique_cols %()% args
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
          ddply_keeping_unique_cols(columns, summarize,
                fit = mean(fit), se.fit = sqrt(sum(se.fit^2))),
          if(collapse) collapse(.) else .,
          labeler)
  }

#label function for each facet
labeler <- function(data) {
  ddply(data, "exp_type", function(data) {
    switch(data$exp_type[[1]],
           numdensity={
             if ("displacement" %in% names(data)) {
               mutate(data, label = sprintf("%s d=%s C=%s",
                              toupper(subject),
                              format(displacement, digits = 2),
                              format(content, digits = 2)))
             } else {
               mutate(data, label = sprintf("%s", toupper(subject)))
             }
           },
           content=mutate(data, label="Content %s", toupper(subject)),
           spacing=mutate(data, label="Spacing %s", tpupper(subject)))
  })
}

#' Try to bin values coming from staircase data into fewer
#' values.
#'
#' Uses individually-counted data, not trial data. (because everyone
#' hates my bubble scatter plots.)  Doing this honestly, i.e. most
#' faithfully representing the difference between model and data,
#' requires doing this with reference to the model (so that
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
      names <- c("p", along, "n")
      quickdf(structure(list(mean(x[[responsevar]]), mean(x[[along]]),
                             nrow(x)),
                        names=names))
    })
  })
}

#bin observations, using an "average" that retains the Pearson
#residual with respect to the model. Depending on your perspective
#this may be a more "honest" depiction of the model's fit to the data.
bin_along_resid <- function(model, data, responsevar, split, along, bins=6) {
  data$fit <- predict(model, data, type="response")
  binned <- ddply(data, split, function(chunk) {
    a <- chunk[[along]]
    chunk$.bin <- floor((order(a) - 1)/length(a) * bins)
    chunk$.pred <- predict(model, newdata=chunk, type="response")
    chunk <- ddply(chunk, .(.bin), function(x) {
      l <- structure(list(mean(x[[along]])), names=along)
      total_obs <- sum(x[[responsevar]])
      total_pred <- sum(x$.pred)
      total_var <- sum(x$.pred * (1-(x$.pred)))
      pearson_resid <- (total_obs - total_pred) / sqrt(total_var)
      quickdf(c(l, list(n = nrow(x), total_obs=total_obs, total_pred=total_pred,
                        total_var = total_var, pearson_resid=pearson_resid)))
    })
    #we want a value for X that leads to the same Pearson residuals as
    #we have observed.
  })
  binned$pred <- predict(model, binned, type="response")
  #this can produce valuce slightly outside [0,1]
  binned <- mutate(binned,
                   p = (n*pred + pearson_resid * sqrt(n*(pred)*(1-pred)))/n,
                   new_resid = n*(p - pred)/sqrt(n*(pred)*(1-pred)) )
  #assert that the residual is the same.
  if (any(with(binned, abs(new_resid - pearson_resid) > 0.01)))
    stop("not binning correctly")
  binned
}



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
  chain(rename(data, replacements),
        mutate(bias = if (folding) 0 else 1))
}

str_match_matching <- function(...) {
  x <- str_match(...)
  x[!is.na(x[,1]), , drop=FALSE]
}

mutate_when_has <- function(data, columns=dots_names(...), ...) {
  if (all(columns %in% names(data))) {
    evalq(function(...) mutate(...), parent.frame())(data, ...)
  } else {
    data
  }
}

binom_se <- function(n, p) sqrt(p*(1-p)/n)

refold <- function(data, fold=TRUE) {
  fold.trial <- with(data, fold & ((content < 0)
                                    | (content == 0 & displacement < 0)))

  refold_me <- function(trials, fold) {
    #how to "fold" the motion energy calculation and other paired fields
    cw_cols <- str_match_matching(sort(colnames(trials)), "(.*)_cw(.*)")
    ccw_cols <- str_match_matching(sort(colnames(trials)), "(.*)_ccw(.*)")
    diff_cols <- str_match_matching(sort(colnames(trials)), "(.*)_diff(.*)")
    if (any(cw_cols[,c(2,3)] != ccw_cols[,c(2,3)])) stop("hmm")
    Map(a=cw_cols[,1], b=ccw_cols[,1],
        f=function(a,b) {trials[fold, c(a,b)] <<- trials[fold, c(b,a)]; NULL})
    lapply(diff_cols, function(c) {trials[fold, c] <<- -trials[fold, c]; NULL})
    trials
  }

  #p <- NA
  chain(data,
        #n_ccw and n_cw were already taken care of...
        mutate_when_has(content = ifelse(fold.trial, -content, content)),
        mutate_when_has(displacement =
                          ifelse(fold.trial, -displacement, displacement)),
        mutate_when_has(response = ifelse(fold.trial, !response, response)),
        mutate_when_has(fit=ifelse(fold.trial, 1-fit, fit)),
        mutate_when_has(p=ifelse(fold.trial, 1-p, p)),
        refold_me(fold.trial))
}

mkrates <- function(data,
                    splits=c("displacement", "content",
                             "spacing", "subject", "exp_type", "bias"),
                    keep=TRUE) {
  counter <- function(s) summarize(s,
    n = length(response), p = mean(response),
    n_cw = sum(response), n_ccw = sum(!response))
  nullcounter <- function(s) mutate(s, n=I(c()), p=I(c()),
                                    n_cw=I(c()), n_ccw=I(c()))
  if (nrow(data) == 0) {
    nullcounter(data)
  } else {
    chain(data,
           (if (keep) ddply_keeping_unique_cols else ddply)(splits, counter),
          arrange(desc(n)))
  }
}

#this is kind of ridiculously slow. Its main purpose it to carry
#around motion energy values which are computed for each unique trial.
ddply_keeping_unique_cols <- function(.data, .columns, .fun, ...) {
  kept.columns <- structure(names(.data), names=names(.data))
  d_ply(.data, .columns, function(df) {
    u <- vapply(kept.columns, function(i)  length(unique(df[[i]])) <= 1, FALSE)
    if (any(!u)) kept.columns <<- kept.columns[names(u[u])]
  })
  summary <- ddply(.data, .columns, .fun, ...)
  kept.columns <- setdiff(kept.columns, names(summary))
  extradata <- ddply(.data, .columns, `[`, 1,  kept.columns)
  cbind(summary, extradata[names(extradata) %-% names(summary)])
}

mkrates <- function(data,
                    splits=c("displacement", "content",
                      "spacing", "subject", "exp_type"), keep_unique=TRUE) {
  counter <- function(s) summarize(s,
    n = length(response), p = mean(response),
    n_cw = sum(response), n_ccw = sum(!response))
  nullcounter <- function(s) mutate(s, n=I(c()), p=I(c()),
                                    n_cw=I(c()), n_ccw=I(c()))
  if (nrow(data) == 0) {
    nullcounter(data)
  } else {
    chain(data,
          (if (keep_unique) ddply_keeping_unique_cols else ddply)(splits, counter),
          arrange(desc(n)))
  }
}

seq_range <- function(range, ...) seq(from=range[[1]], to=range[[2]], ...)
`%++%` <- function(x, y) paste(x, y, sep="")
`%-%` <- setdiff
`%v%` <- union
`%^%` <- intersect
`%call%` <- function(x, y) do.call(x, as.list(y), envir=parent.frame())
mask.na <- function(x, f) `[<-`(x, !f(x), value=NA)

ddply_along <-
  function(.data, .variables, .fun=NULL, ...
           , .progress = "none", .inform = FALSE, .drop = TRUE
           , .parallel = FALSE, .paropts = NULL) {
    split <- plyr:::splitter_d(.data, as.quoted(.variables), drop=.drop)
    rows <- plyr:::splitter_a(attr(split, "split_labels"), 1)
    ldply(  .data=seq_len(length(split))
          , .fun=function(i) .fun(rows[[i]], split[[i]], ...)
          , .progress = .progress, .inform = .inform
          , .parallel = .parallel, .paropts = .paropts)
}

figure <- function(label, ...) {
  #can't use R.devices because of name conflicts

  ## if (devIsOpen(label)){
  ##   devSet(label)
  ## } else {
  ## devNew(...)
  ## devSetLabel(label=label)
  ## }
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

add_energies <- function(data){
  cw_cols <- grep("(.*)_cw\\.\\d$", names(data), value=TRUE)
  ccw_cols <- grep("(.*)_ccw\\.\\d$", names(data), value=TRUE)
  channel_cols <- grep("\\.\\d*$", names(data), value=TRUE)
  chain(data,
        mutate(
          total_e = rowSums(data[c(cw_cols, ccw_cols)]),
          energy_cw = rowSums(data[cw_cols]) / max(total_e),
          energy_ccw = rowSums(data[ccw_cols]) / max(total_e),
          energy_diff = energy_cw - energy_ccw,
          energy_total = total_e / max(total_e),
          norm_diff = energy_diff / energy_total),
        drop_columns(c(cw_cols, ccw_cols, "total_e")),
#        drop_columns(drop_cols),
        rename(c(abs_displacement="displacement",
                 abs_direction_content="content"), warn_missing=FALSE)
        )
}
