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
            if (missing(start)) NULL else list(start=function(x) {print(x); start[names(x)]}))
        })))
    }
  }
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
                             "spacing", "subject", "exp_type", "bias")) {
  counter <- function(s) summarize(s,
    n = length(response), p = mean(response),
    n_cw = sum(response), n_ccw = sum(!response))
  nullcounter <- function(s) mutate(s, n=I(c()), p=I(c()),
                                    n_cw=I(c()), n_ccw=I(c()))
  if (nrow(data) == 0) {
    nullcounter(data)
  } else {
    chain(data,
          ddply_keeping_unique_cols(splits, counter),
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
