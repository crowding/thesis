suppressPackageStartupMessages({
  #library(R.devices)
  library(stringr)
})

load2env <- function(file, env=new.env()) {
  load(file, envir=env)
  env
}

unattr <- function(x) `attributes<-`(x, NULL)

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
  rename(data, replacements)
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

ddply_keeping_unique_cols <- function(.data, .columns, .fun, ...) {
  kept.columns <- structure( rep(TRUE, length(.data))
                            , names=colnames(.data))
  produced.columns <- character(0)
  with_unique_cols <- function(.chunk, .fn, ...) {
    summary <- .fn(.chunk, ...)
    keep <- unlist(lapply(names(.chunk), function(n) {
      if (length(unique(.chunk[[n]])) <= 1) n else NULL
    }))
    kept.columns[names(.chunk)[!names(.chunk) %in% keep]] <<- FALSE
    produced.columns <<- produced.columns %v% names(summary)
    cbind(summary, .chunk[1,keep])
  }
  out <- ddply(.data, .columns, with_unique_cols, .fn=.fun, ...)
  out[produced.columns %v% names(kept.columns[kept.columns])]
}

mkrates <- function(data,
                    splits=c("displacement", "content",
                             "spacing", "subject", "exp_type")) {
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

ddply_keeping_unique_cols <- function(.data, .columns, .fun, ...) {
  kept.columns <- structure( rep(TRUE, length(.data))
                            , names=colnames(.data))
  summary <- ddply(.data, .columns, .fun, ...)
  unique.columns <-
    chain(.data,
          ddply(.columns, colwise(function(x) length(unique(x)) <= 1)),
          `[<-`(unique(as.character(.columns)), value=list()),
          colwise(all)(),
          unlist, .[.], names,
          setdiff(names(summary)))
  extra.data <- ddply(.data, .columns, `[`, 1, unique.columns)
  cbind(summary, extra.data)
}

mkrates <- function(data,
                    splits=c("displacement", "content",
                      "spacing", "subject", "exp_type")) {
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
