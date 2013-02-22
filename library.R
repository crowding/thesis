library(R.devices)
library(stringr)

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

mutate_when_has <- function(data, columns, ...) {
  if (all(columns %in% names(data))) {
    evalq(function(...) mutate(...), parent.frame())(data, ...)
  } else {
    data
  }
}

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
        mutate(content = ifelse(fold.trial, -content, content),
               displacement = ifelse(fold.trial, -displacement, displacement),
               response = ifelse(fold.trial, !response, response)),
        mutate_when_has("fit", fit=ifelse(fold.trial, 1-fit, fit)),
        mutate_when_has("p", p=ifelse(fold.trial, 1-p, p)),
        refold_me(fold.trial))
}


mkrates <- function(data,
                    splits=c("displacement", "content",
                      "spacing", "subject", "exp_type")) {

  kept.columns <- structure( rep(TRUE, length(data))
                           , names=colnames(data))

  counter <- function(s) summarize(s,
    n = length(response), p = mean(response),
    n_cw = sum(response), n_ccw = sum(!response))

  nullcounter <- function(s) mutate(s, n=I(c()), p=I(c()), n_cw=I(c()), n_ccw=I(c()))

  with_unique_cols <- function(summary, data) {
    keep <- unlist(lapply(names(data),
      function(n) if (length(unique(data[[n]])) == 1) n else NULL))
    keep <- keep %-% names(summary)
    kept.columns[names(data)[names(data) %in% keep]] <- FALSE
    cbind(summary, data[1,keep])
  }

  #if there are columns we can propagate, keep them.
  count_and_drop <- function(df)
    chain(df, counter, with_unique_cols(df))

  chain(data
        , if(empty(.)) nullcounter(.) else ddply(.,splits, count_and_drop)
        , arrange(desc(n))
        ) -> out
  out[colnames(out) %-% names(kept.columns)[!kept.columns]]
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
  if (devIsOpen(label)){
    devSet(label)
  } else {
  devNew(...)
  devSetLabel(label=label)
  }
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
