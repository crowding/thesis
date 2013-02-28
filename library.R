#TODO make this ferret out and install packages, if necessary.
#maybe put the package tarball in here? or maybe devtools is nice enough.
suppressPackageStartupMessages({
  library(ptools)
  library(R.devices)
  library(plyr)
})

load2env <- function(file, env=new.env()) {
  load(file, envir=env)
  env
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
  rename(data, replacements)
}

refold <- function(data, fold=TRUE) {
  fold.trial <- with(data, fold & ((content < 0)
                                   | (content == 0 & displacement < 0)))
  response <- NA
  p <- NA
   chain(data,
        mutate(content = ifelse(fold.trial, -content, content),
               displacement = ifelse(fold.trial, -displacement, displacement),
               response = ifelse(fold.trial, !response, response),
               p = ifelse(fold.trial, 1-p, p),
               fit = if(exists("fit")) ifelse(fold.trial, 1-fit, fit)))
}

mkrates <- function(data,
                    splits=c("displacement", "content",
                      "spacing", "subject", "exp_type")) {
  counter <- function(s) summarize(s,
    n = length(response), p = mean(response),
    n_cw = sum(response), n_ccw = sum(!response))
  nullcounter <- function(s) mutate(s, n=I(c()), p=I(c()), n_cw=I(c()), n_ccw=I(c()))
  chain(data
        , if(empty(.)) nullcounter(.) else ddply(.,splits, counter)
        , arrange(desc(n)))
}

seq_range <- function(range, ...) seq(from=range[[1]], to=range[[2]], ...)

`%++%` <- function(x, y) paste(x, y, sep="")

`%-%` <- setdiff

`%v%` <- union

`%^%` <- intersect

`%call%` <- function(x, y) do.call(x, as.list(y), envir=parent.frame())


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
