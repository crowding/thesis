suppressPackageStartupMessages({
  library(ptools)
})

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

mkrates <- function(data,
                    splits=c("displacement", "content",
                      "spacing", "subject", "exp_type")) {
  counter <- function(s) summarize(s,
    n = length(response), p = mean(response),
    n_cw = sum(response), n_ccw = sum(!response))
  nullcounter <- function(s) mutate(s, n=I(c()), p=I(c()), n_cw=I(c()), n_ccw=I(c()))
  browser()
  chain(data
        , if(empty(.)) nullcounter(.) else ddply(.,splits, counter)
        , arrange(desc(n)))
}
 
seq_range <- function(range, ...) seq(from=range[[1]], to=range[[2]], ...)

`%++%` <- function(x, y) paste(x, y, sep="")

`%-%` <- setdiff

`%v%` <- union

`%^%` <- intersect
