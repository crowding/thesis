require(plyr)

latex.percent <- function(x) {
  # a latex-compatible formatter for percentages.
  x <- reshape::round_any(x, ggplot2::precision(x)/10)
  paste(ggplot2::comma(x*100), "\\%", sep="")
}

latex.format <- function (x, ...) {
  ## format a number to look pretty in latex.
  ## based on 'scinot' in the 'emdbook' package
  y <- strsplit(as.character(format(x, ...)), 
                "e")[[1]]
  y[1] <- gsub("^0+", "", y[1])
  y[2] <- ifelse(length(grep("^\\+", y[2])) > 0, gsub("^\\+0+", 
                                                      "", y[2]), gsub("^-0+", "-", y[2]))
  ifelse(is.na(y[2]), y[1],
         ifelse(y[1] == "1",
                paste("10^{", y[2], "}", sep = ""),
                paste(y[1], "\\times 10^{", y[2], "}", sep = "")))
}

pretty.pval <- function(p, alpha=0.0001) {
  ##format scientific notion as TeX
  if (p < alpha) {
    paste("\\ensuremath{p < ", latex.format(alpha), "}", sep="")
  } else {
    paste("\\ensuremath{p = ", latex.format(p, digits=2), "}", sep="")
  }
}

paste.comma <- function(..., comma=",", sep=" ", sep.and = "and", 
                        oxford.comma=TRUE, collapse=FALSE,
                        collapse.sep=sep, collapse.and=sep.and, collapse.comma=comma) {
  args <- lapply(list(...), as.character)
  if (length(args) > 1)
    args[[length(args)]] <- paste(sep.and, args[[length(args)]], sep=sep)
  if (length(args) > 2) {
    if (oxford.comma)
      args[1:(length(args)-1)] <- lapply(args[1:(length(args)-1)], paste, comma, sep="")
    else
      args[1:(length(args)-2)] <- lapply(args[1:(length(args)-2)], paste, comma, sep="")
  }
  pasted <- do.call(paste, c(args, sep=sep))
  if (collapse)
    splat(paste.comma)(pasted, sep=collapse.sep, sep.and=collapse.and, comma=collapse.comma, oxford.comma=oxford.comma)
  else pasted
}
