#!/bin/Rscript
library(ptools)
library(plyr)
source("library.R")

# unpack the CSV files into an R data file, with a column "exp_type"
# for the experiment type
main <- function(outfile="data.RData",
                 filenames=dir(".", "*.series_trials.csv"), ...) {
  filenames <- c(filenames, ...)
  types <- sub('.*?([a-zA_Z]*)_series_trials.*', '\\1', filenames)
  data <- mdply(data.frame(f=I(filenames), t=types),
                function(f, t) cbind(read.csv(f, stringsAsFactors=FALSE), exp_type=t))
  data$f <- NULL
  data$t <- NULL
  save(data, file=outfile)
}

run_as_command()
