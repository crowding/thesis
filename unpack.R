#!/bin/Rscript
source("library.R")

# unpack the CSV files into an R data file, with a column "exp_type"
# for the experiment type
main <- function(outfile="data.RData",
                 filenames=dir(".", "*.series_trials.csv"), ...) {
  filenames <- c(filenames, ...)
  types <- sub('.*?([a-zA_Z]*)_series_trials.*', '\\1', filenames)
  data <- adply(data.frame(f=I(filenames), t=types), 1,
                with, cbind(read.csv(f), exp_type=t))
  save(data, file=outfile)
}

run_as_command()
