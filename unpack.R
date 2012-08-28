#!/bin/Rscript
source("modeling/library.R")

# unpack the CSV files into an R data file.
main <- function(outfile="data.Rd", ...) {
  filenames <- c(...)
  if (length(filenames) == 0) {
    filenames <- dir(".", "*.series.trials.csv");
  }

  data <- structure(  lapply(filenames, read.csv)
                    , names=grep(  '([a-zA_Z]*)_series_trials.*'
                                 , filenames, value=TRUE))
  data$all <- rbind(data$content, data$spacing)

  save(data, file=outfile)
}

run_as_command()
