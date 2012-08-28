#!/bin/Rscript
source(ifelse(file.exists("modeling"), "modeling/library.R", "library.R"))

# unpack the CSV files into an R data file, with a column "exp_type"
# for the experiment type
main <- function(outfile="data.Rd", ...) {
  filenames <- c(...)
  if (length(filenames) == 0) {
    filenames <- dir(".", "*.series.trials.csv");
  }
  types <- sub('.*?([a-zA_Z]*)_series_trials.*', '\\1', filenames)
  data <- mapply(function(f, t) cbind(read.csv(f), exp_type=t)
                 , filenames, types, SIMPLIFY=FALSE)
  data <- do.call(rbind, data)
  
  save(data, file=outfile)
}

run_as_command()
