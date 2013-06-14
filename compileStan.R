library(ptools)
library(rstan)
source("library.R")
set_cppo('fast')

main <- function(infile, outfile) {
  model <- stan_model(file=infile,
                      model_name=strip_extension(infile),
                      verbose=FALSE)
  save(file=outfile, model)
}

run_as_command()
