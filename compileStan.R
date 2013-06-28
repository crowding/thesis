suppressPackageStartupMessages({
  library(ptools)
  library(rstan)
  source("library.R")
})
set_cppo('fast')

infile <- "SlopeModel.stan.R"
outfile <- "SlopeModel.stan.RData"

main <- function(infile, outfile) {
  local({
    source("stanFunctions.R", local=TRUE)
    source(infile, local=TRUE)
    model <- stan_model(model_name=strip_extension(infile), model_code=model_code)
    save(file=outfile, list=ls())
  })
}

run_as_command()
