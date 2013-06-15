library(ptools)
library(rstan)
source("library.R")
set_cppo('fast')

infile <- "SlopeModel.stan.R"
outfile <- "SlopeModel.stan.RData"

main <- function(infile, outfile) {
  local({
    source(infile, local=TRUE)
    save(file=outfile, list=ls())
  })
}

run_as_command()
