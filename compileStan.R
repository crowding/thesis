suppressPackageStartupMessages({
  library(vadr)
  library(plyr)
  library(rstan)
  library(R.cache)
  source("library.R")
})
set_cppo('fast')

infile <- "SlopeModel.stan.R"
outfile <- "SlopeModel.stan.RData"

main <- function(infile, outfile) {
  local({
    base::source("stanFunctions.R", local=TRUE)
    base::source(infile, local=TRUE)

    model <- memoizedCall(stan_model,
                          model_name=strip_extension(infile),
                          model_code=model_code)
    save(file=outfile, list=ls())
  })
}

run_as_command()
