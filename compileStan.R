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

    hash <- list("stan_model",
                 model_name=strip_extension(infile),
                 model_code=model_code,
                 stan_version=packageVersion("rstan"))
    ((model <- loadCache(hash)) %||%
     ammoc(model <- stan_model(model_name=strip_extension(infile),
                               model_code=model_code),
           saveCache(model, hash)))
    save(file=outfile, list=ls())
  })
}

run_as_command()
