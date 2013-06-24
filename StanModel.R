suppressPackageStartupMessages({
  source("library.R")
  library(rstan)
})

infile <- "data.RData"
grid <- "motion_energy.csv"
modelfile <- "SlopeModel.stan.RData"
outfile <- "SlopeModel.fit.RData"

maxll <- function(fit) mkchain(as.data.frame, .[which.max(lp__)[[1]],])

main <- function(infile="data.RData",
                 grid="motion_energy.csv",
                 modelfile="SlopeModel.stan.RData",
                 outfile="SlopeModel.fit.RData",
                 iter=2000
                 ) {

  e <- load2env(modelfile) #contains as least 'model', 'stan_predict'
  d <- load2env(infile)
  menergy <- read.csv(grid)
  chain(d$data, (e$filter_data)(), (e$format_data)(menergy)) -> e$data

  fits <- ddply_along(e$data, e$model_split, function(split, chunk) {
    print(split)
    stan_data <- e$stan_format(chunk)
    fit <- sampling(e$model, data=stan_data, warmup=iter/2, iter=iter)
    print(split)
    print(fit)
    startpoint <- maxll(fit)
    

    quickdf(list(fit = list(fit)))
  })
  e$fits <- asisify(fits)

  save(file=outfile, envir=e, list=ls(e))
}

run_as_command()
