suppressPackageStartupMessages({
  library(rstan)
  library(plyr)
  source("library.R")
  source("stanFunctions.R")
})

infile <- "data.RData"
grid <- "motion_energy.csv"
modelfile <- "SlopeModel.stan.RData"
outfile <- "SlopeModel.fit.RData"

maxll <- mkchain(as.data.frame, .[which.max(.$lp__)[[1]],])

#take a vector of coefficients with the weird names that Stan produces
#and put it into the list data structure that seems to be implied
normalize_coefs <- function(coefs) {
  e <- new.env()
  vars <- lapply(names(coefs), mkchain(parse(text=.), .[[1]]))
  mapply(coefs, vars, FUN=function(c, v) {
    if (length(v) > 1) {
      n <- as.character(v[[2]])
      v[[2]] <- call("$", quote(e), v[[2]])
      if (!n %in% ls(e)) e[[n]] <- c()
    } else {
      v <- call("$", quote(e), v)
    }
    eval(bquote(.(v) <- .(c)))
  })
  as.list(e)
}

main <- function(infile="data.RData",
                 grid="motion_energy.csv",
                 modelfile="SlopeModel.stan.RData",
                 outfile="SlopeModel.fit.RData",
                 iter=2000
                 ) {
  iter <- as.numeric(iter)
  e <- load2env(modelfile) #contains as least 'model', 'stan_predict'
  d <- load2env(infile)
  menergy <- read.csv(grid)
  chain(d$data, (e$filter_data)(), (e$format_data)(menergy)) -> e$data

  fits <- ddply_along(e$data, e$model_split, function(split, chunk) {
    stan_data <- e$stan_format(chunk)
    fit <- sampling(e$model, data=stan_data, warmup=iter/2, iter=iter)
    print(split)
    print(fit)

    # then also get the max-likelihood parameters. And the Hessian? Nah.
    startpoint <- normalize_coefs(maxll(as.data.frame(fit)))
    l = optimizing(e$model, stan_data, init=startpoint)
    optimized = c(l$par, lp__=list(l$value))

    quickdf(list(fit = list(fit), optimized=list(optimized)))
  })
  e$fits <- asisify(fits)

  save(file=outfile, envir=e, list=ls(e))
}

run_as_command()
