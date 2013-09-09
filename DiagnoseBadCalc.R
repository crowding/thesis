suppressPackageStartupMessages({
  library(vadr)
  library(rstan)
  library(plyr)
  library(R.cache)
  source("library.R")
  source("stanFunctions.R")
})

#what we want to do here is...

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

dofit <- function(infile="data.RData",
                 grid="motion_energy.csv",
                 modelfile,
                 outfile,
                 filter_data=e$filter_data,
                 iter=2000
                  ) {
  print(c(infile=infile, modelfile=modelfile, outfile=outfile))
  #interrupting one of these model fits shoud fail fast...
  options(error=NULL)
  iter <- as.numeric(iter)
  e <- load2env(modelfile) #contains as least 'model', 'stan_predict'
  d <- load2env(infile)
  menergy <- read.csv(grid)
  chain(d$data, filter_data, (e$format_data)(menergy)) -> e$data
  fits <- ddply_along(e$data, e$model_split, function(split, chunk) {
    stan_data <- e$stan_format(chunk)
    hash <- list("sampling",
                 e$model@model_name, e$model@model_code, e$model@model_cpp,
                 data=stan_data, warmup=iter/2, iter=iter)
    fit <- (loadCache(hash) %||%
            sampling(e$model, data=stan_data,
                          warmup=iter/2, iter=iter))
    saveCache(fit, hash)
    print(split)
    print(fit)

    # then also get the max-likelihood parameters. And the Hessian? Nah.
    startpoint <- normalize_coefs(maxll(as.data.frame(fit)))
    hash <- list("optimizing", e$model@model_name, e$model@model_code, e$model@model_cpp,
                 stan_data, init=startpoint)
    l = (loadCache(hash) %||% optimizing(e$model, stan_data, init=startpoint))
    saveCache(l, hash)
    optimized = c(l$par, lp__=list(l$value))

    quickdf(list(fit = list(fit), optimized=list(optimized)))
  })
  e$fits <- asisify(fits)

  save(file=outfile, envir=e, list=ls(e))
}

filter_without <- mkchain(
    subset(exp_type %in% c("spacing", "content")),
    match_df(., subset(count(., "subject"), freq>2000), on="subject"))

filter_with <- mkchain(
    subset(exp_type %in% c("spacing", "content", "numdensity")),
    match_df(., subset(count(., "subject"), freq>2000), on="subject"))

main <- function() {
  configs <- expand.grid(
      modelfile = c(manual = "SoftMinModel.stan.RData",
                    subst = "models/d_soft_local_c_global.stan.RData"),
      filter_data = list(segment = filter_with, circle = filter_without),
      stringsAsFactors=FALSE)
  configs <- mutate(configs,
                    outfile = paste0("diagnose/fit_",
                                     names(filter_data), "_",
                                     names(modelfile), ".fit.RData"))
  (Map %<<% configs)(dofit)
}

run_as_command()
