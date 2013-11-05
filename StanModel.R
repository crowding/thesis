suppressPackageStartupMessages({
  library(rstan)
  library(plyr)
  library(R.cache)
  library(foreach)
  library(digest)
  source("library.R")
  source("stanFunctions.R")
  source("stan_predictor.R")
})

infile <- "data.RData"
grid <- "motion_energy.csv"
modelfile <- "SlopeModel.stan.RData"
outfile <- "SlopeModel.fit.RData"

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

getSampling <- function(object, ..., chains=4, init, parallel=FALSE) {
  if (parallel) {
    seeds <- sample.int(.Machine$integer.max, chains)
    if (!missing("init")) {
      l <- llply(.parallel=parallel, seq_len(chains),
                 function(x) with_time(sampling(
                   object, chain_id=x, chains=1,
                   init=list(init), seed=seeds[[x]], ...)))
    } else {
      l <- llply(.parallel=parallel, seq_len(chains),
                 function(x) with_time(sampling(
                   object, chain_id=x, chains=1,
                   seed=seeds[[x]], ...)))
    }
    timesum <- Reduce(`+`, lapply(l, attr, "time"))
    combined <- sflist2stanfit(l)
    attr(combined, "time") <- timesum
    combined
  } else {
    if (!missing("init")) {
      with_time(sampling(object, chains=chains, init=rep(list(init), chains), ...))
    } else {
      with_time(sampling(object, chains=chains, ...))
    }
  }
}

with_time <- function(x) {
  timing <- system.time(result <- x)
  attr(result, "time") <- timing
  result
}

main <- function(infile="data.RData",
                 grid="motion_energy.csv",
                 modelfile="SlopeModel.stan.RData",
                 outfile="SlopeModel.fit.RData",
                 iter=2000,
                 warmup=iter/2,
                 parallel=TRUE
                 ) {
  #interrupting one of these model fits shoud fail fast...
  #options(error=NULL)
  iter <- as.numeric(iter)
  warmup <- as.numeric(warmup)
  if (is.na(as.logical(parallel))) {
    parallel <- as.numeric(parallel)
    if (!is.na(parallel)) {
      options(cores=parallel)
      parallel <- TRUE
    } else {
      parallel <- FALSE
    }
  } else {
    parallel <- as.logical(parallel)
  }
  e <- load2env(modelfile) #contains as least 'model', 'stan_predict'
  d <- load2env(infile)
  menergy <- read.csv(grid)
  chain(d$data, (e$filter_data)(), (e$format_data)(menergy)) -> e$data
  fits <- ddply_along(.parallel=FALSE,
                      e$data, e$model_split, function(split, chunk) {
    stan_data <- e$stan_format(chunk)
    so <- sample_optimize(e, stan_data, iter=iter,
                          warmup=warmup, parallel=parallel)
    bind[sample=fit, optimize=optimized] <- so
    #now if max sampling is too far below optimized, there's a problem
    #and need to repeat sampling.
    overshoot <- ((max(optimized$lp__) - (max(as.data.frame(fit)$lp__))) /
                  sd(as.data.frame(fit)$lp__))
    overshot <- FALSE
    if (overshoot > 0.25 || sd(as.data.frame(fit)$lp__) > 4) {
      overshot <- TRUE
      old.overshoot <- overshoot
      old.lp <- max(optimized$lp__)
      old.sd <- max(sd(as.data.frame(fit)$lp__))
      cat("retrying due to overshoot", overshoot, "sd", old.sd, "\n")
      startpoint <- optimized[ifelse(is.na(e$pars), TRUE, e$pars)]
      so <- sample_optimize(e, stan_data, iter=iter,
                          warmup=warmup, init=startpoint, parallel=parallel)
      bind[sample=fit, optimize=optimized] <- so
      overshoot <- ((max(optimized$lp__) - (max(as.data.frame(fit)$lp__))) /
                    sd(as.data.frame(fit)$lp__))
      lp <- max(optimized$lp__)
      overshoot <- cbind(model_name=e$model@model_name, split,
                         old.overshoot=old.overshoot,
                         old.lp=old.lp, old.sd=old.sd,
                         lp=lp,
                         sd=sd(as.data.frame(fit)$lp__),
                         overshoot=overshoot,
                         overshot=overshot)
    } else {
      overshoot <- cbind(model_name=e$model@model_name, split,
                         sd=sd(as.data.frame(fit)$lp__),
                         lp=max(optimized$lp__),
                         overshoot=overshoot,
                         overshot=overshot)
    }

    print(split)
    print(fit)
    print(overshoot)
    dump("overshoot", file="overshoots.txt",
         append=file.exists("overshoots.txt"))

    attr(fit, "overshoot") <- overshoot

    quickdf(list(fit = list(fit), optimized=list(optimized)))
  })
  e$fits <- asisify(fits)

  save(file=outfile, envir=e, list=ls(e))
}

sample_optimize <- function(e, stan_data, iter, warmup, init, parallel) {
  #Gather samples
  #make sure "parameters" are set correctly...
  if (!missing(init)) {
    hash <- list("sampling",
                 e$model@model_name, e$model@model_code, e$model@model_cpp,
                 data=stan_data, warmup=warmup, iter=iter,
                 pars=e$pars, init=init, version=packageVersion("rstan"))
  } else {
    hash <- list("sampling",
                 e$model@model_name, e$model@model_code, e$model@model_cpp,
                 data=stan_data, warmup=warmup, iter=iter,
                 pars=e$pars, version=packageVersion("rstan"))
  }

  cat("hash", digest(hash), "\n")
  ((fit <- loadCache(hash)) %||%
   ammoc(fit <- getSampling(e$model, data=stan_data,
                            warmup=warmup, iter=iter,
                            parallel=parallel,
                            pars=e$pars, init=init),
         saveCache(fit, hash)))

  # then also get the max-likelihood parameters.
  # And the Hessian at max likelihood? Nah.
  # Note that we don't restrict the params here -- max-ll includes
  # posterior predictions
  startpoint <- normalize_coefs(maxll(as.data.frame(fit)))
  #optimizing occasionally kills everything! So save a temp file.
  hash <- list("optimizing", e$model@model_name, e$model@model_code,
               e$model@model_cpp,
               stan_data, init=startpoint, version=packageVersion("rstan"))
  #filename <- paste0("optimizing_", digest(hash), ".RData")
  #save(file=filename, hash, e, fit)
  ((l <- loadCache(hash)) %||%
   ammoc(l <- optimizing(e$model, stan_data, init=startpoint), saveCache(l, hash)))
  #optimized = c(l$par, lp__=list(l$value))
  #unlink(filename)

  list(sample=fit, optimize=c(l$par, lp__=list(l$value)))
}

getOvershoots <- function(file="overshoots.txt") {
  collection <- data.frame()
  makeActiveBinding(env=environment(), "overshoot",
                    function(value) {
                      collection <<- rbind.fill(collection, value)
                    })
  source(file, local=TRUE)
  collection
}

run_as_command()
