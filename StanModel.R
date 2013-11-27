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
  if (!missing("init")) {
    if (is.character(init)) {
      NULL
    } else if (length(names(init))==0) {
      if (is.na(init) || length(init) == 0)
          init <- "random"
    } else if (!is.list(init) || !is.list(init[[1]])) {
      init <- rep(list(init), chains)
    }
  } else {
    init <- "random"
  }
  if (parallel) {
    if (!is.list(init)) {
      init <- rep(init, chains)
    } else if (!is.list(init[[1]])) {
      init <- rep(list(init), chains)
    }
    seeds <- sample.int(.Machine$integer.max, chains)
    l <- llply(.parallel=parallel, seq_len(chains),
               function(x) with_time(sampling(
                 object, chain_id=x, chains=1,
                 init=init[x], seed=seeds[[x]], ...)))
    timesum <- Reduce(`+`, lapply(l, attr, "time"))
    combined <- sflist2stanfit(l)
    attr(combined, "time") <- timesum
    combined
  } else {
    with_time(sampling(object, chains=chains, init=init, ...))
  }
}

findInit <- function(object, stan_data, parallel,
                     seeds=sample.int(.Machine$integer.max, 10), pars, ...) {
  # start by optimizing a few, random times
  if (missing(pars)) pars <- TRUE
  inits <- lapply(seeds, function(x) optimizing(object, seed=x, ...))
  i <- which.max(vapply(inits, `[[`, 0, "value"))
  print(pars)
  as.list(inits[[i]]$par[pars])
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
    #    if(interactive() && split$subject=="nj") debug(sample_optimize)
    init <- if (exists("init", e)) e$init() else "random"
    so <- sample_optimize(e, stan_data, iter=iter, init=init,
                          warmup=warmup, parallel=parallel)
    bind[sample=fit, optimize=optimized] <- so
    #now if max sampling is too far below optimized, there's a problem
    #and need to repeat sampling.
    overshoot <- ((max(optimized$lp__) - (max(as.data.frame(fit)$lp__))) /
                  sd(as.data.frame(fit)$lp__))
    overshot <- FALSE
    max_overshoot <-
        if (packageVersion("rstan") >= package_version("2.0.0"))
            5 else 0.3
    if (overshoot > max_overshoot || sd(as.data.frame(fit)$lp__) > 4) {
      overshot <- TRUE
      old.overshoot <- overshoot
      old.lp <- max(optimized$lp__)
      old.sd <- max(sd(as.data.frame(fit)$lp__))
      cat("retrying due to overshoot", overshoot, "sd", old.sd, "\n")
      ## print(interactive())
      print(split$subject)
      ## if(interactive() && split$subject=="nj") {
      ##   print("what the hell")
      ##   debug(sample_optimize)
      ## }
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
    #if(interactive() && split$subject=="nj") undebug(sample_optimize)
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
  if (missing(init)) {
    #we'll have to find an init point
    init_hash <- list("findInit", e$model@model_name,
                      e$model@model_code, e$model@model_cpp, data=stan_data,
                      pars=e$pars, version=packageVersion("rstan"))
    #f <- findCache(init_hash)
    #if (!is.null(f)) unlink(f)
    (init <- loadCache(init_hash) %||%
     ammoc(init <- findInit(e$model, data=stan_data,
                            parallel=parallel, pars=e$pars),
           saveCache(init, init_hash)))
    print(init)
  }

  hash <- list("sampling",
               e$model@model_name, e$model@model_code, e$model@model_cpp,
               data=stan_data, warmup=warmup, iter=iter, init=init,
               pars=e$pars, version=packageVersion("rstan"))

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

  startpoint <- normalize_coefs(maxll(as.data.frame(fit))[
    if (is.null(e$pars) || all(is.na(e$pars))) TRUE else e$pars])
  #optimizing occasionally kills everything! So save a temp file.
  hash <- list("optimizing", e$model@model_name, e$model@model_code,
               e$model@model_cpp,
               stan_data, init=startpoint, version=packageVersion("rstan"))
  filename <- paste0("optimizing_", digest(hash), ".RData")
  save(file=filename, hash, e, fit)
  ((l <- loadCache(hash)) %||%
   ammoc(l <- optimizing(e$model, stan_data, init=startpoint, samplefile = "tmp.dump"),
         saveCache(l, hash)))
  #optimized = c(l$par, lp__=list(l$value))
  unlink(filename)

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
