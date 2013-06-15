source("library.R")
library(rstan)

`%+%` <- union

splits <- c("subject", "content",
            "displacement",
            "target_number_all",
            "target_number_shown",
            "spacing", "eccentricity", "bias")

relevant <- splits %+% c("n_cw", "n_obs")

model_split <- "subject"

filter_data <- mkchain(
    subset(exp_type %in% c("spacing", "content"))
    )

format_data <- mkchain[., energy](
    do.rename(folding=FALSE)
    , mutate(bias=1)
    , match_df(., subset(count(., "subject"), freq>2000), on="subject")
    , attach_motion_energy(energy)
    , add_energies
    , mutate(data, displacement=wrap(displacement,spacing))
    , mkrates(splits)
    , recast_data
    )

stan_format <- mkchain[., prediction=.](
                              rbind.fill(mutate(., isdata=TRUE),
                                         mutate(prediction, isdata=FALSE)),
                              subset(select=c(relevant, "isdata")),
                              as.list,
                              factorify,
                              put(names(.), gsub('\\.', '_', names(.))),
                              within({
                                N <- length(subject_ix)
                                N_data <- sum(isdata)
                                i_data <- which(isdata)
                                i_predict <- which(!isdata)
                                rm(isdata)
                              }))

factorify <- mkchain(
    lapply(function(x)
           switch(class(x),
                  character=factorify_col(factor(x)),
                  factor=factorify_col(x),
                  list(x))),
    unlist(recursive=FALSE)
    )

factorify_col <- function(x) {
  list(n = length(x),
       #levels = levels(x),
       ix = as.numeric(x))
}

infile <- "data.RData"
grid <- "motion_energy.csv"
modelfile <- "StanModel.RData"
outfile <- "StanModel.RData"

main <- function(infile="data.RData",
                 grid="motion_energy.csv",
                 modelfile="StanModel_model.RData",
                 outfile="StanModel_fit.RData"
                 ) {

  bind[model=model] <- as.list(load2env(modelfile))

  local({
    e <- load2env(infile)
    menergy <- read.csv(grid)
    chain(e$data, filter_data, format_data(menergy))
  }) -> data

  chunk <- subset(data, subject=="pbm")

  fits <- ddply_along(data, model_split, function(split, chunk) {

    print(split)
    stan_data <- stan_format(chunk)
    fit <- sampling(model, data=stan_data, warmup=1000, iter=2000)
    print(split)
    print(fit)
    quickdf(list(fit = list(fit)))
  })
  fits <- asisify(fits)

  filter_data <- filter_data
  format_data <- format_data
  stan_format <- stan_format
  save(file=outfile, list=ls())
  invisible(fits)
}

run_as_command()
