source("library.R")
library(rstan)

splits <- c("subject", "content",
            "displacement",
            "target_number_all",
            "target_number_shown",
            "spacing", "eccentricity", "bias")

relevant <- splits %v% c("n_cw", "n_obs")

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
modelfile <- "SlopeModel.stan.RData"
outfile <- "SlopeModel.fit.RData"

main <- function(infile="data.RData",
                 grid="motion_energy.csv",
                 modelfile="SlopeModel.stan.RData",
                 outfile="SlopeModel.fit.RData",
                 iter=2000
                 ) {

  load(modelfile) #contains as least 'model', 'stan_predict'

  local({
    e <- load2env(infile)
    menergy <- read.csv(grid)
    chain(e$data, filter_data, format_data(menergy))
  }) -> data

  fits <- ddply_along(data, model_split, function(split, chunk) {
    print(split)
    stan_data <- stan_format(chunk)
    fit <- sampling(model, data=stan_data, warmup=iter/2, iter=iter)
    print(split)
    print(fit)
    quickdf(list(fit = list(fit)))
  })
  fits <- asisify(fits)

  save(file=outfile, list=ls() %-% names(match.call()))
  invisible(fits)
}

run_as_command()
