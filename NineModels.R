suppressPackageStartupMessages({
  library(ptools)
  library(plyr)
  library(rstan)
  library(knitr)
  library(R.cache)
  source("library.R")
})
set_cppo('fast')

#displacement
scenarios <- list(d=list(
    soft_local=list(
        displacement_parameter='
          real <lower=spacing_sensitivity*min(frac_spacing),
                upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;',
        displacement_var = '',
        displacement_computation = '
          displacement_factor <-
             -blur * displacement_sensitivity * spacing_sensitivity
              * log(  exp(- frac_spacing[n] / blur)
                    + exp(- max_sensitivity / spacing_sensitivity
                            * 2 * pi() / blur));',
        displacement_R_computation = alist(
            displacement_factor = (
                -blur * displacement_sensitivity * spacing_sensitivity
                * log(   exp(- frac_spacing / blur)
                      + exp(- max_sensitivity / spacing_sensitivity
                            * 2 * pi / blur))))),
    soft_global=list(
        displacement_parameter='
          real <lower=spacing_sensitivity*min(frac_spacing),
                upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;',
        displacement_var = '
          real inverse_number;',
        displacement_computation = '
          inverse_number <- 2*pi() / target_number_shown[n];
          displacement_factor <-
             -blur * displacement_sensitivity * spacing_sensitivity
              * log(  exp(- inverse_number / blur)
                    + exp(- max_sensitivity / spacing_sensitivity
                            * 2 * pi() / blur));',
        displacement_R_computation = alist(
            inverse_number = 2*pi/target_number_shown,
            displacement_factor = (
                -blur * displacement_sensitivity * spacing_sensitivity
                * log(  exp(- inverse_number / blur)
                      + exp(- max_sensitivity /
                              spacing_sensitivity * 2 * pi / blur)))))
    ,
    soft_windowed=list(
        displacement_parameter='
          real <lower=spacing_sensitivity*min(frac_spacing),
                upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;
          real <lower=0, upper=2*pi()> displacement_field;',
        displacement_var = '
          real inverse_number; // equivalent spacing',
        displacement_computation = '
          //maximum of actual spacing and field width divided by targets shown
          inverse_number <- blur*log(
              exp(displacement_field / target_number_shown[n] / blur)
            + exp(2*pi() / target_number_shown[n] / blur));
          displacement_factor <-
            -blur * displacement_sensitivity * spacing_sensitivity
             * log(  exp(- inverse_number / blur)
                   + exp(- max_sensitivity / spacing_sensitivity
                           * 2 * pi() / blur));',
        displacement_R_computation = alist(
            inverse_number = blur*log(
                  exp(displacement_field / target_number_shown / blur)
                + exp(2*pi() / target_number_shown / blur)),
            displacement_factor = (
                -blur * spacing_sensitivity
                * log(  exp(- inverse_number / blur)
                      + exp(- max_sensitivity / spacing_sensitivity
                            * 2 * pi / blur)))))
    ),
                  #carrier
    c=list(
        global=list(
            carrier_parameter = '',
            carrier_var = '',
            carrier_computation=
                'carrier_factor <- target_number_all[n] * carrier_sensitivity;',
            carrier_R_computation=alist(
                carrier_factor = target_number_all * carrier_sensitivity)),
        local=list(
            carrier_parameter = '',
            carrier_var = '',
            carrier_computation =
                'carrier_factor <- 2*pi()*frac_spacing[n] * carrier_sensitivity;',
            carrier_R_computation=alist(
                carrier_factor = 2*pi* frac_spacing * carrier_sensitivity)),
        windowed=list(
            carrier_parameter = 'real<lower=0, upper=2*pi()> carrier_field;',
            carrier_var = '
                real frac_shown;
                real frac_in_carrier_field;
                ',
            carrier_computation = '
                frac_shown <- (target_number_shown[n]+0.0) / target_number_all[n];
                frac_in_carrier_field <- -blur * log(
                    exp(-frac_shown/blur) + exp(-2*pi()*carrier_field/blur));
                carrier_factor <-
                    target_number_all[n] * frac_in_carrier_field * carrier_sensitivity;
                ',
            carrier_R_computation=alist(
                frac_shown = target_number_shown[n] / target_number_all[n],
                frac_in_carrier_field = -blur * log(
                    exp(-frac_shown/blur) + exp(-2*pi*carrier_field/blur)),
                carrier_factor =
                    target_number_all * frac_in_carrier_field * carrier_sensitivity))
        ))

modelTemplate <- '
data {
  int<lower=0> N;
  int target_number_all[N];
  int target_number_shown[N];
  real displacement[N];
  real content[N];
  real spacing[N];
  real eccentricity[N];
  int n_cw[N];
  int n_obs[N];
  real blur;
  real lapse_limit;
}
transformed data {
  real frac_spacing[N];
  for (n in 1:N)
      frac_spacing[n] <- 2*pi()./target_number_all[n];
}
parameters {
  real <lower=0, upper=lapse_limit> lapse;
  real bias;
  real displacement_sensitivity;
  real carrier_sensitivity;
  real <lower=0>spacing_sensitivity;
  {{displacement_parameter}}
  {{carrier_parameter}}
  real repulsion;
  real nonlinearity;
}
model {
  {{displacement_var}}
  {{carrier_var}}
  real carrier_factor;
  real displacement_factor;
  real link_displacement;
  real link_repulsion;
  real link_summation;
  real link;
  for (n in 1:N) {
    {{displacement_computation}}
    {{carrier_computation}}
    link_displacement <- displacement[n] * displacement_factor;
    link_repulsion <- (repulsion * content[n]
                       + nonlinearity * (content[n] * abs(content[n])));
    link_summation <- content[n] * carrier_factor;
    link <- bias + link_displacement + link_repulsion + link_summation;
    n_cw[n] ~ binomial( n_obs[n],
      inv_logit( link ) .* (1-lapse) + lapse/2);
  }
}'

predictorTemplate <- quote(stan_predict <- mkchain[., coefs](
    mutate(frac_spacing = 2*pi/target_number_all)
    , with(coefs, summarize(
        .
        , ...(displacement_computation)
        , ...(carrier_computation)
        , link_displacement = (displacement_factor * displacement)
        , link_repulsion = (repulsion * content
                            + nonlinearity * (content * abs(content)))
        , link_content = (content * carrier_sensitivity)
        , link = bias + link_displacement + link_repulsion + link_summation
        , response = plogis(link) * (1-lapse) + lapse/2))))

otherFunctions <- quote({
  filter_data <- mkchain(
      subset(exp_type %in% c("spacing", "content", "numdensity"))
      , match_df(., subset(count(., "subject"), freq>2000), on="subject")
      )

  splits <- c("subject", "content",
              "displacement",
              "target_number_all",
              "target_number_shown",
              "spacing", "eccentricity", "bias")

  relevant <- splits %v% c("n_cw", "n_obs")

  blur <- 0.2
  lapse_limit <- 0.05

  stan_format <- mkchain(
      subset(select=relevant),
      as.list,
      factorify,
      put(names(.), gsub('\\.', '_', names(.))),
      within({
        N <- length(subject_ix)
        blur <- blur
        lapse_limit <- lapse_limit
      }))
})

makeModelEnv <- function(selection=lapply(scenarios, mkchain(names, .[[1]])),
                         scenarios=parent.env(environment())$scenarios,
                         envir = new.env(parent=globalenv())) {
  #substitute bits and pieces into the
  envir <- as.environment(envir)
  with(envir, {
    source("stanFunctions.R", local=TRUE)
  })
  substituting.env <- new.env(parent=envir)
  mapply(names(selection), selection,
         FUN=function(name, value) {
           assignments <- scenarios[[name]][[value]]
           mapply(names(assignments), assignments,
                  FUN = assign,
                  MoreArgs = list(envir=substituting.env)
                  )
         })
  envir$model_code <- do.call(knit_expand,
                              list(text=modelTemplate), envir=substituting.env)
  envir$model_name <- chain(selection, rbind(names(.),.),
                            paste(collapse="_"))
  envir$stan_predict <- eval(do.call(template,
                                     list(predictorTemplate,
                                          substituting.env)),
                             envir)
  eval(otherFunctions, envir)
  envir
}

using <-function(conn, fn) {
  fn(ammoc(conn, on.exit(close(conn))))
}

compileModelEnv <- function(envir, listing_file) {
  filename <- paste0("models/", envir$model_name, ".stan.RData")
  envir$model <- memoizedCall(stan_model,
                              model_name=envir$model_name,
                              model_code=envir$model_code)
  save(file=filename, list=ls(envir), envir=envir)
  message(envir$model_name)
  writeLines(filename, listing_file)
}

main <- function(outfile='NineModels.list') {
  using(file(outfile, open="w"), function(outconn){
    chain(scenarios,
          lapply(names),
          do.call(expand.grid, .),
          as.matrix,
          alply(1, makeModelEnv),
          lapply(compileModelEnv, outconn))
  })
}

run_as_command()

