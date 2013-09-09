suppressPackageStartupMessages({
  library(vadr)
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
             -blur * spacing_sensitivity
              * log(  exp(- frac_spacing[n] / blur)
                    + exp(- max_sensitivity / spacing_sensitivity
                            * 2 * pi() / blur));',
        displacement_R_computation = alist(
            displacement_factor <- (
                -blur * spacing_sensitivity
                * log(   exp(- frac_spacing / blur)
                      + exp(- max_sensitivity / spacing_sensitivity
                            * 2 * pi / blur))))),
    soft_local_boosted=list(
        displacement_parameter='
          real <lower=spacing_sensitivity*min(frac_spacing),
                upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;
          real <lower=0,upper=max_sensitivity> min_sensitivity;
          ',
        displacement_var = '',
        displacement_computation = '
          displacement_factor <-
             -blur * spacing_sensitivity
              * log(  exp(- frac_spacing[n] / blur)
                    + exp(- max_sensitivity / spacing_sensitivity
                            * 2 * pi() / blur));
          displacement_factor <- ( (displacement_factor/max_sensitivity)
                                  *(max_sensitivity-min_sensitivity)
                                  + min_sensitivity);
        ',
        displacement_R_computation = alist(
            displacement_factor <- (
                -blur * spacing_sensitivity
                * log(   exp(- frac_spacing / blur)
                      + exp(- max_sensitivity / spacing_sensitivity
                            * 2 * pi / blur))))),
    soft_local_endpoints=list(
        displacement_parameter='
          real <lower=spacing_sensitivity*min(frac_spacing),
                upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;
          real<lower=0, upper=max_sensitivity> endpoint_sensitivity;
',
        displacement_var = '
          real endpoints;
          real raw_spacing_sensitivity;
          real envelope_interior;
          real envelope_norm;
          ',
        displacement_computation = '
          raw_spacing_sensitivity <- -blur * spacing_sensitivity
              * log(  exp(- frac_spacing[n] / blur)
                    + exp(- max_sensitivity / spacing_sensitivity
                            * 2 * pi() / blur));
          if (target_number_shown[n] == target_number_all[n])
            endpoints <- 0;
          else endpoints <- 2;
          envelope_interior <- pow(raw_spacing_sensitivity, 2) * target_number_shown[n];
          envelope_norm <- (  endpoints * endpoint_sensitivity
                              + target_number_shown[n] * raw_spacing_sensitivity);
          displacement_factor <- envelope_interior / envelope_norm;
        ',
        displacement_R_computation = alist(
            endpoints <- ifelse(target_number_shown == target_number_all, 0, 2),
            raw_spacing_sensitivity <- (
                -blur * spacing_sensitivity
                * log(  exp(- frac_spacing / blur)
                      + exp(- max_sensitivity / spacing_sensitivity
                            * 2 * pi / blur))),
            envelope_interior <- raw_spacing_sensitivity^2 * target_number_shown,
            envelope_norm <- (endpoints * endpoint_sensitivity
                              + target_number_shown * raw_spacing_sensitivity),
            displacement_factor <- envelope_interior / envelope_norm
            )),
    soft_global=list(
        displacement_parameter='
          real <lower=spacing_sensitivity*min(frac_spacing),
                upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;',
        displacement_var = '
          real inverse_number;',
        displacement_computation = '
          inverse_number <- 2*pi() / target_number_shown[n];
          displacement_factor <-
             -blur * spacing_sensitivity
              * log(  exp(- inverse_number / blur)
                    + exp(- max_sensitivity / spacing_sensitivity
                            * 2 * pi() / blur));',
        displacement_R_computation = alist(
            inverse_number <- 2*pi/target_number_shown,
            displacement_factor <- (
                -blur *  spacing_sensitivity
                * log(  exp(- inverse_number / blur)
                      + exp(- max_sensitivity /
                            spacing_sensitivity * 2 * pi / blur)))))
    ,
    soft_windowed=list(
        displacement_parameter='
          real <lower=spacing_sensitivity*min(frac_spacing),
                upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;
          real <lower=min(frac_spacing), upper=2*pi()> displacement_field;',
        displacement_var = '
          real inverse_number; // equivalent spacing',
        displacement_computation = '
          //maximum of actual spacing and field width divided by targets shown
          inverse_number <- blur*log(
              exp(displacement_field / target_number_shown[n] / blur)
            + exp(2*pi() / target_number_shown[n] / blur));
          displacement_factor <-
            -blur * spacing_sensitivity
             * log(  exp(- inverse_number / blur)
                   + exp(- max_sensitivity / spacing_sensitivity
                           * 2 * pi() / blur));',
        displacement_R_computation = alist(
            inverse_number <- blur*log(
                exp(displacement_field / target_number_shown / blur)
                + exp(2*pi / target_number_shown / blur)),
            displacement_factor <- (
                -blur * spacing_sensitivity
                * log(  exp(- inverse_number / blur)
                      + exp(- max_sensitivity / spacing_sensitivity
                            * 2 * pi / blur))))))
                  ,
                  #CARRIER
    c=list(
        global=list(
            carrier_parameter = '',
            carrier_var = '',
            carrier_computation=
                'carrier_factor <- target_number_shown[n] * carrier_sensitivity;',
            carrier_R_computation=alist(
                carrier_factor <- target_number_shown * carrier_sensitivity)),
         endpoints=list(
            carrier_parameter = '
                real <lower=0, upper=max_sensitivity> endpoint_carrier_weight;',
            carrier_var = '
                real carrier_weight_inside;
                real carrier_norm;
                real n_endpoints;',
            carrier_computation= '
                if (target_number_shown[n] == target_number_all[n])
                    n_endpoints <- 0;
                else n_endpoints <- 2;
                carrier_weight_inside <- target_number_shown[n] * displacement_factor;
                carrier_norm <-
                   ( target_number_shown[n] * displacement_factor
                    + n_endpoints * endpoint_carrier_weight ) / (n_endpoints + target_number_shown[n]);
                carrier_factor <- carrier_sensitivity *
                                  carrier_weight_inside / carrier_norm;
                ',
            carrier_R_computation=alist(
                endpoints <- ifelse(target_number_shown == target_number_all, 0, 2),
                carrier_weight_inside <- target_number_shown * displacement_factor,
                carrier_norm <-
                   ( target_number_shown * displacement_factor
                    + endpoints * endpoint_carrier_weight) / (endpoints + target_number_shown),
                carrier_factor <- carrier_sensitivity *
                carrier_weight_inside / carrier_norm)),
        repulsive_endpoints=list(
            carrier_parameter = '
                real <lower=0, upper=max_sensitivity> endpoint_carrier_weight;
                // a factor multiplied by the summation going on elsewhere...
                real <lower=-10, upper=10> endpoint_carrier_factor;',
            carrier_var = '
                real carrier_weight_inside;
                real carrier_norm;
                real n_endpoints;',
            carrier_computation= '
                if (target_number_shown[n] == target_number_all[n])
                    n_endpoints <- 0;
                else n_endpoints <- 2;
                carrier_weight_inside <- target_number_shown[n] * displacement_factor;
                carrier_norm <-
                   ( carrier_weight_inside
                    + n_endpoints * endpoint_carrier_weight ) / (n_endpoints + target_number_shown[n]);
                carrier_factor <- (
                      carrier_sensitivity * target_number_shown[n] * displacement_factor
                    + carrier_weight_inside * n_endpoints *
                      endpoint_carrier_factor * carrier_weight_inside
                    ) / carrier_norm;
                ',
            carrier_R_computation=alist(
                endpoints <- ifelse(target_number_shown == target_number_all, 0, 2),
                carrier_weight_inside <- target_number_shown * displacement_factor,
                carrier_norm <-
                   ( target_number_shown * displacement_factor
                    + endpoints * endpoint_carrier_weight) / (endpoints + target_number_shown),
                carrier_factor <- (  (carrier_sensitivity * carrier_weight_inside)
                                   + (endpoint_carrier_factor * carrier_weight_inside)
                                  ) / carrier_norm
                )),
        local=list(
            carrier_parameter = '',
            carrier_var = '',
            carrier_computation =
                'carrier_factor <- 2*pi()*frac_spacing[n] * carrier_sensitivity;',
            carrier_R_computation=alist(
                carrier_factor <- 2*pi* frac_spacing * carrier_sensitivity)),
        windowed=list(
            carrier_parameter = 'real<lower=min(frac_spacing), upper=2*pi()> carrier_field;',
            carrier_var = '
                real frac_shown;
                real frac_in_carrier_field;
                ',
            carrier_computation = '
                frac_shown <- (target_number_shown[n]+0.0) / target_number_all[n];
                frac_in_carrier_field <- -blur * log(
                    exp(-2*pi()*frac_shown/blur) + exp(-carrier_field/blur));
                carrier_factor <-
                    target_number_all[n] * frac_in_carrier_field
                        * carrier_sensitivity / carrier_field;
                ',
            carrier_R_computation=alist(
                frac_shown <- target_number_shown / target_number_all,
                frac_in_carrier_field <- -blur * log(
                    exp(-2*pi*frac_shown/blur) + exp(-carrier_field/blur)),
                carrier_factor <-
                    target_number_all * frac_in_carrier_field
                * carrier_sensitivity / carrier_field))
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
    with(coefs, within(., {
      frac_spacing <- 2*pi/target_number_all
      ...(displacement_R_computation)
      ...(carrier_R_computation)
    }))
    , with(coefs, with(., within(list(), {
      link_displacement <- (displacement_factor * displacement)
      link_repulsion <- (repulsion * content
                         + nonlinearity * (content * abs(content)))
      link_carrier <- (content * carrier_factor)
      link <- bias + link_displacement + link_repulsion + link_carrier
      response <- plogis(link) * (1-lapse) + lapse/2
    })))
    , as.data.frame))

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

  summarizeData <- function(.data, ...) {
    env <- list2env(.data, parent = parent.frame())
    cols <- eval(substitute(alist(...)))
    for (col in names(cols)) {
      env[[col]] <- eval(cols[[col]], env)
    }
    dput(mget(ls(env), env))
  }

  mutateList <- function (.data, ...) {
    env <- list2env(.data, parent = parent.frame())
    cols <- eval(substitute(alist(...)))
    for (col in names(cols)) {
      env[[col]] <- eval(cols[[col]], env)
    }
    mget(ls(env), env)
  }

})

#these models screw up in some way and I omit them from taking up computer time.
losers <- c(
    "d_soft_windowed_c_local",
    "d_soft_local_c_local",
    "d_soft_local_endpoints_c_local",
    "d_soft_global_c_endpoints",
    "d_soft_global_c_local",
    "d_soft_global_c_windowed",
    "d_soft_windowed_c_windowed",
    "d_soft_windowed_c_local",
    "d_soft_global_c_global",
    "d_soft_windowed_c_global",
    "d_soft_global_c_global",
    "d_soft_windowed_c_global",
    "d_soft_global_c_endpoints",
    "d_soft_windowed_c_endpoints")

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
          Filter(f=function(x) !x$model_name %in% losers),
          lapply(compileModelEnv, outconn))
  })
}

run_as_command()

