suppressPackageStartupMessages({
  library(vadr)
  library(plyr)
  library(rstan)
  library(knitr)
  library(R.cache)
  source("library.R")
  set_cppo('fast')
})

#displacement
scenarios <- list(
  #endpoint adjustments
  d=list(
    ## local=list( #for whatever reason spacing_sensitivity explodes
    ##     displacement_parameter='
    ##       real<lower=0.01, upper=100> max_sensitivity;',
    ##     displacement_parameter_name="max_sensitivity",
    ##     displacement_var = '',
    ##     displacement_computation = '
    ##       displacement_factor <-
    ##          (2 - 2/(1+exp(-1/spacing_sensitivity/frac_spacing[n])))
    ##          * max_sensitivity;',
    ##     displacement_R_computation = alist(
    ##       displacement_factor <-
    ##         (2 - 2/(1+exp(-1/spacing_sensitivity/frac_spacing)))
    ##* max_sensitivity)
    ##     ) ,
    soft_local=list(
      displacement_parameter='
          real <lower=spacing_sensitivity*min(frac_spacing),
                upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;',
      displacement_parameter_name="max_sensitivity",
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
    ## soft_local_boosted=list(
    ##     displacement_parameter='
    ##       real <lower=spacing_sensitivity*min(frac_spacing),
    ##             upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;
    ##       real <lower=0,upper=max_sensitivity> min_sensitivity;
    ##       ',
    ##     displacement_parameter_name=c("max_sensitivity", "min_sensitivity"),
    ##     displacement_var = '',
    ##     displacement_computation = '
    ##       displacement_factor <-
    ##          -blur * spacing_sensitivity
    ##           * log(  exp(- frac_spacing[n] / blur)
    ##                 + exp(- max_sensitivity / spacing_sensitivity
    ##                         * 2 * pi() / blur));
    ##       displacement_factor <- ( (displacement_factor/max_sensitivity)
    ##                               *(max_sensitivity-min_sensitivity)
    ##                               + min_sensitivity);
    ##     ',
    ##     displacement_R_computation = alist(
    ##         displacement_factor <- (
    ##             -blur * spacing_sensitivity
    ##             * log(   exp(- frac_spacing / blur)
    ##                   + exp(- max_sensitivity / spacing_sensitivity
    ##                         * 2 * pi / blur))))),
##     soft_local_endpoints=list(
##       displacement_parameter='
##           real <lower=spacing_sensitivity*min(frac_spacing),
##                 upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;
##           real<lower=0, upper=max_sensitivity> endpoint_sensitivity;
## ',
##       displacement_parameter_name=c("max_sensitivity", "endpoint_sensitivity"),
##       displacement_var = '
##           real endpoints;
##           real raw_spacing_sensitivity;
##           real envelope_interior;
##           real envelope_norm;
##           ',
##       displacement_computation = '
##         raw_spacing_sensitivity <- -blur * spacing_sensitivity
##             * log(  exp(- frac_spacing[n] / blur)
##                   + exp(- max_sensitivity / spacing_sensitivity
##                           * 2 * pi() / blur));
##         if (target_number_shown[n] == target_number_all[n])
##           endpoints <- 0;
##         else endpoints <- 2;
##         envelope_interior <- pow(raw_spacing_sensitivity, 2)
##                              * target_number_shown[n];
##         envelope_norm <- (  endpoints * endpoint_sensitivity
##                             + target_number_shown[n] * raw_spacing_sensitivity);
##         displacement_factor <- envelope_interior / envelope_norm;
##         ',
##       displacement_R_computation = alist(
##         endpoints <- ifelse(target_number_shown == target_number_all, 0, 2),
##         raw_spacing_sensitivity <- (
##           -blur * spacing_sensitivity
##           * log(  exp(- frac_spacing / blur)
##                 + exp(- max_sensitivity / spacing_sensitivity
##                       * 2 * pi / blur))),
##         envelope_interior <- raw_spacing_sensitivity^2 * target_number_shown,
##         envelope_norm <- (endpoints * endpoint_sensitivity
##                           + target_number_shown * raw_spacing_sensitivity),
##         displacement_factor <- envelope_interior / envelope_norm
##         )),
    soft_global=list(
      displacement_parameter='
          real <lower=spacing_sensitivity*min(frac_spacing),
                upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;',
      displacement_parameter_name=c("max_sensitivity"),
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
                      spacing_sensitivity * 2 * pi / blur))))),
    soft_hemi=list(
      displacement_parameter='
        real <lower=spacing_sensitivity*min(frac_spacing),
              upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;',
      displacement_parameter_name=c("max_sensitivity"),
      displacement_var = '
        real extent;
        real covered;
        real eff_number;',
      displacement_computation = '
        extent <- (0.0 + target_number_shown[n]) / target_number_all[n]; //[0,1]
        covered <- fmin(extent, 0.5);
        //one shown equiv to two in full circle:
        eff_number <- 4*pi() * covered / frac_spacing[n];
        displacement_factor <-
           -blur * spacing_sensitivity
            * log(  exp(- eff_number / blur)
                  + exp(- max_sensitivity / spacing_sensitivity
                          * 2 * pi() / blur));',
      displacement_R_computation = alist(
        extent <- (0.0 + target_number_shown) / target_number_all,
        covered <- pmin(extent, 0.5),
        eff_number <- 4*pi*covered / frac_spacing,
        displacement_factor <- (
          -blur *  spacing_sensitivity
          * log(  exp(- eff_number / blur)
                + exp(- max_sensitivity /
                      spacing_sensitivity * 2 * pi / blur)))))
    ,
    soft_windowed=list(
      displacement_parameter='
          real <lower=spacing_sensitivity*min(frac_spacing),
                upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;
          real <lower=min(frac_spacing), upper=2*pi()> displacement_field;',
      displacement_parameter_name=c("max_sensitivity", "displacement_field"),
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
                  ##CARRIER
    c=list(
      global=list(
        carrier_parameter = '',
        carrier_parameter_name=c(),
        carrier_var = '',
        carrier_computation=
        'carrier_factor <- target_number_shown[n] * carrier_sensitivity;',
        carrier_R_computation=alist(
          carrier_factor <- target_number_shown * carrier_sensitivity)),
      hemi=list(
        carrier_parameter = '',
        carrier_parameter_name=c(),
        carrier_var = '
          real c_extent;
          real c_covered;
          real c_n_eff;
        ',
        carrier_computation= '
          c_extent <- (0.0 + target_number_shown[n]) / target_number_all[n]; //[0,1]
          c_covered <- fmin(c_extent, 0.5);
          c_n_eff <- pi() * c_covered / frac_spacing[n];
          //one shown equiv to two in full circle:
          carrier_factor <- c_n_eff * carrier_sensitivity;',
        carrier_R_computation=alist(
          c_extent <- (0.0 + target_number_shown) / target_number_all,
          c_covered <- pmin(c_extent, 0.5),
          c_n_eff <- pi * c_covered / frac_spacing,
          carrier_factor <- c_n_eff * carrier_sensitivity)),
      ## endpoints=list(
      ##   carrier_parameter = '
      ##           real <lower=0, upper=max_sensitivity> endpoint_carrier_weight;',
      ##   carrier_parameter_name=c("endpoint_carrier_weight"),
      ##   carrier_var = '
      ##           real carrier_weight_inside;
      ##           real carrier_norm;
      ##           real n_endpoints;',
      ##       carrier_computation= '
      ##         if (target_number_shown[n] == target_number_all[n])
      ##             n_endpoints <- 0;
      ##         else n_endpoints <- 2;
      ##         carrier_weight_inside <- target_number_shown[n]
      ##                                  * displacement_factor;
      ##         carrier_norm <-
      ##            ( target_number_shown[n] * displacement_factor
      ##             + n_endpoints * endpoint_carrier_weight ) /
      ##            (n_endpoints + target_number_shown[n]);
      ##         carrier_factor <- carrier_sensitivity *
      ##                           carrier_weight_inside / carrier_norm;
      ##         ',
      ##       carrier_R_computation=alist(
      ##         endpoints <- ifelse(target_number_shown == target_number_all,
      ##                             0, 2),
      ##         carrier_weight_inside <- (target_number_shown
      ##                                   * displacement_factor),
      ##         carrier_norm <- (
      ##           target_number_shown * displacement_factor
      ##           + endpoints * endpoint_carrier_weight
      ##           ) / (endpoints + target_number_shown),
      ##         carrier_factor <- (carrier_sensitivity
      ##                            * carrier_weight_inside) / carrier_norm)),
      ##   repulsive_endpoints=list(
      ##     carrier_parameter = '
      ##       real <lower=0, upper=max_sensitivity> endpoint_carrier_weight;
      ##       // a factor multiplied by the summation going on elsewhere...
      ##       real <lower=-10, upper=10> endpoint_carrier_factor;',
      ##     carrier_parameter_name=c("endpoint_carrier_weight",
      ##                              "endpoint_carrier_factor"),
      ##     carrier_var = '
      ##       real carrier_weight_inside;
      ##       real carrier_norm;
      ##       real n_endpoints;',
      ##     carrier_computation= '
      ##       if (target_number_shown[n] == target_number_all[n])
      ##           n_endpoints <- 0;
      ##       else n_endpoints <- 2;
      ##       carrier_weight_inside <- target_number_shown[n]
      ##                                * displacement_factor;
      ##       carrier_norm <- (
      ##           carrier_weight_inside + n_endpoints * endpoint_carrier_weight
      ##         ) / (
      ##           n_endpoints + target_number_shown[n]
      ##         );
      ##       carrier_factor <- (
      ##         (
      ##           carrier_sensitivity * carrier_weight_inside
      ##         ) + (
      ##           n_endpoints * endpoint_carrier_factor *
      ##           carrier_sensitivity * endpoint_carrier_weight
      ##         )
      ##       ) / carrier_norm;
      ##       ',
      ##     carrier_R_computation=alist(
      ##       endpoints <- ifelse(target_number_shown == target_number_all, 0, 2),
      ##       carrier_weight_inside <- target_number_shown * displacement_factor,
      ##       carrier_norm <- (
      ##         carrier_weight_inside + endpoints * endpoint_carrier_weight) / (
      ##           endpoints + target_number_shown),
      ##       carrier_factor <- (
      ##         (
      ##           carrier_sensitivity * carrier_weight_inside
      ##         ) + (
      ##           endpoints * endpoint_carrier_factor
      ##           * carrier_sensitivity * endpoint_carrier_weight
      ##         )
      ##       ) / carrier_norm
      ##       )),
        local=list(
          carrier_parameter = '',
          carrier_parameter_name=c(),
          carrier_var = '',
          carrier_computation = '
            carrier_factor <- 2*pi()*frac_spacing[n] * carrier_sensitivity;',
          carrier_R_computation=alist(
            carrier_factor <- 2*pi* frac_spacing * carrier_sensitivity)),
      windowed=list(
        carrier_parameter = '
          real<lower=min(frac_spacing), upper=2*pi()> carrier_field;',
        carrier_parameter_name=c("carrier_field"),
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
          * carrier_sensitivity / carrier_field))),
  #ENDPOINT
  e = list(
    none=list(
      endpoint_parameter='',
      endpoint_parameter_name=c(),
      endpoint_var = '',
      endpoint_computation='',
      endpoint_R_computation=''),
    B=list(
      #bias
      endpoint_parameter = '
          real<lower=-5, upper=5> endpoint_intercept;',
      endpoint_parameter_name=c("endpoint_intercept"),
      endpoint_var = '
          real n_endpoints;',
      endpoint_computation= '
          if (target_number_shown[n] == target_number_all[n])
            n_endpoints <- 0;
          else n_endpoints <- 2;
          offset <- offset + endpoint_intercept * n_endpoints;',
      endpoint_R_computation=alist(
        n_endpoints <- ifelse(target_number_shown == target_number_all, 0, 2),
        offset <- offset + endpoint_intercept * n_endpoints)),
    R=list(
      #repulsive
      endpoint_parameter = '
          real<lower=-5, upper=5> endpoint_repulsion;',
      endpoint_parameter_name=c("endpoint_repulsion"),
      endpoint_var = '
          real n_endpoints;',
      endpoint_computation= '
          if (target_number_shown[n] == target_number_all[n])
            n_endpoints <- 0;
          else n_endpoints <- 2;
          carrier_factor <- carrier_factor + n_endpoints * endpoint_repulsion;',
      endpoint_R_computation=alist(
        n_endpoints <- ifelse(target_number_shown == target_number_all, 0, 2),
        carrier_factor <- carrier_factor + (n_endpoints * endpoint_repulsion))),
    RE=list(
      #repulsive with extent factor
      endpoint_parameter = '
          real<lower=-5, upper=5> endpoint_repulsion;
          real<lower=-5, upper=5> extent_summation;',
      endpoint_parameter_name=c("endpoint_repulsion", "extent_summation"),
      endpoint_var = '
          real n_endpoints;
          real e_extent;',
      endpoint_computation= '
          e_extent <- (0.0 + target_number_shown[n]) / target_number_all[n];
          if (target_number_shown[n] == target_number_all[n])
            n_endpoints <- 0;
          else n_endpoints <- 2;
          carrier_factor <- carrier_factor + n_endpoints * endpoint_repulsion
            + e_extent * extent_summation * target_number_shown[n] * content[n];',
      endpoint_R_computation=alist(
        e_extent <- target_number_shown / target_number_all,
        n_endpoints <- ifelse(target_number_shown == target_number_all, 0, 2),
        carrier_factor <- carrier_factor + (
          n_endpoints * endpoint_repulsion
          + e_extent * extent_summation * target_number_shown * content))),
    A=list(
      #additive
      endpoint_parameter = '
          real<lower=-5, upper=5> endpoint_summation;',
      endpoint_parameter_name=c("endpoint_summation"),
      endpoint_var = '
          real n_endpoints;',
      endpoint_computation= '
          if (target_number_shown[n] == target_number_all[n])
            n_endpoints <- 0;
          else n_endpoints <- 2;
          carrier_factor <- carrier_factor
            + n_endpoints * endpoint_summation * target_number_shown[n];',
      endpoint_R_computation=alist(
        n_endpoints <- ifelse(target_number_shown == target_number_all, 0, 2),
        carrier_factor <- carrier_factor + (
          n_endpoints * endpoint_summation * target_number_shown))),
    AB=list(
      #additive biased
      endpoint_parameter = '
          real<lower=-5, upper=5> endpoint_summation;
          real<lower=-5, upper=5> endpoint_intercept;',
      endpoint_parameter_name=c("endpoint_summation", "endpoint_intercept"),
      endpoint_var = '
          real n_endpoints;',
      endpoint_computation= '
          if (target_number_shown[n] == target_number_all[n])
            n_endpoints <- 0;
          else n_endpoints <- 2;
          carrier_factor <- carrier_factor
            + n_endpoints * endpoint_summation * target_number_shown[n];
          offset <- offset + endpoint_intercept * n_endpoints;',
      endpoint_R_computation=alist(
        n_endpoints <- ifelse(target_number_shown == target_number_all, 0, 2),
        carrier_factor <- carrier_factor + (
          n_endpoints * endpoint_summation * target_number_shown),
        offset <- offset + endpoint_intercept * n_endpoints)),
    AE=list(
      #additive with extent factor
      endpoint_parameter = '
          real<lower=-5, upper=5> endpoint_summation;
          real<lower=-10, upper=10> extent_summation;',
      endpoint_parameter_name=c("endpoint_summation", "extent_summation"),
      endpoint_var = '
          real n_endpoints;
          real e_extent;',
      endpoint_computation= '
          e_extent <- (0.0 + target_number_shown[n]) / target_number_all[n]; //[0,1]
          if (target_number_shown[n] == target_number_all[n])
            n_endpoints <- 0;
          else n_endpoints <- 2;
          carrier_factor <- carrier_factor
            + n_endpoints * endpoint_summation * target_number_shown[n]
            + e_extent * extent_summation * target_number_shown[n] * content[n];',
      endpoint_R_computation=alist(
        e_extent <- target_number_shown / target_number_all,
        n_endpoints <- ifelse(target_number_shown == target_number_all, 0, 2),
        carrier_factor <- carrier_factor + (
          n_endpoints * endpoint_summation * target_number_shown
          + e_extent * extent_summation * target_number_shown * content))),
    RA=list(
      #repulsive and
      endpoint_parameter = '
          real<lower=-5, upper=5> endpoint_repulsion;
          real<lower=-5, upper=5> endpoint_summation;',
      endpoint_parameter_name=c("endpoint_summation", "endpoint_repulsion"),
      endpoint_var = '
          real n_endpoints;',
      endpoint_computation= '
          if (target_number_shown[n] == target_number_all[n])
            n_endpoints <- 0;
          else n_endpoints <- 2;
          carrier_factor <-
            target_number_shown[n] * carrier_sensitivity
            + n_endpoints * endpoint_summation * target_number_shown[n]
            + n_endpoints * endpoint_repulsion;',
      endpoint_R_computation=alist(
        n_endpoints <- ifelse(target_number_shown == target_number_all, 0, 2),
        carrier_factor <- (
          target_number_shown * carrier_sensitivity
          + n_endpoints * endpoint_repulsion
          + n_endpoints * endpoint_summation * target_number_shown))),
    RAE=list(
      #repulsive and additive with extent factor
      endpoint_parameter = '
          real<lower=-5, upper=5> endpoint_repulsion;
          real<lower=-5, upper=5> endpoint_summation;
          real<lower=-5, upper=5> extent_summation;',
      endpoint_parameter_name=c(
        "extent_repulsion", "endpoint_summation", "extent_summation"),
      endpoint_var = '
          real n_endpoints;
          real e_extent;',
      endpoint_computation= '
          e_extent <- (0.0 + target_number_shown[n]) / target_number_all[n]; //[0,1]
          if (target_number_shown[n] == target_number_all[n])
            n_endpoints <- 0;
          else n_endpoints <- 2;
          carrier_factor <- carrier_factor
            + n_endpoints * endpoint_summation * target_number_shown[n]
            + e_extent * extent_summation * target_number_shown[n] * content[n]
            + n_endpoints * endpoint_repulsion;',
      endpoint_R_computation=alist(
        e_extent <- target_number_shown / target_number_all,
        n_endpoints <- ifelse(target_number_shown == target_number_all, 0, 2),
        carrier_factor <- carrier_factor + (
          n_endpoints * endpoint_summation * target_number_shown
          + n_endpoints * endpoint_repulsion
          + e_extent * extent_summation * target_number_shown * content)))))

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
  int trial_id[N];
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
  real intercept;
  real carrier_sensitivity;
  real <lower=0>spacing_sensitivity;
  {{displacement_parameter}}
  {{carrier_parameter}}
  {{endpoint_parameter}}
  real repulsion;
  real nonlinearity;
}
model {
  {{displacement_var}}
  {{carrier_var}}
  {{endpoint_var}}
  real carrier_factor;
  real displacement_factor;
  real offset;
  real link_displacement;
  real link_repulsion;
  real link_summation;
  real link_offset;
  real link;
  for (n in 1:N) {
    offset <- 0;
    {{displacement_computation}}
    {{carrier_computation}}
    {{endpoint_computation}}
    link_displacement <- displacement[n] * displacement_factor;
    link_repulsion <- (repulsion * content[n]
                       + nonlinearity * (content[n] * fabs(content[n])));
    link_summation <- content[n] * carrier_factor;
    link_offset <- intercept + offset;
    link <- link_offset + link_displacement + link_repulsion + link_summation;
    n_cw[n] ~ binomial( n_obs[n],
      inv_logit( link ) .* (1-lapse) + lapse/2);
  }
}
generated quantities {
  real displacement_sens[N];
  real carrier_sens[N];
  real link_displacement[N];
  real link_repulsion[N];
  real link_summation[N];
  real link_offset[N];
  real trial_id_stan[N];
  real link[N];
  real fit[N];
  for (n in 1:N) {
    {{displacement_var}}
    {{carrier_var}}
    {{endpoint_var}}
    real carrier_factor;
    real displacement_factor;
    real offset;
    offset <- 0;
    {{displacement_computation}}
    {{carrier_computation}}
    {{endpoint_computation}}
    displacement_sens[n] <- displacement_factor;
    carrier_sens[n] <- carrier_factor;
    link_displacement[n] <- displacement[n] * displacement_factor;
    link_repulsion[n] <- (repulsion * content[n]
                          + nonlinearity * (content[n] * fabs(content[n])));
    link_summation[n] <- content[n] * carrier_factor;
    link_offset[n] <- intercept + offset;
    link[n] <- link_offset[n] + link_displacement[n]
               + link_repulsion[n] + link_summation[n];
    trial_id_stan[n] <- trial_id[n];
    fit[n] <- inv_logit( link[n] ) * (1-lapse) + lapse/2;
  }
}'


predictorTemplate <- quote(stan_predict <- mkchain[., coefs](
    with(coefs, within(., {
      frac_spacing <- 2*pi/target_number_all
      offset <- 0
      ...(displacement_R_computation)
      ...(carrier_R_computation)
      ...(endpoint_R_computation)
    }))
    , with(coefs, with(., within_df(list(), {
      link_displacement <- (displacement_factor * displacement)
      link_repulsion <- (repulsion * content
                         + nonlinearity * (content * abs(content)))
      link_carrier <- (content * carrier_factor)
      link_offset <- (intercept+offset)
      link <- link_offset + link_displacement + link_repulsion + link_carrier
      response <- plogis(link) * (1-lapse) + lapse/2
    })))
    , as.data.frame))

otherFunctions <- quote({
  filter_data <- mkchain(
      subset(exp_type %in% c("spacing", "content", "numdensity"))
      , match_df(., subset(count(., "subject"), freq>2000), on="subject")
    )

  within_df <- function(df, expr, enclos=parent.frame()) {
    if (is.null(names(df))) {
      names(df) <- rep("", length(df))
    }
    e <- list2env(df, parent=enclos)
    eval(substitute(expr), e)
    as.data.frame(as.list(e))
  }

  splits <- c("subject", "content",
              "displacement",
              "target_number_all",
              "target_number_shown",
              "spacing", "eccentricity", "bias")

  pars <- qe(c("lapse", "intercept",
               "spacing_sensitivity", "carrier_sensitivity",
               "repulsion", "nonlinearity",
               ..(specifications$displacement_parameter_name),
               ..(specifications$carrier_parameter_name),
               ..(specifications$endpoint_parameter_name)))

  relevant <- splits %v% c("n_cw", "n_obs")

  blur <- 0.2
  lapse_limit <- 0.05

  format_data0 <- format_data
  format_data <- mkchain[., menergy](format_data0(menergy),
                                     mutate(format_data,
                                            trial_id=seq_along(displacement)))

  stan_format <- mkchain(
    .[c(relevant, "trial_id")],
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
#...
losers <- c(
  "d_(?:soft_windowed|soft_hemi|soft_global)_c_.*_e_(?:RE|AE|RAE|RR)",
  "d_.*_c_(?:windowed|local|global)_e_(?:RE|AE|RAE|RR)",
  "d_soft_global_c_windowed_e.*",

  #memory a splode
  "d_soft_local_c_hemi_e_B",

  "d_soft_local_c_hemi_e_RA",
  "d_soft_global_c_global_e_A",

  "d_soft_windowed_c_endpoints_e.*",
  "d_soft_windowed_c_global_e.*",
  "d_soft_windowed_c_global_e.*",
  "d_soft_windowed_c_local_e.*",
  "d_soft_windowed_c_windowed_e.*",
  "d_soft_windowed_c_local_e.*",

  "d_soft_local_endpoints_c_local.*",
  "d_soft_windowed_c_repulsive_endpoints.*",
  "d_soft_global_c_endpoints.*",
  "d_soft_global_c_endpoints.*",

  "d_soft_windowed_c_hemi_e_RA.*" #memory explodes during optimization???
#    "d_soft_global_c_global", #these suck but they need to suck for comparison.
#    "d_soft_global_c_local",
#    "d_soft_local_c_local",
  )

makeModelEnv <- function(selection=lapply(scenarios, mkchain(names, .[[1]])),
                         scenarios=parent.env(environment())$scenarios,
                         envir = new.env(parent=globalenv())) {
  #substitute bits and pieces into the
  envir <- as.environment(envir)
  with(envir, {
    base::source("stanFunctions.R", local=TRUE)
  })
  substituting.env <- new.env(parent=envir)
  envir$specifications <- substituting.env
  mapply(names(selection), selection,
         FUN=function(name, value) {
           assignments <- scenarios[[name]][[value]]
           mapply(names(assignments), assignments,
                  FUN = assign,
                  MoreArgs = list(envir=substituting.env)
                  )
         })
  envir$model_code <- do.call(knit_expand,
                              list(text=modelTemplate),
                              envir=substituting.env)
  envir$model_name <- chain(selection, rbind(names(.),.),
                            paste(collapse="_"))
  envir$stan_predict <- eval(do.call(vadr::qq,
                                     list(predictorTemplate),
                                     envir=substituting.env),
                             envir)
  eval(otherFunctions, envir)
  envir
}

using <- function(conn, fn) {
  fn(ammoc(conn, on.exit(close(conn))))
}

compileModelEnv <- function(envir, listing_file) {
  filename <- paste0("models/", envir$model_name, ".stan.RData")
  hash <- list("stan_model",
               model_name=envir$model_name,
               model_code=envir$model_code,
               stan_version=packageVersion("rstan"))
  ((model <- loadCache(hash)) %||%
   ammoc(model <- stan_model(model_name=envir$model_name,
                             model_code=envir$model_code),
         saveCache(model, hash)))
  envir$model <- model
  save(file=filename, list=ls(envir), envir=envir)
  message(envir$model_name)
  filename
}

main <- function(outfile='NineModels.list') {
  using(file(outfile, open="w"), function(outconn) {
    chain(scenarios,
          lapply(names),
          do.call(expand.grid, .),
          as.matrix,
          alply(1, makeModelEnv),
          filter_models,
          llply(compileModelEnv, outconn, .parallel=TRUE),
          unlist,
          writeLines(outconn))
  })
  unlink("overshoots.txt")
}

filter_models <- mkchain(
  models=.,
  vapply(function(x)x$model_name, ""),
  lapply(losers, grepl, .),
  Reduce(`|`, .),
  models[!.],
  {print(vapply(., function(x)x$model_name, "")); .}
  )

run_as_command()
