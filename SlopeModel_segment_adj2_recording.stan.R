filter_data <- mkchain(
  subset(exp_type %in% c("spacing", "content", "numdensity"))
  , match_df(., subset(count(., "subject"), freq>2000), on="subject"))

splits <- c("subject", "content",
            "displacement",
            "target_number_all",
            "target_number_shown",
            "spacing", "eccentricity", "bias")

relevant <- splits %v% c("n_cw", "n_obs")

lapse_limit <- 0.05

stan_format <- mkchain(
    .[c(relevant, "trial_id")],
    as.list,
    factorify,
    put(names(.), gsub('\\.', '_', names(.))),
    within({
      N <- length(subject_ix)
      lapse_limit <- lapse_limit
    }))

format_data0 <- format_data
format_data <- mkchain[., menergy](format_data0(menergy),
                                   mutate(format_data,
                                          trial_id=seq_along(displacement)))

model_code <- '
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
  real lapse_limit;
}

transformed data {
  real frac_spacing[N];

  for (n in 1:N)
      frac_spacing[n] <- 2*pi()./target_number_all[n];
}

parameters {
  real beta_dx;
  real<lower=0, upper=lapse_limit> lapse;
  real intercept;

  real<lower=0,upper=2*pi()> cs;
  real summation;
  real repulsion;
  real nonlinearity;
  real<lower=-100, upper=100> endpoint_summation;
  real<lower=-100, upper=100> endpoint_repulsion;
}

model {
  real crowdedness;
  real link_displacement;
  real link_repulsion;
  real link_summation;
  real link_endpoint;
  real link[N];
  real endpoints;
  real fit[N];

  for (n in 1:N) {
    if (target_number_shown[n] == target_number_all[n])
      endpoints <- 0;
    else endpoints <- 2;

    crowdedness <- 2 - 2/(1+exp(-cs/frac_spacing[n]));
    link_displacement <- beta_dx * displacement[n] * crowdedness;
    link_repulsion <- (repulsion * content[n]
                       + nonlinearity * (content[n] * fabs(content[n])));
    link_endpoint <- endpoints * (
          endpoint_repulsion * content[n]
          + endpoint_summation * content[n] * target_number_shown[n]);
    link_summation <- target_number_shown[n] * content[n] * summation;
    link[n] <- intercept + link_displacement + link_repulsion
               + link_summation + link_endpoint;
    fit[n] <- inv_logit( link[n] ) .* (1-lapse) + lapse/2;
    n_cw[n] ~ binomial( n_obs[n], fit[n]);
  }
}

generated quantities {
  real fit[N];
  real link[N];
  real id[N];
  real link_repulsion[N];
  real link_nonlinearity[N];
  real trial_id_stan[N];
  for (n in 1:N) {
    real crowdedness;
    real link_displacement;
    real link_summation;
    real link_endpoint;
    real endpoints;

    if (target_number_shown[n] == target_number_all[n])
      endpoints <- 0;
    else endpoints <- 2;

    crowdedness <- 2 - 2/(1+exp(-cs/frac_spacing[n]));
    link_displacement <- beta_dx * displacement[n] * crowdedness;
    link_nonlinearity[n] <- nonlinearity * (content[n] * fabs(content[n]));
    link_repulsion[n] <- (repulsion * content[n]);
    link_endpoint <- endpoints * (
          endpoint_repulsion * content[n]
          + endpoint_summation * content[n] * target_number_shown[n]);
    link_summation <- target_number_shown[n] * content[n] * summation;
    link[n] <- intercept + link_displacement + link_repulsion[n]
             + link_summation + link_endpoint + link_nonlinearity[n];
    fit[n] <- inv_logit( link[n] ) .* (1-lapse) + lapse/2;
    trial_id_stan[n] <- trial_id[n];
  }
}'

pars <- c("beta_dx", "lapse", "intercept", "cs",
          "summation", "repulsion", "nonlinearity",
          "endpoint_summation", "endpoint_repulsion")

stan_predict <- mkchain[., coefs](
  mutate(frac_spacing = 2*pi/target_number_all,
         endpoints = ifelse(target_number_shown == target_number_all, 0, 2))
  , with(coefs, summarize(
    .
    , link_displacement = (beta_dx * displacement
                           * (2 - 2/(1+exp(-cs/frac_spacing))))
    , link_nonlinearity = nonlinearity * content * abs(content)
    , link_repulsion = repulsion * content
    , link_summation = (target_number_shown * content * summation)
    , link_endpoint = endpoints * (
        endpoint_repulsion * content +
        endpoint_summation * content * target_number_shown)
    , link = (intercept + link_displacement + link_repulsion
              + link_nonlinearity
              + link_summation + link_endpoint)
    , response = plogis(link) * (1-lapse) + lapse/2
    )))


