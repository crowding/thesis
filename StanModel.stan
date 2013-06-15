data {
  int<lower=0> N;
  int subject_n;
  int subject_ix[N];
  int target_number_all[N];
  int target_number_shown[N];
  real displacement[N];
  real content[N];
  real spacing[N];
  real eccentricity[N];
  int n_cw[N];
  int n_obs[N];
}

transformed data {
  real frac_spacing[N];

  for (n in 1:N)
      frac_spacing[n] <- 2*pi()./target_number_all[n];
}

parameters {
  real beta_dx;
  real<lower=0, upper=0.1> lapse;
  real bias;

  real<lower=0,upper=2*pi()> cs;
  real content_global_summation;
  real content_repulsion;
  real content_nonlinearity;
}

model {
  real crowdedness;

  for (n in 1:N) {
    crowdedness <- 2 - 2/(1+exp(-cs/frac_spacing[n]));

    n_cw[n] ~ binomial(
      n_obs[n],
      inv_logit(
        beta_dx * displacement[n] * crowdedness
        + target_number_all[n] * content[n] * content_global_summation
        + content_repulsion * content[n]
        + content_nonlinearity * (content[n] * abs(content[n]))
        + bias
        ) .* (1-lapse) + lapse/2);
  }
}
