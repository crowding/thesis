data {
  int<lower=0> N;
  int subject_n;
  int subject_ix[N];
  real displacement[N];
  real content[N];
  real spacing[N];
  real eccentricity[N];
  int n_cw[N];
  int n_obs[N];
}

parameters {
  real beta_dx;
  real<lower=0, upper=0.1> lapse;
  real bias;

  //  real<lower=0,upper={{2*pi*max(radius)}}> cs;
  //  real global_content;
  //  real content_repulsion;
  //  real content_nonlin;
}

model {
  for (n in 1:N) {
    n_cw[n] ~ binomial(
      n_obs[n],
      inv_logit(beta_dx * displacement[n] + bias) .* (1-lapse) + lapse/2);
  }
}
