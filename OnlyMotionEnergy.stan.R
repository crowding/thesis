splits <- c("subject", "content",
            "displacement",
            "target_number_all",
            "target_number_shown",
            "spacing", "eccentricity", "bias")

relevant <- splits %v% c("n_cw", "n_obs", "energy_diff", "norm_diff")

stan_format <- mkchain(
    add_energies,
    subset(select=relevant),
    as.list,
    factorify,
    put(names(.), gsub('\\.', '_', names(.))),
    within({
      N <- length(subject_ix)
    }))

model_code <- '
data {
  int<lower=0> N;
  int target_number_all[N];
  int target_number_shown[N];
  real displacement[N];
  real content[N];
  real spacing[N];
  real eccentricity[N];
  real energy_diff[N];
  real norm_diff[N];
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
  real energy_summation;
  real energy_repulsion;
  real energy_nonlinearity;
}

model {
  real crowdedness;
  real local_energy;
  real link_displacement;
  real link_energy_repulsion;
  real link_energy_summation;
  real link;

  for (n in 1:N) {
    crowdedness <- 2 - 2/(1+exp(-cs/frac_spacing[n]));
    local_energy <- norm_diff[n] / target_number_shown[n];
    link_displacement <- beta_dx * displacement[n] * crowdedness;
    link_energy_repulsion <- (energy_repulsion * local_energy
                              + energy_nonlinearity * local_energy * abs(local_energy));
    link_energy_summation <- norm_diff[n] * energy_summation;
    link <- bias + link_displacement
            + link_energy_repulsion
            + link_energy_summation;
    n_cw[n] ~ binomial( n_obs[n],
      inv_logit( link ) .* (1-lapse) + lapse/2);
  }
}'

stan_predict <- mkchain[., coefs](
    mutate(frac_spacing = 2*pi/target_number_all,
           local_energy = norm_diff / target_number_shown)
  , with(coefs, summarize(
      .
    , link_displacement = (beta_dx * displacement
                           * (2 - 2/(1+exp(-cs/frac_spacing))))
    , link_energy_repulsion = (energy_repulsion * local_energy
                               + energy_nonlinearity * (
                                 local_energy * abs(local_energy)))
    , link_energy_summation = norm_diff * energy_summation
    , link = (bias + link_displacement
              + link_repulsion + link_energy_repulsion
              + link_summation + link_energy_summation)
    , response = plogis(link) * (1-lapse) + lapse/2
    )))
