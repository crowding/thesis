splits <- c("subject", "content",
            "displacement",
            "target_number_all",
            "target_number_shown",
            "spacing", "eccentricity", "bias")

relevant <- splits %v% c("n_cw", "n_obs", "normalized_energy")

# the constant 0..... determined by regression on norm_diff
# me <- chain("motion_energy.csv", read.csv, add_energies)
# lm( norm_diff ~ I(content*target_number_all) - 1, data=data)$coef
motion_energy_scale <- 0.05136

stan_format <- mkchain(
    add_energies,
    mutate(normalized_energy =
           (norm_diff / motion_energy_scale)),
    subset(select=relevant),
    as.list,
    factorify,
    put(names(.), gsub('\\.', '_', names(.))),
    within({
      N <- length(subject_ix)
      motion_energy_scale <- motion_energy_scale
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
  real normalized_energy[N];
  int n_cw[N];
  int n_obs[N];
  real motion_energy_scale;
}

transformed data {
  real frac_spacing[N];

  for (n in 1:N)
      frac_spacing[n] <- 2*pi()./target_number_all[n];

}

parameters {
  real beta_dx;
  real<lower=0, upper=0.05> lapse;
  real bias;

  real<lower=0,upper=2*pi()> cs;
  real energy_weight; // how much "energy" to add to content
  real summation;
  real repulsion;
  real nonlinearity;
}

model {
  real crowdedness;
  real local_energy;
  real global_energy;
  real link_displacement;
  real link_repulsion;
  real link_summation;
  real link;

  //lapse ~ beta(1.5, 40); //what the shit, this explodes the bias

  for (n in 1:N) {
    local_energy <- normalized_energy[n] / target_number_shown[n] * energy_weight
                    + content[n] * (1-energy_weight);
    global_energy <- normalized_energy[n] * energy_weight
                    + content[n] * target_number_shown[n] * (1-energy_weight);
    crowdedness <- 2 - 2/(1+exp(-cs/frac_spacing[n]));
    link_displacement <- beta_dx * displacement[n] * crowdedness;
    link_repulsion <- (repulsion * local_energy
                       + nonlinearity * local_energy * abs(local_energy));
    link_summation <- global_energy * summation;
    link <- bias + link_displacement
            + link_repulsion
            + link_summation;
    n_cw[n] ~ binomial( n_obs[n],
      inv_logit( link ) .* (1-lapse) + lapse/2);
    }
}'

stan_predict <- mkchain[., coefs](
    mutate(frac_spacing = 2*pi/target_number_all,
           normalized_energy = norm_diff / motion_energy_scale,
           local_energy = normalized_energy / target_number_shown * energy_weight
                          + content * (1-energy_weight),
           global_energy = normalized_energy * energy_weight
                           + content * target_number_shown * (1-energy_weight)
           )
  , with(coefs, summarize(
      .
    , link_displacement = (beta_dx * displacement
                           * (2 - 2/(1+exp(-cs/frac_spacing))))
    , link_repulsion = (energy_repulsion * local_energy
                               + energy_nonlinearity * (
                                 local_energy * abs(local_energy)))
    , link_summation = norm_diff * energy_summation
    , link = (bias + link_displacement + link_repulsion + link_summation)
    , response = plogis(link) * (1-lapse) + lapse/2
    )))
