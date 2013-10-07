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
lapse_limit <- 0.05
blur <- 0.2

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
      blur <- blur
      motion_energy_scale <- motion_energy_scale
      lapse_limit <- lapse_limit
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
  real blur;
  real motion_energy_scale;
  real lapse_limit;
 }

transformed data {
  real frac_spacing[N];

  for (n in 1:N)
      frac_spacing[n] <- 2*pi()./target_number_all[n];
}

parameters {
  real spacing_sensitivity;
  real<lower=0, upper=lapse_limit> lapse;
  real intercept;

  real<lower=0,upper=2*pi()> cs;
  real summation;
  real repulsion;
  real nonlinearity;
}

model {
 real beta_dx;
  real local_energy;
  real global_energy;
  real link_displacement;
  real link_repulsion;
  real link_summation;
  real link;

  //lapse ~ beta(1.5, 40); //what the shit, this explodes the bias

  for (n in 1:N) {
    beta_dx <- -blur * spacing_sensitivity
                * log(  exp(- frac_spacing[n] / blur)
                      + exp(- cs / blur));
    local_energy <- normalized_energy[n] / target_number_shown[n];
    global_energy <- normalized_energy[n];
    link_displacement <- beta_dx * displacement[n];
    link_repulsion <- ( repulsion * local_energy
                             + nonlinearity * local_energy * fabs(local_energy));
    link_summation <- global_energy * summation;
    link <- intercept
            + link_displacement
            + link_repulsion
            + link_summation;
    n_cw[n] ~ binomial( n_obs[n],
      inv_logit( link ) .* (1-lapse) + lapse/2);
  }
}'

stan_predict <- mkchain[., coefs](
    mutate(frac_spacing = 2*pi/target_number_all,
           normalized_energy = (norm_diff / motion_energy_scale),
           local_energy = normalized_energy / target_number_shown,
           global_energy = normalized_energy
           )
  , with(coefs, summarize(
      .
   , beta_dx = (-blur * spacing_sensitivity
                * log(  exp(- frac_spacing / blur)
                      + exp(- cs / blur)))
   , link_displacement = (beta_dx * displacement)
    , link_repulsion = ( repulsion * local_energy
                       + nonlinearity * (local_energy * abs(local_energy))
                       )
    , link_summation = global_energy * summation
    , link = ( intercept
             + link_displacement
             + link_repulsion
             + link_summation
             )
    , response = plogis(link) * (1-lapse) + lapse/2
    )))
