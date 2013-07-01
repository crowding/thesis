splits <- c("subject", "content",
            "displacement",
            "target_number_all",
            "target_number_shown",
            "spacing", "eccentricity", "bias")

relevant <- splits %v% c("n_cw", "n_obs")

stan_format <- mkchain(
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
  real<lower=0, upper=0.05> lapse;
  real bias;

  real<lower=0,upper=2*pi()> cs;
  real summation;
  real repulsion;
  real nonlinearity;
}

model {
  real crowdedness;
  real link_displacement;
  real link_repulsion;
  real link_summation;
  real link;

  //lapse ~ beta(1.5, 40); //what the shit, this explodes the bias

  for (n in 1:N) {
    crowdedness <- 2 - 2/(1+exp(-cs/frac_spacing[n]));
    link_displacement <- beta_dx * displacement[n] * crowdedness;
    link_repulsion <- (repulsion * content[n]
                       + nonlinearity * (content[n] * abs(content[n])));
    link_summation <- target_number_shown[n] * content[n] * summation;
    link <- bias + link_displacement + link_repulsion + link_summation;
    n_cw[n] ~ binomial( n_obs[n],
      inv_logit( link ) .* (1-lapse) + lapse/2);
  }

}'

 stan_predict <- mkchain[., coefs](
     mutate(frac_spacing = 2*pi/target_number_all)
   , with(coefs, summarize(
       .
     , link_displacement = (beta_dx * displacement
                            * (2 - 2/(1+exp(-cs/frac_spacing))))
     , link_repulsion = (repulsion * content
                         + nonlinearity * (content * abs(content)))
     , link_summation = (target_number_all * content * summation)
     , link = bias + link_displacement + link_repulsion + link_summation
     , response = plogis(link) * (1-lapse) + lapse/2
     )))
