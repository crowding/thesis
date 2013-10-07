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

 stan_predict <- mkchain[., coefs](
     mutate(frac_spacing = 2*pi/target_number_all)
   , with(coefs, summarize(
       .
     , beta_dx = (-blur * spacing_sensitivity
                  * log(  exp(- frac_spacing / blur)
                        + exp(- max_sensitivity / spacing_sensitivity * 2 * pi / blur)))
     , link_displacement = (beta_dx * displacement)
     , link_repulsion = (repulsion * content
                         + nonlinearity * (content * abs(content)))
     , link_summation = (target_number_all * content * summation)
     , link = intercept + link_displacement + link_repulsion + link_summation
     , response = plogis(link) * (1-lapse) + lapse/2
     )))

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
  real blur;
  real lapse_limit;
}

transformed data {
  real frac_spacing[N];

  for (n in 1:N)
      frac_spacing[n] <- 2*pi()./target_number_all[n];
}

parameters {
  real <lower=0>spacing_sensitivity;
  real<lower=0, upper=lapse_limit> lapse;
  real intercept;

  real<lower=0,upper=spacing_sensitivity*max(frac_spacing)> max_sensitivity;
  real summation;
  real repulsion;
  real nonlinearity;
}

model {
  real beta_dx;
  real link_displacement;
  real link_repulsion;
  real link_summation;
  real link;

  //lapse ~ beta(1.5, 40); //what the shit, this explodes the bias

  for (n in 1:N) {
    beta_dx <- -blur * spacing_sensitivity
                * log(  exp(- frac_spacing[n] / blur)
                      + exp(- max_sensitivity / spacing_sensitivity * 2 * pi() / blur));
    link_displacement <- beta_dx * displacement[n];
    link_repulsion <- (repulsion * content[n]
                       + nonlinearity * (content[n] * fabs(content[n])));
    link_summation <- target_number_shown[n] * content[n] * summation;
    link <- intercept + link_displacement + link_repulsion + link_summation;
    n_cw[n] ~ binomial( n_obs[n],
      inv_logit( link ) .* (1-lapse) + lapse/2);
  }

}'

