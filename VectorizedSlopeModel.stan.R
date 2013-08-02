splits <- c("subject", "content",
            "displacement",
            "target_number_all",
            "target_number_shown",
            "spacing", "eccentricity", "bias")

relevant <- splits %v% c("n_cw", "n_obs")

lapse_limit <- 0.05

stan_format <- mkchain(
    subset(select=relevant),
    as.list,
    factorify,
    put(names(.), gsub('\\.', '_', names(.))),
    within({
      N <- length(subject_ix)
      lapse_limit <- lapse_limit
    }))

model_code <- '
data {
  int<lower=0> N;
  vector[N] target_number_all;
  vector[N] target_number_shown;
  vector[N] displacement;
  vector[N] content;
  vector[N] spacing;
  vector[N] eccentricity;
  int n_cw[N];
  int n_obs[N];
  real lapse_limit;
}

transformed data {
  vector[N] frac_spacing;
  vector[N] ones;
  vector[N] content_nonlinear;
  frac_spacing <- (2 * pi() * target_number_all);
  for (n in 1:N) {
    ones[n] <- 1.0;
    content_nonlinear[n] <- content[n] * abs(content[n]);
  }
}

parameters {
  real beta_dx;
  real<lower=0, upper=lapse_limit> lapse;
  real bias;

  real<lower=0,upper=2*pi()> cs;
  real summation;
  real repulsion;
  real nonlinearity;
}

model {
  vector[N] crowdedness;
  vector[N] link_displacement;
  vector[N] link_repulsion;
  vector[N] link_summation;
  vector[N] link;
  vector[N] response;

  crowdedness <- 2 - (ones*2) ./ (1+exp((ones*cs) ./ frac_spacing));
  link_displacement <- beta_dx * displacement .* crowdedness;
  link_repulsion <- (repulsion * content + nonlinearity * content_nonlinear);
  link_summation <- summation * content .* target_number_shown;
  link <- bias + link_displacement + link_repulsion + link_summation;
  response <- (ones ./ (1+ exp(-link))) * (1-lapse) + lapse/2;
  n_cw ~ binomial( n_obs, response );
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
