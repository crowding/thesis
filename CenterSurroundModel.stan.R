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
  int target_number_all[N];
  int target_number_shown[N];
  real displacement[N];
  real content[N];
  real spacing[N];
  real eccentricity[N];
  int n_cw[N];
  int n_obs[N];
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
  real bias;

  real<lower=0,upper=2*pi()> cs;
  real<lower=0> summation;
  real<upper=0> repulsion;
  real<lower=0, upper=pi()> summation_width;
  real<lower=0, upper=pi()> repulsion_width;
  //real center_nonlinearity;
}

model {
  real crowdedness;
  real link_displacement;
  real link_summation;
  real link_repulsion;
  real link;
  real response;

  for (n in 1:N) {
    link_summation <- 0;
    link_repulsion <- 0;
    crowdedness <- 2 - 2/(1+exp(-cs/frac_spacing[n]));
    link_displacement <- beta_dx * displacement[n] * crowdedness;

    //compute the summation by Gaussian center-surrounds
    for (t in 1:target_number_shown[n]) {
        real loc_t;
        loc_t <- spacing[n] * t;
        for (q in 1:target_number_shown[n]) {
            real loc_q;
            loc_q <- spacing[n] * q;
            link_summation <-
               (link_summation + repulsion *
                   exp(normal_log(loc_t, loc_q, summation_width)));
            link_repulsion <-
               (link_repulsion + repulsion *
                   exp(normal_log(loc_t, loc_q, repulsion_width)));
        }
    }
    link <- bias + link_displacement + link_summation;
    response <- inv_logit( link ) .* (1-lapse) + lapse/2;
    n_cw[n] ~ binomial( n_obs[n], response);
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
