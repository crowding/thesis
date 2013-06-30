splits <- c("subject", "content",
            "displacement",
            "target_number_all",
            "target_number_shown",
            "spacing", "eccentricity", "bias")

dsharpness <- 10

relevant <- splits %v% c("n_cw", "n_obs")

stan_format <- mkchain(
    subset(select=relevant),
    as.list,
    factorify,
    put(names(.), gsub('\\.', '_', names(.))),
    within({
      N <- length(subject_ix)
      dsharpness <- dsharpness
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
  real dsharpness;
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
  real link_displacement;
  real link_repulsion;
  real link_summation;
  real link;

  for (n in 1:N) {

    crowdedness <- -log(              //soft min of
          exp(-1*dsharpness)          //sensitivity 1 (uncrowded or)
        + exp(-frac_spacing[n]/cs/2*dsharpness) //sensitivity limited by crowding
    )/dsharpness;

    link_displacement <- beta_dx * displacement[n] * crowdedness;
    link_repulsion <- (content_repulsion * content[n]
                       + content_nonlinearity * (content[n] * abs(content[n])));
    link_summation <- target_number_all[n] * content[n] * content_global_summation;
    link <- bias + link_displacement + link_repulsion + link_summation;
    n_cw[n] ~ binomial( n_obs[n],
      inv_logit( link ) .* (1-lapse) + lapse/2);
  }

}'


stan_predict <- mkchain[., coefs](
    mutate(frac_spacing = 2*pi/target_number_all)
  , with(coefs, summarize(
      .
    , crowdedness= -log(                                #soft min of
                        exp(-1*dsharpness)              #sensitivity 1 (uncrowded) or
                        + exp(-spacing/cs/2*dsharpness) #sensitivity limited by crowding
                        )/dsharpness
    , link_displacement = (beta_dx * displacement
                           * (2 - 2/(1+exp(-cs/frac_spacing))))
    , link_repulsion = (content_repulsion * content
                        + content_nonlinearity * (
                          content * abs(content)))
    , link_summation = (target_number_all
                        * content
                        * content_global_summation)
    , link = bias + link_displacement + link_repulsion + link_summation
    , response = plogis(link) * (1-lapse) + lapse/2
    )))
