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
  real nonlinearity;
}

model {
  for (n in 1:N) {
    real crowdedness;
    real link_displacement;
    real link_summation;
    real link_repulsion;
    real link_nonlinearity;
    real link;
    real response;

    link_summation <- 0;
    link_repulsion <- 0;
    crowdedness <- 2 - 2/(1+exp(-cs/frac_spacing[n]));
    link_displacement <- beta_dx * displacement[n] * crowdedness;

    for (d in 0:target_number_shown[n]-1) {
        real distance;
        real multiplicity;
        //add up center-surround influence between every pair of elements
        //separated by d spaces
        real dd;
        if ( d > target_number_all[n] - d) {
            dd <- target_number_all[n] - d;
        } else {
            dd <- d;
        };
        distance <- frac_spacing[n] * dd;
        if (d == 0) {
            //an element only influences itself once
            multiplicity <- target_number_shown[n];
        } else {
            //separated elements mutually influence, so twice
            multiplicity <- 2 * (target_number_shown[n] - dd);
        }
        link_summation <- link_summation
          + summation * multiplicity * content[n]
            * exp(normal_log(distance, 0, summation_width));
        link_repulsion <- link_repulsion
          + repulsion * multiplicity * content[n]
            * exp(normal_log(distance, 0, repulsion_width));
    }
    //Average the influence over each element
    link_summation <- link_summation / target_number_shown[n];
    link_repulsion <- link_repulsion / target_number_shown[n];
    link_nonlinearity <- link_nonlinearity * content[n] * abs(content[n]);

    link <- bias + link_displacement + link_summation + link_repulsion + link_nonlinearity;
    response <- inv_logit( link ) .* (1-lapse) + lapse/2;
    n_cw[n] ~ binomial( n_obs[n], response);
  }
}'

 stan_predict <- mkchain[., coefs](
     subset(., select=intersect(relevant, names(.)))
   , mutate(frac_spacing = 2*pi/target_number_all)
   , mutate(., .n=1:nrow(.))
   , ddply(.,c("target_number_shown", "target_number_all"), function(group) {
     d <- 0:(group$target_number_shown[1] - 1)
     d <- pmin(d, group$target_number_all[1] - d)
     distance <- d * group$frac_spacing[1]
     multiplicity <- ifelse(
       d==0,
       group$target_number_shown[1], 2*(group$target_number_shown[1] - d))
     raw_summation <- (
         rep(1, length(multiplicity))
         %*% matrix(dnorm(distance, 0,
                          rep(coefs$summation_width, length(distance)))
                    * multiplicity,
                    nrow=length(distance)))
     raw_repulsion <- (
         rep(1, length(multiplicity))
         %*% matrix(dnorm(distance, 0,
                          rep(coefs$repulsion_width, length(distance)))
                    * multiplicity,
                    nrow=length(distance)))
     cbind(group, data.frame(
         link_summation =
             (group$content * coefs$summation
              * as.vector(raw_summation)/group$target_number_shown)
       , link_repulsion =
             (group$content * coefs$repulsion
              * as.vector(raw_repulsion)/group$target_number_shown)
         ))
   })
   , .[order(.$.n), ]
   , with(coefs, summarize(
       .
     , link_displacement = beta_dx * displacement
                             * (2 - 2/(1+exp(-cs/frac_spacing)))
     , link = bias + link_displacement + link_repulsion + link_summation
     , response = plogis(link) * (1-lapse) + lapse/2
       )))
