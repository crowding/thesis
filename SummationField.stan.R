filter_data <- mkchain(
    subset(exp_type %in% c("spacing", "content", "numdensity"))
    , match_df(., subset(count(., "subject"), freq>2000), on="subject")
    )

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
  real frac_extent[N];
  for (n in 1:N) {
      frac_spacing[n] <- 2*pi()./target_number_all[n];
      frac_extent[n] <- spacing[n] * target_number_shown[n]/eccentricity[n];
  }
}

parameters {
  real beta_dx;
  real<lower=0, upper=lapse_limit> lapse;
  real bias;

  real<lower=0,upper=2*pi()> cs;
  real<lower=0, upper=2*pi()> field;
  real summation;
  //put limits on because unidentifiable if no number/density experiment
  real<lower=-1, upper=1> inside_outside;
  real repulsion;
  real nonlinearity;
}

model {
  real crowdedness;
  real link_displacement;
  real link_repulsion;
  real link_summation_inside;
  real link_summation_outside;
  real link_summation;
  real link;
  real content_in_field;
  real content_outside_field;


  //lapse ~ beta(1.5, 40); //what the shit, this explodes the bias

  for (n in 1:N) {
    real frac_covered;

    frac_covered <- fmin(frac_extent[n], field);
    content_in_field <- content[n] / spacing[n] * frac_covered;
    content_outside_field <- content[n] / spacing[n] * (frac_extent[n]-frac_covered);
    crowdedness <- 2 - 2/(1+exp(-cs/frac_spacing[n]));
    link_displacement <- beta_dx * displacement[n] * crowdedness;
    link_repulsion <- (repulsion * content[n]
                       + nonlinearity * (content[n] * abs(content[n])));
    link_summation_inside <- content_in_field * (summation + inside_outside);
    link_summation_outside <- content_outside_field * (summation - inside_outside);
    link_summation <- link_summation_inside + link_summation_outside;
    link <- bias + link_displacement + link_repulsion + link_summation;
    n_cw[n] ~ binomial( n_obs[n],
      inv_logit( link ) .* (1-lapse) + lapse/2);
  }

}'

stan_predict <- mkchain[., coefs](
     mutate(frac_spacing = 2*pi/target_number_all,
            frac_extent = spacing * target_number_shown / eccentricity,
)
   , with(coefs, summarize(
       .
     , frac_covered = pmin(frac_extent, field)
     , content_in_field = content / spacing * frac_covered
     , content_outside_field = content / spacing * (frac_extent-frac_covered)
     , link_displacement = (beta_dx * displacement
                            * (2 - 2/(1+exp(-cs/frac_spacing))))
     , link_repulsion = (repulsion * content
                         + nonlinearity * (content * abs(content)))
     , link_summation_inside = content_in_field * (summation + inside_outside)
     , link_summation_outside = content_outside_field * (summation - inside_outside)
     , link_summation = link_summation_inside + link_summation_outside
     , link = bias + link_displacement + link_repulsion + link_summation
     , response = plogis(link) * (1-lapse) + lapse/2
     )))
