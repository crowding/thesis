data {
  int<lower=0> N;
  int<lower=0> N_data;
  int i_data[N_data];
  int i_predict[N - N_data];

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
  int n; 

  for (id in 1:N_data) {
    n <- i_data[id];
    crowdedness <- 2 - 2/(1+exp(-cs/frac_spacing[n]));
    link_displacement <- beta_dx * displacement[n] * crowdedness;
    link_repulsion <- (content_repulsion * content[n]
                       + content_nonlinearity * (content[n] * abs(content[n])));
    link_summation <- target_number_all[n] * content[n] * content_global_summation;
    link <- link_displacement + link_repulsion + link_summation;

    n_cw[n] ~ binomial( n_obs[n],
      inv_logit( link + bias ) .* (1-lapse) + lapse/2);
  }
}

generated quantities {
  real predict_link[N - N_data];
  real predict_displacement[N - N_data];
  real predict_repulsion[N - N_data];
  real predict_link[N - N_data];

  for (ip in 1:N - N_data) {
    n <- i_predict[ip];
    crowdedness <- 2 - 2/(1+exp(-cs/frac_spacing[n]));
    predict_link_displacement[ip] <- beta_dx * displacement[n] * crowdedness;
    predict_link_repulsion[ip] <-
        (content_repulsion * content[n]
         + content_nonlinearity * (content[n] * abs(content[n])));
    predict_link_summation[ip] <-
        target_number_all[n] * content[n] * content_global_summation;
    predict_link[ip] <-
        predict_link_displacement[ip] +
            predict_link_repulsion[ip] + predict_link_summation[ip];
  }
}
