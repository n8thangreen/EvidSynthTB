// estimate ltbi prevalence

functions {
#include /include/distributions.stan
}

data {
  int<lower=0> ltbi;        // counts positive
  int<lower=1> n_ethn;      // number of ethnic groups
  int<lower=1> N[n_ethn];   // subsample sizes
}

parameters {
  real alpha;
  real[n_ethn] b_ethn;
}

transformed parameters {
  real mu;
}

model {
  // priors
  alpha ~ normal(0, 10);

  // fixed effect
  for (i in 1:n_ethn) {
    b_ethn[i] ~ normal(0, 10);
  }

  // incidence model (ltbi prevalence)
  for (i in 1:n_ethn) {
    y[i] ~ binomial(mu[i], N[i]);
    logit(mu[i]) = alpha + b_ethn[i]
  }
}

generated quantities {
}
