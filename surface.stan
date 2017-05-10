data {
  vector[48] y;
  vector[48] conc1;
  vector[48] conc2;
}

parameters {
  // mean
  real<lower=0> E0;
  real b11;
  real b12;
  real b13;
  real b21;
  real b22;
  real b23;
  
  // variance
  real<lower=0> sigma_lognormal;
}

model {
  vector[48] mu;
  vector[48] mu_lognormal;
  real IC50;
  real M;
  real conc;
  
  // Diffused variance prior
  sigma_lognormal ~ cauchy(0, 1);

  E0 ~ lognormal(4.66, 0.08);
  
  // Control
  for (i in 1:3) {
    mu[i] = E0;
    mu_lognormal[i] = log(mu[i]) - 0.5 * (sigma_lognormal^2);
    y[i] ~ lognormal(mu_lognormal[i], sigma_lognormal);
  }
  
  // Drug
  for (i in 4:48) {
    IC50 = exp(b11 * conc1[i] + b12 * conc2[i] + b13 * conc1[i] * conc2[i]);
    M = exp(b21 * conc1[i] + b22 * conc2[i] + b23 * conc1[i] * conc2[i]);
    conc = conc1[i] + conc2[i];
    mu[i] = E0/(1 + (conc/IC50)^M);
    mu_lognormal[i] = log(mu[i]) - 0.5 * (sigma_lognormal^2);
    y[i] ~ lognormal(mu_lognormal[i], sigma_lognormal);
  }
}