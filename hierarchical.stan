data {
  vector[36] y;
  vector[36] conc;
}

parameters {
  real mu_logE0;
  real mu_logIC50;
  real mu_logm;
  real xi_logE0;
  real tau_eta_logE0;
  real xi_logIC50;
  real tau_eta_logIC50;
  real xi_logm;
  real tau_eta_logm;
  real tau_lognormal;
}

transformed parameters {
  real sigma_logE0;
  real tau_logE0;
  real sigma_logIC50;
  real tau_logIC50;
  real sigma_logm;
  real tau_logm;
 
  sigma_logE0 = abs(xi_logE0) * tau_eta_logE0^0.5;
  tau_logE0 = sigma_logE0^(-2);

  sigma_logIC50 = abs(xi_logIC50) * tau_eta_logIC50^0.5;
  tau_logIC50 = sigma_logIC50^(-2);

  sigma_logm = abs(xi_logm) * tau_eta_logm^0.5;
  tau_logm = sigma_logm^(-2);
}

model {
  vector[3] eta_logE0;
  vector[3] eta_logIC50;
  vector[3] eta_logm;
  vector[3] theta_logE0;
  vector[3] logE0;
  vector[3] theta_logIC50;
  vector[3] logIC50;
  vector[3] theta_logm;
  vector[3] logm;

  matrix[3, 4] mu;
  matrix[3, 4] mu_lognormal;

  tau_lognormal ~ gamma(0.01, 0.01);

  mu_logE0 ~ normal(0, 10);
  xi_logE0 ~ normal(0, 0.1);  
  tau_eta_logE0 ~ gamma(0.5, 0.5);

  mu_logIC50 ~ normal(0, 10);
  xi_logIC50 ~ normal(0, 0.1);
  tau_eta_logIC50 ~ gamma(0.5, 0.5);

  mu_logm ~ normal(0, 10);
  xi_logm ~ normal(0, 0.1);
  tau_eta_logm ~ gamma(0.5, 0.5);

  for (i in 1:3) {
    eta_logE0[i] ~ normal(0, tau_eta_logE0);
    eta_logIC50[i] ~ normal(0, tau_eta_logIC50);
    eta_logm[i] ~ normal(0, tau_eta_logm);
  }

  for (i in 1:3) {
    theta_logE0[i] = mu_logE0 + xi_logE0 * eta_logE0[i];
    logE0[i] ~ normal(theta_logE0[i], tau_logE0);
  }

  for (i in 1:3) {
    theta_logIC50[i] = mu_logIC50 + xi_logIC50 * eta_logIC50[i];
    logIC50[i] ~ normal(theta_logIC50[i], tau_logIC50);
  }

  for (i in 1:3) {
    theta_logm[i] = mu_logm + xi_logm * eta_logm[i];
    logm[i] ~ normal(theta_logm[i], tau_logm);
  }

  for (i in 1:3) {
    for (j in 1:4) {
      mu[i, j] = exp(logE0[i])/(1 + (conc[j]/exp(logIC50[i]))^(exp(logm[i])));
      mu_lognormal[i, j] = log(mu[i, j]) - (1/(2 * tau_lognormal));
      for (k in 1:3) {
        y[(i - 1) * 12 + (j -1) * 3 + k] ~ lognormal(mu_lognormal[i, j], tau_lognormal);
      }
    }
  }
}
