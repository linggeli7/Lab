data {
  vector[18] y1;
  vector[6] conc1;
  vector[18] y2;
  vector[6] conc2;
  vector[18] y3;
  vector[6] conc3;
  int n1;
  int n2;
  int n3;
}

parameters {
  // mean
  real<lower=0> E;
  real<lower=0> C1;
  real<lower=0> M1;
  real<lower=0> C2;
  real<lower=0> M2;
  real<lower=0> C3;
  real<lower=0> M3;

  // variance
  real<lower=0> sigma;
}

model {
  vector[6] mu1;
  vector[6] mu2;
  vector[6] mu3;

  // diffused priors
  sigma ~ cauchy(0, 1);
  C1 ~ lognormal(3.9, 0.2);
  M1 ~ lognormal(0, 0.4);
  C2 ~ lognormal(3.9, 0.2);
  M2 ~ lognormal(0, 0.4);
  C3 ~ lognormal(3.9, 0.2);
  M3 ~ lognormal(0, 0.4);
  E ~ lognormal(4.6, 0.1);
  
  for (i in 1:n1) {
    mu1[i] = E/(1 + (conc1[i]/C1)^M1);
    for (j in 1:3) {
      y1[(i - 1) * 3 + j] ~ normal(mu1[i], sigma);
    }
  }
  
  for (i in 1:n2) {
    mu2[i] = E/(1 + (conc2[i]/C2)^M2);
    for (j in 1:3) {
      y2[(i - 1) * 3 + j] ~ normal(mu2[i], sigma);
    }
  }
  
  for (i in 1:n3) {
    mu3[i] = E/(1 + (conc3[i]/C3)^M3);
    for (j in 1:3) {
      y3[(i - 1) * 3 + j] ~ normal(mu3[i], sigma);
    }
  }
}

generated quantities{
  real IC20_1;
  real IC80_1;
  real IC20_2;
  real IC80_2;
  real IC20_3;
  real IC80_3;
  real L50;
  real L20;
  real L80;

  IC20_1 = C1 * 0.25^(1/M1);
  IC80_1 = C1 * 4.0^(1/M1);
  IC20_2 = C2 * 0.25^(1/M2);
  IC80_2 = C2 * 4.0^(1/M2);
  IC20_3 = C3 * 0.25^(1/M3);
  IC80_3 = C3 * 4.0^(1/M3);
  
  L50 = 0.5 * C3/C1 + 0.5 * C3/C2;
  L20 = 0.5 * IC20_3/IC20_1 + 0.5 * IC20_3/IC20_2;
  L80 = 0.5 * IC80_3/IC80_1 + 0.5 * IC80_3/IC80_2;
}

