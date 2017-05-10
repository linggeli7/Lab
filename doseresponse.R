library(ggplot2)
library(grid)
library(gridExtra)

################################# Dose-response curves #################################

# Emax model
emax <- function(C, E0, IC50, M) {
  return(E0/(1+(C/IC50)^M))
}

# Drug A similar to Rosner et al
E0 <- 100
IC50_A <- 40
M_A <- 0.6

C <- seq(0, 120, 0.1)
EY_A <- sapply(C, emax, E0, IC50_A, M_A)
qplot(C, EY_A, geom='path', xlab='Concentration', ylab='Cell count', main='Drug A', ylim=c(0,100))

# Drug B with similar shape
IC50_B <- 100
M_B <- 0.5

C <- seq(0, 512, 0.1)
EY_B <- sapply(C, emax, 3, IC50_B, M_B)
qplot(C, (1-EY_B/E0), geom='path', xlab='Concentration', ylab='Fraction Affected', main='Drug B')

# Mixture
IC50_M <- 50
M_M <- 0.4

C <- seq(0, 512, 0.1)
EY_M <- sapply(C, emax, E0=3, IC50_M, M_M)
qplot(C, (1-EY_M/E0), geom='path', xlab='Concentration', ylab='Fraction Affected', main='Mixture')

# Different shapes
C <- seq(0, 256, 0.1)
EY <- sapply(C, emax, 3, 40, 0.1)
qplot(C, (1-EY/E0), geom='path', xlab='Concentration', ylab='Fraction Affected', main='Steep')

EY <- sapply(C, emax, 3, 100, 3)
qplot(C, (1-EY/E0), geom='path', xlab='Concentration', ylab='Fraction Affected', main='Sigmoidal')

# IC calculation
IC <- function(x, IC50, M) {
  return(IC50*(1/(1-x)-1)^(1/M))
}

# Isobologram
Isobole <- function(x) {
  IC_A <- IC(x, IC50_A, M_A)
  IC_B <- IC(x, IC50_B, M_B)
  IC_M <- IC(x, IC50_M, M_M)
  qplot(c(IC_A, IC_M*0.5, 0), c(0, IC_M*0.5, IC_B), geom='point', xlab='Drug A', ylab='Drug B') +
    geom_polygon(data=data.frame(x=c(0,0,IC_A), y=c(0,IC_B,0)), aes(x, y), fill='blue', alpha=0.2) + 
    geom_polygon(data=data.frame(x=c(0,IC_A,max(IC_A,IC_M*0.5),max(IC_A,IC_M*0.5)), 
                                 y=c(IC_B,0,0,max(IC_B,IC_M*0.5))), aes(x, y), fill='red', alpha=0.2)
}

# Loewe index
Loewe <- function(x) {
  IC_A <- IC(x, IC50_A, M_A)
  IC_B <- IC(x, IC50_B, M_B)
  IC_M <- IC(x, IC50_M, M_M)
  return(0.5*(IC_M/IC_A + IC_M/IC_B))
}

################################# Data generation #################################

# E0 distribution
# It should not change for different drugs?
E <- rlnorm(1000, 1.09, 0.04)
qplot(E, geom="density", fill=1)

# Drug A
IC50_A <- 40
M_A <- 0.6

y_A <- rep(NA, 72)
x_A <- rep(NA, 72)
# Serial 2-fold dilution
conc_A <- c(0, 4, 8, 16, 32, 64, 128, 256)
# 3 experiments (i = 3) 
# 8 concentration levels (j = 8) 
# 3 duplicates (k = 3)
for (i in 1:3) {
  E0_i <- rlnorm(1, 1.09, 0.04)
  for (j in 1:8) {
    C_ij <- conc_A[j]
    mu_ij <- emax(C_ij, E0_i, IC50_A, M_A)
    for (k in 1:3) {
      # Heteroskedasticity
      epsilon_ijk <- rnorm(1, 0, 0.2)
      Y_ijk <- mu_ij + mu_ij*epsilon_ijk
      id <- 24*(i-1) + 3*(j-1) + k
      y_A[id] <- max(Y_ijk, 0)
      x_A[id] <- C_ij
    }
  }
}

# Scatterplot of raw data plus smoother
C <- seq(0, 256, 0.1)
qplot(x_A, y_A, xlab='Concentration', ylab='Observed') + geom_smooth(size = 0.6) + 
  geom_path(data=data.frame(C, EY_A), aes(x=C, y=EY_A), colour='red', size = 0.8)

# Drug B
IC50_B <- 100
M_B <- 0.5

y_B <- rep(NA, 72)
x_B <- rep(NA, 72)
# Serial 2-fold dilution
conc_B <- c(0, 8, 16, 32, 64, 128, 256, 512)
for (i in 1:3) {
  E0_i <- rlnorm(1, 1.09, 0.04)
  for (j in 1:8) {
    C_ij <- conc_B[j]
    mu_ij <- emax(C_ij, E0_i, IC50_B, M_B)
    for (k in 1:3) {
      # Heteroskedasticity
      epsilon_ijk <- rnorm(1, 0, 0.2)
      Y_ijk <- mu_ij + mu_ij*epsilon_ijk
      id <- 24*(i-1) + 3*(j-1) + k
      y_B[id] <- max(Y_ijk, 0)
      x_B[id] <- C_ij
    }
  }
}

C <- seq(0, 512, 0.1)
qplot(x_B, y_B, xlab='Concentration', ylab='Observed') + geom_smooth(size = 0.6) + 
  geom_path(data=data.frame(C, EY_B), aes(x=C, y=EY_B), colour='red', size = 0.8)

# Mixture
IC50_M <- 50
M_M <- 0.4

y_M <- rep(NA, 72)
x_M <- rep(NA, 72)
# Serial 2-fold dilution
conc_M <- c(0, 8, 16, 32, 64, 128, 256, 512)
for (i in 1:3) {
  E0_i <- rlnorm(1, 1.09, 0.04)
  for (j in 1:8) {
    C_ij <- conc_M[j]
    mu_ij <- emax(C_ij, E0_i, IC50_M, M_M)
    for (k in 1:3) {
      # Heteroskedasticity
      epsilon_ijk <- rnorm(1, 0, 0.2)
      Y_ijk <- mu_ij + mu_ij*epsilon_ijk
      id <- 24*(i-1) + 3*(j-1) + k
      y_M[id] <- max(Y_ijk, 0)
      x_M[id] <- C_ij
    }
  }
}

C <- seq(0, 512, 0.1)
qplot(x_M, y_M, xlab='Concentration', ylab='Observed') + geom_smooth(size = 0.6) + 
  geom_path(data=data.frame(C, EY_M), aes(x=C, y=EY_M), colour='red', size = 0.8)

################################# Hierarchical model #################################

# Diffused priors
IC50 <- rlnorm(2000, 4.0, 0.4)
qplot(IC50, geom="density", fill=I('red'), main='Drug A')

IC50 <- rlnorm(2000, 4.6, 0.5)
qplot(IC50, geom="density", fill=1, main='Drug B')

IC50 <- rlnorm(2000, 4.2, 0.5)
qplot(IC50, geom="density", fill=1, main='Mixture')

m <- rlnorm(1000, 0, 0.4)
qplot(m, geom="density", fill=1)

# Fit model in Stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fit <- stan(file='/Users/linggeli/Documents/2016/Drug Synergy/Naive.stan', iter=10000)

# Posteriors
print(fit, c('IC50_A', 'M_A', 'IC50_B', 'M_B', 'IC50_M', 'M_M', 'mu_logE0', 'sigma_logE0'))
stan_dens(fit, c('IC50_A', 'M_A'))

# Loewe index
print(fit, c('L50', 'L20', 'L80'))
stan_dens(fit, c('L50', 'L20', 'L80'))

# Well mixed?
traceplot(fit, c('IC50_A', 'M_A', 'mu_logE0', 'sigma_logE0'), inc_warmup=FALSE)

# Autocorrelation

################################# Median effect model #################################

# Pre-process the data
prep <- function(y, conc) {
  ratio <- rep(NA, 21)
  dose <- rep(NA, 21)
  for (i in 1:3) {
    control <- mean(y[((i-1)*24+1):((i-1)*24+3)])
    for (j in 1:7) {
      # THIS AVERAGING DOES NOT MAKE SENSE
      drug <- mean(y[((i-1)*24+1+j*3):((i-1)*24+3+j*3)])
      ratio[(i-1)*7+j] <- (control-drug)/drug
      dose[(i-1)*7+j] <- conc[j+1]
      if ((control-drug)/control < 0 || (control-drug)/control > 1) {
        print('dammit')
      }
    }
  }
  return(data.frame(ratio, dose))  
}

medA <- prep(y_A, conc_A)

# Log-linear fit
qplot(log(medA$dose), log(medA$ratio)) + geom_smooth(method='lm')
fit_logA <- lm(log(medA$ratio)~log(medA$dose))
summary(fit_logA)

# Confidence intervals
loglinear_CI <- function(fit) {
  b0 <- as.numeric(fit$coefficients[1])
  b1 <- as.numeric(fit$coefficients[2])
  M_hat <- b1
  se1 <- as.numeric(coef(summary(fit))[, 'Std. Error'][2])
  CI_M <- c(b1-1.96*se1, b1+1.96*se1)
  grad <- as.matrix(c(-1/b1, b0/(b1^2)))
  se2 <- (t(grad) %*% vcov(fit) %*% grad)^0.5
  IC50_hat <- exp(-b0/b1)
  CI_logIC50 <- c(-b0/b1-1.96*se2, -b0/b1+1.96*se2)
  CI_IC50 <- exp(CI_logIC50)
  goodstuff <- list(M_hat, CI_M, IC50_hat, CI_IC50)
  names(goodstuff) <- c('M_hat', 'CI_M', 'IC50_hat', 'CI_IC50')
  return(goodstuff)
}

loglinear_CI(fit_logA)

# More drugs
medB <- prep(y_B, conc_B)

qplot(log(medB$dose), log(medB$ratio)) + geom_smooth(method='lm')
fit_logB <- lm(log(medB$ratio)~log(medB$dose))

loglinear_CI(fit_logB)

medM <- prep(y_M, conc_M)

qplot(log(medM$dose), log(medM$ratio)) + geom_smooth(method='lm')
fit_logM <- lm(log(medM$ratio)~log(medM$dose))

loglinear_CI(fit_logM)

# Loewe index confidence intervals not as good
Loewe_CI <- function(fitA, fitB, fitM, x) {
  b0_A <- as.numeric(fitA$coefficients[1])
  b1_A <- as.numeric(fitA$coefficients[2])
  b0_B <- as.numeric(fitB$coefficients[1])
  b1_B <- as.numeric(fitB$coefficients[2])
  b0_M <- as.numeric(fitM$coefficients[1])
  b1_M <- as.numeric(fitM$coefficients[2])
  Cx <- 1/(1/x-1)
  logIC_A <- (log(Cx)-b0_A)/b1_A
  IC_A <- exp(logIC_A)
  logIC_B <- (log(Cx)-b0_B)/b1_B
  IC_B <- exp(logIC_B)
  logIC_M <- (log(Cx)-b0_M)/b1_M
  IC_M <- exp(logIC_M)
  L <- (IC_M/IC_A + IC_M/IC_B) * 0.5
  grad <- as.matrix(c(0.5*IC_M/(IC_A*b1_A), 0.5*IC_M*(log(Cx)-b0_A)/(IC_A*b1_A*b1_A),
                      0.5*IC_M/(IC_B*b1_B), 0.5*IC_M*(log(Cx)-b0_B)/(IC_B*b1_B*b1_B),
                      -L/b1_M, -L*(log(Cx)-b0_M)/(b1_M*b1_M)))/L
  Sigma <- matrix(0, nrow=6, ncol=6)
  Sigma[1:2, 1:2] <- vcov(fitA)
  Sigma[3:4, 3:4] <- vcov(fitB)
  Sigma[5:6, 5:6] <- vcov(fitM)
  se <- (t(grad) %*% Sigma %*% grad)^0.5
  CI_logL <- c(log(L)-1.96*se, log(L)+1.96*se)
  goodstuff <- list(L, exp(CI_logL))
  names(goodstuff) <- c('L_hat', 'CI_L')
  return(goodstuff)
}

Loewe_CI(fit_logA, fit_logB, fit_logM, 0.5)

################################# Experimentation #################################
# No random effects but fewer concentrations
# Less variation with pooled variance
E0 <- 100

################################# Drug A #################################
IC50_A <- 36
M_A <- 0.61

y_A <- rep(NA, 18)
x_A <- rep(NA, 18)
# fewer concentration levels
conc_A <- c(0, 10, 20, 40, 60, 100)
for (j in 1:6) {
  C_j <- conc_A[j]
  mu_j <- emax(C_j, E0, IC50_A, M_A)
  for (k in 1:3) {
    # Heteroskedasticity
    epsilon_jk <- rnorm(1, 0, 0.1)
    Y_jk <- mu_j + mu_j*epsilon_jk
    id <- 3*(j-1) + k
    y_A[id] <- max(Y_jk, 0)
    x_A[id] <- C_j
  }
}

# Scatterplot of raw data
qplot(x_A, y_A, xlab='Concentration', ylab='Observed', ylim=c(0,120)) + 
  geom_path(data=data.frame(C, EY_A), aes(x=C, y=EY_A), colour='blue', size = 0.8)

medA <- prep(y_A, conc_A)

# Log-linear fit
qplot(log(medA$dose), log(medA$ratio)) + geom_smooth(method='lm')
fit_logA <- lm(log(medA$ratio)~log(medA$dose))
summary(fit_logA)
loglinear_CI(fit_logA)

################################# Drug B #################################
IC50_B <- 92
M_B <- 0.47

y_B <- rep(NA, 18)
x_B <- rep(NA, 18)
# fewer concentration levels
conc_B <- c(0, 30, 60, 100, 150, 200)
for (j in 1:6) {
  C_j <- conc_B[j]
  mu_j <- emax(C_j, E0, IC50_B, M_B)
  for (k in 1:3) {
    # Heteroskedasticity
    epsilon_jk <- rnorm(1, 0, 0.04)
    Y_jk <- mu_j + mu_j*epsilon_jk
    id <- 3*(j-1) + k
    y_B[id] <- max(Y_jk, 0)
    x_B[id] <- C_j
  }
}

# HACK for control
y_B[1:3] <- y_A[1:3]

# Scatterplot of raw data
qplot(x_B, y_B, xlab='Concentration', ylab='Observed', ylim=c(0,120))

medB <- prep(y_B, conc_B)

# Log-linear fit
qplot(log(medB$dose), log(medB$ratio)) + geom_smooth(method='lm')
fit_logB <- lm(log(medB$ratio)~log(medB$dose))
summary(fit_logB)
loglinear_CI(fit_logB)

################################# Mixture #################################
IC50_M <- 45
M_M <- 0.42

y_M <- rep(NA, 18)
x_M <- rep(NA, 18)
# fewer concentration levels
conc_M <- c(0, 15, 30, 50, 80, 120)
for (j in 1:6) {
  C_j <- conc_M[j]
  mu_j <- emax(C_j, E0, IC50_M, M_M)
  for (k in 1:3) {
    # Heteroskedasticity
    epsilon_jk <- rnorm(1, 0, 0.04)
    Y_jk <- mu_j + mu_j*epsilon_jk
    id <- 3*(j-1) + k
    y_M[id] <- max(Y_jk, 0)
    x_M[id] <- C_j
  }
}

# HACK for control
y_M[1:3] <- y_A[1:3]

# Scatterplot of raw data
qplot(x_M, y_M, xlab='Concentration', ylab='Observed', ylim=c(0,120))

medM <- prep(y_M, conc_M)

# Log-linear fit
qplot(log(medM$dose), log(medM$ratio)) + geom_smooth(method='lm')
fit_logM <- lm(log(medM$ratio)~log(medM$dose))
summary(fit_logM)
loglinear_CI(fit_logM)

Loewe_CI(fit_logA, fit_logB, fit_logM, 0.2)
Loewe_CI(fit_logA, fit_logB, fit_logM, 0.5)
Loewe_CI(fit_logA, fit_logB, fit_logM, 0.8)

# NO AVERAGING
prep <- function(y, conc) {
  ratio <- rep(NA, 15)
  dose <- rep(NA, 15)
  control <- mean(y[1:3])
  for (j in 1:5) {
    drug <- y[(1+j*3):(3+j*3)]
    ratio[(1+(j-1)*3):(3+(j-1)*3)] <- (control-drug)/drug
    dose[(1+(j-1)*3):(3+(j-1)*3)] <- conc[j+1]
    if ((control-drug)/control < 0 || (control-drug)/control > 1) {
      print('dammit')
    }
  }
  return(data.frame(ratio, dose))  
}

################################# Prior buster #################################
P <- rlnorm(2000, 3.9, 0.2)
qplot(P, geom="density", fill=I('#009E73'))
summary(P)

# Bayesian model
fit <- stan(file='/Users/linggeli/Documents/2016/Drug Synergy/Trial.stan')

# Posteriors
print(fit, c('IC50_A', 'M_A', 'IC50_B', 'M_B', 'IC50_M', 'M_M', 'E0'))
stan_dens(fit, c('IC50_A', 'M_A'))

# Loewe index
print(fit, c('L20', 'L50', 'L80'))
stan_dens(fit, c('L20', 'L50', 'L80'))

################################# Surface #################################
y <- rep(NA, 48)
y[1:18] <- y_A[1:18]
y[19:33] <- y_B[4:18]
y[34:48] <- y_M[4:18]

conc1 <- rep(NA, 48)
conc2 <- rep(NA, 48)

conc1[1:18] <- x_A
conc2[1:18] <- 0
conc1[19:33] <- 0
conc2[19:33] <- x_B[4:18]
conc1[34:48] <- 0.5 * x_M[4:18]
conc2[34:48] <- 0.5 * x_M[4:18]

fit <- stan(file='/Users/linggeli/Documents/2016/Drug Synergy/Surface.stan')

print(fit, c('b11', 'b12', 'b13', 'b21', 'b22', 'b23', 'E0'))

# induced priors

######################## Posterior band #########################
fit <- stan(file='/Users/linggeli/Documents/2016/Drug Synergy/Trial.stan')

# Posteriors
print(fit, c('IC50_A', 'M_A', 'IC50_B', 'M_B', 'IC50_M', 'M_M', 'E0'))
print(fit, c('L50', 'L20', 'L80'))
stan_dens(fit, c('L50', 'L20', 'L80'))

sampleIC <- function(ra, n) {
  index <- sample(seq(1,4000,1), n)
  IC_A <- IC50$IC50_A[index]
  IC_B <- IC50$IC50_B[index]
  IC_M <- IC50$IC50_M[index]
  b1 <- IC_A
  b2 <- IC_B
  b3 <- 4*IC_M-2*(IC_A+IC_B)
  IC <- b1*ra+b2*(1-ra)+b3*ra*(1-ra)
  return(list(x=IC*ra, y=IC*(1-ra)))
}

D <- sampleIC(0.5, 20)
qplot(D$x, D$y, xlim=c(0,100), ylim=c(0,200), alpha=0.02)

A <- rep(10000)
B <- rep(10000)
for (i in 1:200) {
  ra <- (i-1)*0.005
  D <- sampleIC(ra, 50)
  A[((i-1)*50+1):((i-1)*50+50)] <- D$x
  B[((i-1)*50+1):((i-1)*50+50)] <- D$y
}

IC50_A <- 38
IC50_B <- 98
IC50_M <- 46

qplot(A, B, alpha=I(0.2), stroke=I(0)) +
  geom_polygon(data=data.frame(x=c(0,0,IC50_A), y=c(0,IC50_B,0)), aes(x, y), fill='blue', alpha=0.2) + 
  geom_polygon(data=data.frame(x=c(0,0,IC50_A,max(40,IC50_M*0.5),max(40,IC50_M*0.5)), 
                               y=c(100,80,0,0,max(100,IC50_M*0.5))), aes(x, y), fill='red', alpha=0.2)
