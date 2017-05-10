library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## parameter specification
E0 <- 100
IC50 <- 50
M <- 1.5

## design search space
dose <- seq(20, 200, 20)
design <- matrix(nrow=400, ncol=3)
n <- 1
for (i in 1:10) {
  for (j in 1:10) {
    if (i < j) {
      design[n, 1] <- 0
      design[n, 2] <- dose[i]
      design[n, 3] <- dose[j]
      n <- n + 1
    }
  }
}
design <- design[1:(n-1), ]

## compile model
sm <- stan_model(file='/Users/linggeli/Documents/2016/Drug Synergy/Dose.stan')

results <- rep(NA, 45)
## loop through design space
for (i in 1:45) {
  conc <- design[i, ]
  X <- makeX(conc, 3)
  candy <- rep(NA, 100)
  ## sample data and fit model
  for (j in 1:100) {
    y <- generateY(X, E0, IC50, M, 0, 4)
    fit <- sampling(sm, data=list(y=y, conc=conc), chains=1, show_messages=FALSE)
    #distE <- extract(fit, c('E'))$E
    distC <- extract(fit, c('C'))$C
    #distM <- extract(fit, c('M'))$M
    distIC20 <- extract(fit, c('IC20'))$IC20
    distIC80 <- extract(fit, c('IC80'))$IC80
    ## D-optimality
    candy[j] <- var(distC) + var(distIC20) + var(distIC80)
  }
  results[i] <- mean(candy)
}

design[which(results==min(results)), ]

## two D-optimal support points are 20 and 80 when shape parameter is 1.5
## then sequentially they tend to go up due to gain in shape
## they become 20 and 140 when shape parameter is 0.5
## now with concentrated priors 20 and 120
## A-optimal are 40 and 100

## add dosing points
dose <- seq(0, 300, 20)
results <- rep(NA, 16)
for (i in 1:16) {
  conc <- c(0, 20, 120, 120, 120, dose[i])
  X <- makeX(conc, 3)
  candy <- rep(NA, 100)
  ## sample data and fit model
  for (j in 1:100) {
    y <- generateY(X, E0, IC50, M, 0, 4)
    fit <- sampling(sm, data=list(y=y, conc=conc), chains=1, show_messages=FALSE)
    #distE <- extract(fit, c('E'))$E
    distC <- extract(fit, c('C'))$C
    #distM <- extract(fit, c('M'))$M
    distIC20 <- extract(fit, c('IC20'))$IC20
    distIC80 <- extract(fit, c('IC80'))$IC80
    ## D-optimality
    #candy[j] <- var(distC) + var(distM) + var(distE)
    candy[j] <- var(distC) + var(distIC20) + var(distIC80)
  }
  results[i] <- mean(candy)
}

results

## optimal design for synergy
dose <- seq(20, 200, 40)
design <- matrix(nrow=125, ncol=3)
n <- 1
for (i in 1:5) {
  for (j in 1:5) {
    for (k in 1:5) {
      design[n, 1] <- dose[i]
      design[n, 2] <- dose[j]
      design[n, 3] <- dose[k]
      n <- n + 1
    }
  }
}

dose <- seq(10, 160, 10)
n <- length(dose)
results <- rep(NA, n)
for (i in 1:n) {
  #conc1 <- c(0, 100, 30, 70, 110, 140)
  conc1 <- c(0, 20, 40, 120, 160, 140)
  #conc1 <- c(0, 60, 40, 160, 180, 200)
  X1 <- makeX(conc1, 3)
  #conc2 <- c(0, 100, 30, 20, 160, dose[i])
  conc2 <- c(0, 20, 40, 120, 160, 140)
  #conc2 <- c(0, 60, 60, 140, 180, 180)
  X2 <- makeX(conc2, 3)
  #conc3 <- c(0, 100, 20, 40, 140, dose[i])
  #conc3 <- c(0, 20, 40, 70, 80, 80)
  conc3 <- c(0, 20, 40, 80, 100, 120)
  X3 <- makeX(conc3, 3)
  candy <- rep(NA, 100)
  for (j in 1:100) {
    y1 <- generateY(X1, 100, 30, 1.5, 0, 4)
    y2 <- generateY(X2, 100, 35, 1.5, 0, 4)
    y3 <- generateY(X3, 100, 25, 1.5, 0, 4)
    fit <- sampling(sm, data=list(y1=y1, conc1=conc1, y2=y2, conc2=conc2, y3=y3, conc3=conc3, n1=6, n2=6, n3=6), 
                  chains=1, show_messages=FALSE)
    distL50 <- extract(fit, c('L50'))$L50
    distL20 <- extract(fit, c('L20'))$L20
    distL80 <- extract(fit, c('L80'))$L80
    candy[j] <- var(distL50) + var(distL20) + var(distL80)
  }
  results[i] <- mean(candy)
}
results
plot(dose, results)
dose[which(results==min(results))]