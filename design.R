emax <- function(C, E0, IC50, M) {
  return(E0/(1+(C/IC50)^M))
}

makeX <- function(conc, k) {
  X <- rep(conc, k)
  X <- sort(X)
  return(X)
}

generateY <- function(X, E0, IC50, M, l, s) {
  EY <- sapply(X, emax, E0, IC50, M)
  epsilon <- rnorm(length(X), 0, s)
  Y <- EY+epsilon*(EY^l)
  return(Y)
}

prep <- function(Y, E0) {
  ratio <- (E0-Y)/Y
  ratio[ratio<0] <- NA
  return(ratio)  
}

# current design
X1 <- makeX(c(10, 20, 40, 80, 160), 3)
Y1 <- generateY(X1, 100, 50, 0.6, 3)
Yt1 <- prep(Y1, 100)

qplot(X1, Y1, xlab='Concentration', ylab='Observed')
qplot(log(X1), log(Yt1), xlab='log(Conc)', ylab='Fa/(1-Fa)') + geom_smooth(method='lm')

fit1 <- lm(log(Yt1)~log(X1))
summary(fit1)
loglinear_CI(fit1)

# new design
X2 <- makeX(c(10, 40, 160), 5)
Y2 <- generateY(X2, 100, 50, 0.6, 3)
Yt2 <- prep(Y2, 100)

qplot(X2, Y2, xlab='Concentration', ylab='Observed')
qplot(log(X2), log(Yt2), xlab='log(Conc)', ylab='Fa/(1-Fa)') + geom_smooth(method='lm')

fit2 <- lm(log(Yt2)~log(X2))
summary(fit2)
loglinear_CI(fit2)

# compare asymptotic efficiency
n <- 1000
s <- 2
CI1 <- matrix(nrow=n, ncol=2)
CI2 <- matrix(nrow=n, ncol=2)
X1 <- makeX(c(10, 20, 40, 80, 160), 4)
X2 <- makeX(c(20, 40, 80, 160), 4)
for (i in 1:n) {
  e <- rnorm(8, 0, s)
  control <- 100+e
  Y1 <- generateY(X1, 100, 50, 0.6, s)
  Yt1 <- prep(Y1, mean(control[1:4]))
  fit1 <- lm(log(Yt1)~log(X1))
  CI1[i, ] <- loglinear_CI(fit1)$CI_IC50
  Y2 <- generateY(X2, 100, 50, 0.6, s)
  Yt2 <- prep(Y2, mean(control))
  fit2 <- lm(log(Yt2)~log(X2))
  CI2[i, ] <- loglinear_CI(fit2)$CI_IC50
}

# coverage
mean(CI1[ ,1]<50 & CI1[ ,2]>50)
mean(CI2[ ,1]<50 & CI2[ ,2]>50)

# margin of error
mean(CI1[ ,2]-CI1[ ,1])
mean(CI2[ ,2]-CI2[ ,1])

library(nlme)
conc <- makeX(c(0, 10, 20, 40, 80, 160), 3)
cell <- generateY(conc, 100, 50, 0.6, 1, 0.1)
plot(cell~conc)

Plate <- data.frame(conc, cell)
Plate$ID <- 1

# constant variance
fit1 <- nls(cell~E0/(1+(conc/IC50)^m), data=Plate, start=list(E0=100, IC50=50, m=1))
summary(fit1)

# power of mean
fit2 <- nlme(cell~E0/(1+(conc/IC50)^m),
            fixed=E0+IC50+m~1,
            groups=~ID,
            data=Plate, 
            start=c(100, 50, 1),
            weights=varPower(form=~cell))

summary(fit2)

# sanity check
fit3 <- nlme(cell~E0/(1+(conc/IC50)^m),
             fixed=E0+IC50+m~1,
             groups=~ID,
             data=Plate, 
             start=c(100, 50, 1))

summary(fit3)

# lognormal
library(ggplot2)
P <- rlnorm(2000, 3, 0.01)
qplot(P, geom="density", fill=I('#009E73'))

# more simulations
n <- 500
s <- 2
CI1 <- matrix(nrow=n, ncol=2)
CI2 <- matrix(nrow=n, ncol=2)
conc <- makeX(c(0, 10, 20, 40, 80, 160), 4)
cell <- generateY(conc, 100, 50, 0.6, 1, 0.1)
for (i in 1:n) {
  cell <- generateY(conc, 100, 50, 0.6, 1, 0.1)
  
  Yt1 <- prep(Y1, mean(control[1:4]))
  fit1 <- lm(log(Yt1)~log(X1))
  CI1[i, ] <- loglinear_CI(fit1)$CI_IC50
  Y2 <- generateY(X2, 100, 50, 0.6, s)
  Yt2 <- prep(Y2, mean(control))
  fit2 <- lm(log(Yt2)~log(X2))
  CI2[i, ] <- loglinear_CI(fit2)$CI_IC50
}

# coverage
mean(CI1[ ,1]<50 & CI1[ ,2]>50)
mean(CI2[ ,1]<50 & CI2[ ,2]>50)

# margin of error
mean(CI1[ ,2]-CI1[ ,1])
mean(CI2[ ,2]-CI2[ ,1])
summary(P)