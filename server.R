library(shiny)
library(ggplot2)

emax <- function(C, E0, IC50, M) {
  return(E0/(1+(C/IC50)^M))
}

generateY <- function(X, E0, IC50, M, l, s) {
  EY <- sapply(X, emax, E0, IC50, M)
  epsilon <- rnorm(length(X), 0, s)
  Y <- EY+epsilon*(EY^l)
  return(Y)
}

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
  names(goodstuff) <- c('m_hat', 'CI_m', 'IC50_hat', 'CI_IC50')
  return(goodstuff)
}

shinyServer(function(input, output) {
   
  output$distPlot <- renderPlot({
    X <- seq(0, input$max, 1)
    EY <- sapply(X, emax, input$E0, input$IC50, input$m)
    if (input$fa) {
      qplot(X, (input$E0-EY)/input$E0, geom='path', lwd=I(1), col=I('red'), 
            xlab='Concentration', ylab='Fraction affected', ylim=c(0, 1))
    } else {
      qplot(X, EY, geom='path', lwd=I(1), col=I('red'),
            xlab='Concentration', ylab='Cell count', ylim=c(0, input$E0))
    }
  })
  
  dat <- reactiveValues(X = c(), Y = c())
  
  observeEvent(input$add, {
    current <- dat$X
    dat$X <- c(current, rep(input$conc, input$k))
  })
  
  observeEvent(input$clear, {
    dat$X <- c()
  })
  
  output$wells <- renderText(dat$X)
  
  observeEvent(input$cook, {
    Y <- generateY(dat$X, input$E0, input$IC50, input$m, input$lambda, input$sigma)
    dat$Y <- Y
    X <- dat$X
    output$sample <- renderPlot({
      qplot(X, Y, xlab='Concentration', ylab='Observed cell count')
    })
  })
  
  observeEvent(input$crunch, {
    treatment <- dat$Y[dat$X != 0]
    control <- dat$Y[dat$X == 0]
    Xt <- dat$X[dat$X != 0]
    Yt <- (mean(control)-treatment)/treatment
    Yt[Yt < 0] <- NA
    fit <- lm(log(Yt)~log(Xt))
    
    output$trans <- renderPlot({
      qplot(log(Xt), log(Yt), xlab='log(Concentration)', ylab='Fa/(1-Fa)') + geom_smooth(method='lm')
    })
    
    output$conf <- renderPrint({
      loglinear_CI(fit)
    })
    
    output$summary <- renderPrint({
      summary(fit)$r.squared
    })
  })
  
  observeEvent(input$run, {
    CI <- matrix(nrow=500, ncol=2)
    for (i in 1:500) {
      Y <- generateY(dat$X, input$E0, input$IC50, input$m, input$lambda, input$sigma)
      treatment <- Y[dat$X != 0]
      control <- Y[dat$X == 0]
      Xt <- dat$X[dat$X != 0]
      Yt <- (mean(control)-treatment)/treatment
      Yt[Yt < 0] <- NA
      fit <- lm(log(Yt)~log(Xt))
      CI[i, ] <- loglinear_CI(fit)$CI_IC50
    }
    goodstuff <- c(mean(CI[ ,1]<input$IC50 & CI[ ,2]>input$IC50), (mean(CI[ ,2]-CI[ ,1])/3.92)^2)
    names(goodstuff) <- c('Coverage', 'Variance')
    output$performance <- renderPrint({
      goodstuff
    })
    output$density <- renderPlot({
      hat <- 0.5*(CI[ ,1]+CI[ ,2])
      cover <- CI[ ,1]<input$IC50 & CI[ ,2]>input$IC50
      dens <- data.frame(hat, cover)
      ggplot(dens, aes(hat, fill=cover)) + geom_dotplot(dotsize=0.4, stackdir="centerwhole")
    })
  })
})
