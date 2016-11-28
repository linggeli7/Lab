library(shiny)
library(ggplot2)
library(gridExtra)
library(grid)

emax <- function(C, E0, IC50, M) {
  return(E0 / (1 + (C / IC50) ^ M))
}

IC <- function(x, IC50, M) {
  return(IC50 * (1 / (1 - x) - 1) ^ (1 / M))
}

generateY <- function(X, E0, IC50, M, l, s) {
  EY <- sapply(X, emax, E0, IC50, M)
  epsilon <- rnorm(length(X), 0, s)
  Y <- EY + epsilon * (EY ^ l)
  return(Y)
}

loglinear_CI <- function(fit) {
  b0 <- as.numeric(fit$coefficients[1])
  b1 <- as.numeric(fit$coefficients[2])
  M_hat <- b1
  se1 <- as.numeric(coef(summary(fit))[, 'Std. Error'][2])
  CI_M <- c(b1 - 1.96 * se1, b1 + 1.96 * se1)
  grad <- as.matrix(c(-1 / b1, b0 / (b1 ^ 2)))
  se2 <- (t(grad) %*% vcov(fit) %*% grad) ^ 0.5
  IC50_hat <- exp(-b0 / b1)
  CI_logIC50 <- c(-b0 / b1 - 1.96 * se2,-b0 / b1 + 1.96 * se2)
  CI_IC50 <- exp(CI_logIC50)
  goodstuff <- list(M_hat, CI_M, IC50_hat, CI_IC50)
  names(goodstuff) <- c('m_hat', 'CI_m', 'IC50_hat', 'CI_IC50')
  return(goodstuff)
}

Loewe_CI <- function(fitA, fitB, fitM, x) {
  b0_A <- as.numeric(fitA$coefficients[1])
  b1_A <- as.numeric(fitA$coefficients[2])
  b0_B <- as.numeric(fitB$coefficients[1])
  b1_B <- as.numeric(fitB$coefficients[2])
  b0_M <- as.numeric(fitM$coefficients[1])
  b1_M <- as.numeric(fitM$coefficients[2])
  Cx <- 1 / (1 / x - 1)
  logIC_A <- (log(Cx) - b0_A) / b1_A
  IC_A <- exp(logIC_A)
  logIC_B <- (log(Cx) - b0_B) / b1_B
  IC_B <- exp(logIC_B)
  logIC_M <- (log(Cx) - b0_M) / b1_M
  IC_M <- exp(logIC_M)
  L <- (IC_M / IC_A + IC_M / IC_B) * 0.5
  grad <-
    as.matrix(c(
      0.5 * IC_M / (IC_A * b1_A),
      0.5 * IC_M * (log(Cx) - b0_A) / (IC_A * b1_A * b1_A),
      0.5 * IC_M / (IC_B * b1_B),
      0.5 * IC_M * (log(Cx) - b0_B) / (IC_B * b1_B * b1_B),-L / b1_M,
      -L * (log(Cx) - b0_M) / (b1_M * b1_M)
    )) / L
  Sigma <- matrix(0, nrow = 6, ncol = 6)
  Sigma[1:2, 1:2] <- vcov(fitA)
  Sigma[3:4, 3:4] <- vcov(fitB)
  Sigma[5:6, 5:6] <- vcov(fitM)
  se <- (t(grad) %*% Sigma %*% grad) ^ 0.5
  CI_logL <- c(log(L) - 1.96 * se, log(L) + 1.96 * se)
  goodstuff <- list(L, exp(CI_logL))
  names(goodstuff) <- c('L_hat', 'CI_L')
  return(goodstuff)
}

shinyServer(function(input, output) {
  output$distPlot <- renderPlot({
    highest <- max(input$IC50A, input$IC50B, input$IC50M)
    X1 <- seq(0, highest * 4, 1)
    EY1 <- sapply(X1, emax, 100, input$IC50A, input$mA)
    if (input$fa) {
      p1 <- qplot(
        X1,
        (100 - EY1) / 100,
        geom = 'path',
        lwd = I(1),
        col = I('red'),
        xlab = 'Concentration',
        ylab = 'Fraction affected',
        ylim = c(0, 1)
      ) + theme(aspect.ratio = 1) + ggtitle(label='Drug A')
    } else {
      p1 <- qplot(
        X1,
        EY1,
        geom = 'path',
        lwd = I(1),
        col = I('red'),
        xlab = 'Concentration',
        ylab = 'Cell count',
        ylim = c(0, 100)
      ) + theme(aspect.ratio = 1) + ggtitle(label='Drug A')
    }
    
    X2 <- seq(0, highest * 4, 1)
    EY2 <- sapply(X2, emax, 100, input$IC50B, input$mB)
    if (input$fa) {
      p2 <- qplot(
        X2,
        (100 - EY2) / 100,
        geom = 'path',
        lwd = I(1),
        col = I('red'),
        xlab = 'Concentration',
        ylab = 'Fraction affected',
        ylim = c(0, 1)
      ) + theme(aspect.ratio = 1) + ggtitle(label='Drug B')
    } else {
      p2 <- qplot(
        X2,
        EY2,
        geom = 'path',
        lwd = I(1),
        col = I('red'),
        xlab = 'Concentration',
        ylab = 'Cell count',
        ylim = c(0, 100)
      ) + theme(aspect.ratio = 1) + ggtitle(label='Drug B')
    }
    
    X3 <- seq(0, highest * 4, 1)
    EY3 <- sapply(X3, emax, 100, input$IC50M, input$mM)
    if (input$fa) {
      p3 <- qplot(
        X3,
        (100 - EY3) / 100,
        geom = 'path',
        lwd = I(1),
        col = I('red'),
        xlab = 'Concentration',
        ylab = 'Fraction affected',
        ylim = c(0, 1)
      ) + theme(aspect.ratio = 1) + ggtitle(label='Drug M')
    } else {
      p3 <- qplot(
        X3,
        EY3,
        geom = 'path',
        lwd = I(1),
        col = I('red'),
        xlab = 'Concentration',
        ylab = 'Cell count',
        ylim = c(0, 100)
      ) + theme(aspect.ratio = 1) + ggtitle(label='Drug M')
    }
    
    grid.arrange(p1, p2, p3, ncol = 3, 
                 bottom=textGrob('Dose-response curves', gp=gpar(fontsize=16, fontface=2)))
  })
  
  output$iso <- renderPlot({
    IC_A1 <- IC(0.2, input$IC50A, input$mA)
    IC_B1 <- IC(0.2, input$IC50B, input$mB)
    IC_M1 <- IC(0.2, input$IC50M, input$mM)
    p1 <-
      qplot(
        c(IC_A1, IC_M1 * 0.5, 0),
        c(0, IC_M1 * 0.5, IC_B1),
        geom = 'point',
        xlab = 'Drug A',
        ylab = 'Drug B'
      ) +
      geom_polygon(
        data = data.frame(x = c(0, 0, IC_A1), y = c(0, IC_B1, 0)),
        aes(x, y),
        fill = 'blue',
        alpha = 0.2
      ) +
      geom_polygon(
        data = data.frame(
          x = c(0, 0, max(IC_A1, IC_M1 * 0.5), max(IC_A1, IC_M1 * 0.5), IC_A1),
          y = c(IC_B1, max(IC_B1, IC_M1 * 0.5), max(IC_B1, IC_M1 *0.5), 0, 0)
        ),
        aes(x, y),
        fill = 'red',
        alpha = 0.2
      ) +
      theme(aspect.ratio = 1) + ggtitle('20% effect')
    
    IC_A2 <- IC(0.5, input$IC50A, input$mA)
    IC_B2 <- IC(0.5, input$IC50B, input$mB)
    IC_M2 <- IC(0.5, input$IC50M, input$mM)
    p2 <-
      qplot(
        c(IC_A2, IC_M2 * 0.5, 0),
        c(0, IC_M2 * 0.5, IC_B2),
        geom = 'point',
        xlab = 'Drug A',
        ylab = 'Drug B'
      ) +
      geom_polygon(
        data = data.frame(x = c(0, 0, IC_A2), y = c(0, IC_B2, 0)),
        aes(x, y),
        fill = 'blue',
        alpha = 0.2
      ) +
      geom_polygon(
        data = data.frame(
          x = c(0, 0, max(IC_A2, IC_M2 * 0.5), max(IC_A2, IC_M2 * 0.5), IC_A2),
          y = c(IC_B2, max(IC_B2, IC_M2 * 0.5), max(IC_B2, IC_M2 *0.5), 0, 0)
        ),
        aes(x, y),
        fill = 'red',
        alpha = 0.2
      ) +
      theme(aspect.ratio = 1) + ggtitle('50% effect')
    
    IC_A3 <- IC(0.8, input$IC50A, input$mA)
    IC_B3 <- IC(0.8, input$IC50B, input$mB)
    IC_M3 <- IC(0.8, input$IC50M, input$mM)
    p3 <-
      qplot(
        c(IC_A3, IC_M3 * 0.5, 0),
        c(0, IC_M3 * 0.5, IC_B3),
        geom = 'point',
        xlab = 'Drug A',
        ylab = 'Drug B'
      ) +
      geom_polygon(
        data = data.frame(x = c(0, 0, IC_A3), y = c(0, IC_B3, 0)),
        aes(x, y),
        fill = 'blue',
        alpha = 0.2
      ) +
      geom_polygon(
        data = data.frame(
          x = c(0, 0, max(IC_A3, IC_M3 * 0.5), max(IC_A3, IC_M3 * 0.5), IC_A3),
          y = c(IC_B3, max(IC_B3, IC_M3 * 0.5), max(IC_B3, IC_M3 *0.5), 0, 0)
        ),
        aes(x, y),
        fill = 'red',
        alpha = 0.2
      ) +
      theme(aspect.ratio = 1) + ggtitle('80% effect')
    
    grid.arrange(p1, p2, p3, ncol = 3,
                 bottom=textGrob('Isobolograms', gp=gpar(fontsize=16, fontface=2)))
  })
  
  dat <-
    reactiveValues(
      X1 = c(),
      X2 = c(),
      X3 = c(),
      Y1 = c(),
      Y2 = c(),
      Y3 = c()
    )
  
  observeEvent(input$add, {
    if (input$drug == 'A') {
      current <- dat$X1
      dat$X1 <- c(current, rep(input$conc, input$k))
    }
    if (input$drug == 'B') {
      current <- dat$X2
      dat$X2 <- c(current, rep(input$conc, input$k))
    }
    if (input$drug == 'Mixture') {
      current <- dat$X3
      dat$X3 <- c(current, rep(input$conc, input$k))
    }
  })
  
  observeEvent(input$clear, {
    dat$X1 <- c()
    dat$X2 <- c()
    dat$X3 <- c()
  })
  
  output$wells1 <- renderText(dat$X1)
  
  output$wells2 <- renderText(dat$X2)
  
  output$wells3 <- renderText(dat$X3)
  
  observeEvent(input$cook, {
    X1 <- dat$X1
    Y1 <-
      generateY(X1, 100, input$IC50A, input$mA, input$lambda, input$sigma)
    dat$Y1 <- Y1
    X2 <- dat$X2
    Y2 <-
      generateY(X2, 100, input$IC50B, input$mB, input$lambda, input$sigma)
    dat$Y2 <- Y2
    X3 <- dat$X3
    Y3 <-
      generateY(X3, 100, input$IC50M, input$mM, input$lambda, input$sigma)
    dat$Y3 <- Y3
    
    output$sample <- renderPlot({
      p1 <-
        qplot(X1, Y1, xlab = 'Concentration', ylab = 'Observed cell count', ylim=c(0, 120)) + 
        theme(aspect.ratio = 1) + ggtitle(label='Drug A')
      p2 <-
        qplot(X2, Y2, xlab = 'Concentration', ylab = 'Observed cell count', ylim=c(0, 120)) + 
        theme(aspect.ratio = 1) + ggtitle(label='Drug B')
      p3 <-
        qplot(X3, Y3, xlab = 'Concentration', ylab = 'Observed cell count', ylim=c(0, 120)) + 
        theme(aspect.ratio = 1) + ggtitle(label='Drug M')
      
      grid.arrange(p1, p2, p3, ncol=3,
                   bottom=textGrob('Simulated data', gp=gpar(fontsize=16, fontface=2)))
    })
  })
  
  # Chou-Talalay fit
  observeEvent(input$crunch, {
    # pool control
    allY <- c(dat$Y1, dat$Y2, dat$Y3)
    allX <- c(dat$X1, dat$X2, dat$X3)
    control <- allY[allX == 0]
    
    # process data for Drug A
    treatment1 <- dat$Y1[dat$X1 != 0]
    Xt1 <- dat$X1[dat$X1 != 0]
    Yt1 <- (mean(control) - treatment1) / treatment1
    Yt1[Yt1 < 0] <- NA
    fit1 <- lm(log(Yt1) ~ log(Xt1))
    
    # process data for Drug B
    treatment2 <- dat$Y2[dat$X2 != 0]
    Xt2 <- dat$X2[dat$X2 != 0]
    Yt2 <- (mean(control) - treatment2) / treatment2
    Yt2[Yt2 < 0] <- NA
    fit2 <- lm(log(Yt2) ~ log(Xt2))
    
    # process data for Drug M
    treatment3 <- dat$Y3[dat$X3 != 0]
    Xt3 <- dat$X3[dat$X3 != 0]
    Yt3 <- (mean(control) - treatment3) / treatment3
    Yt3[Yt3 < 0] <- NA
    fit3 <- lm(log(Yt3) ~ log(Xt3))
    
    # plots
    output$trans <- renderPlot({
      p1 <-
        qplot(log(Xt1), log(Yt1), xlab = 'log(Concentration)', ylab = 'Fa/(1-Fa)') +
        geom_smooth(method = 'lm') + theme(aspect.ratio = 1) + ggtitle(label='Drug A')
      p2 <-
        qplot(log(Xt2), log(Yt2), xlab = 'log(Concentration)', ylab = 'Fa/(1-Fa)') +
        geom_smooth(method = 'lm') + theme(aspect.ratio = 1) + ggtitle(label='Drug B')
      p3 <-
        qplot(log(Xt3), log(Yt3), xlab = 'log(Concentration)', ylab = 'Fa/(1-Fa)') +
        geom_smooth(method = 'lm') + theme(aspect.ratio = 1) + ggtitle(label='Drug M')
      
      grid.arrange(p1, p2, p3, ncol = 3,
                   bottom=textGrob('Transformed data', gp=gpar(fontsize=16, fontface=2)))
    })
    
    # summary
    output$conf1 <- renderPrint({
      c(loglinear_CI(fit1), r_squared = summary(fit1)$r.squared)
    })
    
    output$conf2 <- renderPrint({
      c(loglinear_CI(fit2), r_squared = summary(fit2)$r.squared)
    })
    
    output$conf3 <- renderPrint({
      c(loglinear_CI(fit3), r_squared = summary(fit3)$r.squared)
    })
    
    output$index1 <- renderPrint({
      Loewe_CI(fit1, fit2, fit3, 0.2)
    })
    
    output$index2 <- renderPrint({
      Loewe_CI(fit1, fit2, fit3, 0.5)
    })
    
    output$index3 <- renderPrint({
      Loewe_CI(fit1, fit2, fit3, 0.8)
    })
    
  })
  
  # Asymptotic results
  observeEvent(input$run, {
    CI <- matrix(nrow = 500, ncol = 2)
    for (i in 1:500) {
      Y <-
        generateY(dat$X,
                  input$E0,
                  input$IC50,
                  input$m,
                  input$lambda,
                  input$sigma)
      treatment <- Y[dat$X != 0]
      control <- Y[dat$X == 0]
      Xt <- dat$X[dat$X != 0]
      Yt <- (mean(control) - treatment) / treatment
      Yt[Yt < 0] <- NA
      fit <- lm(log(Yt) ~ log(Xt))
      CI[i,] <- loglinear_CI(fit)$CI_IC50
    }
    goodstuff <-
      c(mean(CI[, 1] < input$IC50 &
               CI[, 2] > input$IC50), (mean(CI[, 2] - CI[, 1]) / 3.92) ^ 2)
    names(goodstuff) <- c('Coverage', 'Variance')
    output$performance <- renderPrint({
      goodstuff
    })
    
    #output$density <- renderPlot({
      #IC50_hat <- 0.5 * (CI[, 1] + CI[, 2])
      #cover <- CI[, 1] < input$IC50 & CI[, 2] > input$IC50
      #dens <- data.frame(IC50_hat, cover)
      #ggplot(dens, aes(IC50_hat, fill = cover)) + 
        #geom_dotplot(dotsize = 0.4, stackdir = "centerwhole")
    #})
  })
})
