multiinitialize = function(){
  
  
  ################################# GIBBS SAMPLING  ###################################################
  
  
  
  ## HYPER PRIORS
  ## Hyper parameters of the DP
  shape.alpha <- 2
  rate.alpha <- 1
 
  
  source('rchinese.R')
  alpha  = rgamma(1, shape = shape.alpha, rate = rate.alpha )
  c <-  rchinese(N,alpha)
  f <- table(factor(c, levels = 1:max(c)))
  
  #Sparsity controlling hyperparameter of the BAYESIAN LASSO MODEL
  r =1
  si = 1.78
  
  
  
  ### LETS MAKE A LIST "gmmx" to store parameters/hyperprameters for X and "regy" to store paameters for Regression Y
  ## For the First Data Set
  gmmx1 <- list(0)
  gmmx1$epsilon <-  as.vector(apply(Y1,2,mean))
  gmmx1$W <- diag(diag(cov(Y1)))
  gmmx1$mu <- matrix(data = NA, nrow = K, ncol = D1)
  gmmx1$S <-  array(data = NA, dim =c(K,D1,D1))
  gmmx1$ro <- 0.5
  gmmx1$beta <- D1 +1
  
  
  
  
  regy1 <- list(0)
  regy1$lambda2 <- numeric(K)
  regy1$tau2 = matrix(data = NA, nrow = K, ncol = D1)
  regy1$betahat = matrix(data = NA, nrow = K, ncol = D1)
  regy1$sigma2 <- rep(NA, K)
  regy1$beta0 <- rep(NA, K)
  
  ## For the second data set
  gmmx2 <- list(0)
  gmmx2$epsilon <-  as.vector(apply(Y2,2,mean))
  gmmx2$W <- diag(diag(cov(Y2)))
  gmmx2$mu <- matrix(data = NA, nrow = K, ncol = D2)
  gmmx2$S <-  array(data = NA, dim =c(K,D2,D2))
  gmmx2$ro <- 0.5
  gmmx2$beta <- D2 +1
  
  
  
  regy2 <- list(0)
  regy2$lambda2 <- numeric(K)
  regy2$tau2 = matrix(data = NA, nrow = K, ncol = D2)
  regy2$betahat = matrix(data = NA, nrow = K, ncol = D2)
  regy2$sigma2 <- rep(NA, K)
  regy2$beta0 <- rep(NA, K)
  
  
  ###### To initialize the parameters for all the data sets
  That <-  time
  ####### We can use a simple Linear Model to get some estimates of the variance##########
  Yg <- cbind(Y1,Y2)
  Dg <- (D1 + D2)
  
  ## Fitting a linear model to the whole model
  Ysc <- scale(Yg[1:N,1:Dg], center = TRUE, scale =TRUE)
  lm.data <- lm(time ~ Ysc)
  sig2.dat <-  var(lm.data$residuals)
  
  
  ## Set Some Initial Values for the Cluster Parameters
  source('multiinit.R')
  
  ## For the first data set
  cont1 <- multiinit(Y1,c, gmmx1$beta, gmmx1$W, gmmx1$epsilon, gmmx1$ro, r, si,N,D1, sig2.dat)
  gmmx1$mu <- cont1$mu
  gmmx1$S <- cont1$S
  regy1$lambda2 <- cont1$lambda2
  regy1$tau2 <- cont1$tau2
  regy1$betahat <- cont1$betahat
  regy1$sigma2 <- cont1$sigma2
  regy1$beta0 <- cont1$beta0
  
  ## For the second data set
  cont2 <- multiinit(Y2,c, gmmx1$beta, gmmx2$W, gmmx2$epsilon, gmmx2$ro, r, si,N,D2, sig2.dat)
  gmmx2$mu <- cont2$mu
  gmmx2$S <- cont2$S
  regy2$lambda2 <- cont2$lambda2
  regy2$tau2 <- cont2$tau2
  regy2$betahat <- cont2$betahat
  regy2$sigma2 <- cont2$sigma2
  regy2$beta0 <- cont2$beta0
  
  
  surv.obj <-  Surv(time,censoring)
  ## Initialization part for the parmaters of AFT Model with k-means and Bayesian Lasso and Normal Bayesian Regression
  source('multikmeansBlasso.R')
  km <- multikmeansBlasso(c,Y1,Y2,D1,D2,That,K, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2,surv.obj )
  c <- km$c
  c.kmeans <<- c
  
  gmmx1 <- km$gmmx1
  gmmx2 <- km$gmmx2 
  regy1 <- km$regy1
  regy2 <- km$regy2
  
  
  ## Adjusted Initial Rand INDEX measure
 ## randindexi <<- adjustedRandIndex(c.true,as.factor(c))
  
  ## Initial Likelihood
  source('multilikelihood.R')
  likli.int <<- multiloglikelihood( c,Y1,Y2,D1,D2,That,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  
  
  
  
  
 
  
  assign("time", time, envir = .GlobalEnv)
  assign("r", r, envir = .GlobalEnv)
  assign("si", si, envir = .GlobalEnv)
  assign("shape.alpha", shape.alpha, envir = .GlobalEnv)
  assign("rate.alpha", shape.alpha, envir = .GlobalEnv)
  assign("alpha", alpha, envir = .GlobalEnv)
   assign("That", That, envir = .GlobalEnv)
  assign("c", c, envir = .GlobalEnv)
  assign("gmmx1", gmmx1, envir = .GlobalEnv)
  assign("gmmx2", gmmx2, envir = .GlobalEnv)
  assign("regy1", regy1, envir = .GlobalEnv)
  assign("regy2", regy2, envir = .GlobalEnv)
  assign("sig2.dat", sig2.dat, envir = .GlobalEnv)
  
}
