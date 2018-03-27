Pinitialize = function(Y,D,N,K,time,censoring){
  
  Time <- cbind(time, censoring) 
  ## HYPER PRIORS
  ## Hyper parameters of the DP
  shape.alpha <- 2
  rate.alpha <- 1
  ## Hyperparameters for the GMM
  beta  = D+1
  ro = 0.5
  
  
  source('rchinese.R')
  alpha  = rgamma(1, shape = shape.alpha, rate = rate.alpha )
  c <-  rchinese(N,alpha)
  f <- table(factor(c, levels = 1:max(c)))
  
  ## Empirical Bayes Estimate of the Hyperparameters
  epsilon = as.vector(apply(Y,2,mean))
  W = diag(diag(cov(Y)))
  
  
  ## Initialization of the parameters for Gaussian Mixture
  mu = matrix(data = NA, nrow = K, ncol = D)
  S = array(data = NA, dim =c(K,D,D))
  
  
  #Sparsity controlling hyperparameter of the BAYESIAN LASSO MODEL
  r =1
  si = 1.78
  
  
  ##Actual parameters
  lambda2 <- numeric(K)
  tau2 = matrix(data = NA, nrow = K, ncol = D)
  betahat = matrix(data = NA, nrow = K, ncol = D)
  sigma2 <- rep(NA, K)
  beta0 <- rep(NA, K)
  That <-  numeric(N)
  
  ## Fitting a linear model to the whole model
  Ysc <- scale(Y[1:N,1:D], center = TRUE, scale =TRUE)
  lm.data <- lm(time ~ Ysc)
  sig2.dat <-  var(lm.data$residuals)
  
  
  ## Set Some Initial Values for the Cluster Parameters
  
  source('priorPARAMETERS.R')
  disclass <- table(factor(c, levels = 1:K))
  activeclass <- which(disclass!=0)
  for ( j in 1:length(activeclass)){
    
    priorone <- priordraw(beta, W, epsilon, ro, r, si, N, D, sig2.dat)  
    mu[activeclass[j],] <- (priorone$mu) 
    S[activeclass[j],1:D,1:D]  <- priorone$Sigma  
    beta0[activeclass[j]] <- priorone$beta0 
    sigma2[activeclass[j]] <- priorone$sigma2
    betahat[activeclass[j],1:D] <- priorone$betahat 
    lambda2[activeclass[j]] <- priorone$lambda2 
    tau2[activeclass[j], 1:D] <- priorone$tau2
  }
  
  # The Time has to be initialized
  source('posteriorCensoredTime.R')
  ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
  That <- ti$time
  
  
  ################# K-Means BLASSO INITIALIZATION ############################################
  G <- F
  
  k.data <- kmeans(Y,F,nstart =10)
  
  c <- k.data$cluster
  
  c.kmeans <- c
  #### Under special cases
  ###### c <- c.true
  
  
  prior.numclust <- table(factor(c, levels = 1:K))
  prior.activeclass<- which(prior.numclust!=0)
  
  ### The means are set using the k-means
  for ( i in 1:length(prior.activeclass)){
    mu[prior.activeclass[i],1:D] <-  k.data$centers[i,1:D] 
    
    S[prior.activeclass[i],1:D,1:D] <-  priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)$Sigma
    
    lclust <- which(c == prior.activeclass[i])
    
    reg.blas <- 0
    
    sum <- c(0)
    
    coeff <- 0
    
    Ytemp <-  matrix(NA, nrow = length(lclust), ncol = D)
    
    Ytemp <- scale(Y[lclust,1:D], center = TRUE, scale = TRUE)
    
    
    ### Part where I use the MONOMVN PACKAGE
    
    Ttemp <- as.vector(That[lclust])
    
    ntemp <- length(lclust)
    
    reg.blas <- blasso(Ytemp, Ttemp, T = 300,thin =  50, RJ = TRUE, mprior = 0.0 ,normalize = TRUE, verb = 0)
    
    sum <- summary(reg.blas, burnin= 100)
    
    ## Selecting those features which are relevant
    
    coeff <- unlist(lapply(strsplit(sum$coef[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
    
    beta0[prior.activeclass[1]] <- coeff[1]
    
    indexplusone <- D+1
    
    ind <- 2:indexplusone
    
    betahat[prior.activeclass[i], ] <- coeff[ind]
    
    ta <- unlist(lapply(strsplit(sum$tau2i[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
    
    tau2[prior.activeclass[i],] <- ta
    
    sigma2[prior.activeclass[i]] <- sum$s2[3]
    
    lambda2[prior.activeclass[i]] <- sum$lambda2[3]
    
  }
  
  ## Deleting those values which are no longer relevant
  g <- table(factor(c, levels = 1:K))
  inactive <- which(g==0)
  
  for ( i in 1:length(inactive)){
    mu[inactive[i],1:D]  <- NA 
    S[inactive[i],1:D,1:D]  <- NA  
    beta0[inactive[i]] <- NA 
    sigma2[inactive[i]] <- NA
    betahat[inactive[i],1:D] <- NA 
    lambda2[inactive[i]] <- NA
    tau2[inactive[i], 1:D] <- NA
  }
  
  
  
  return(list(Time = 'Time', 'alpha' = alpha,'That' = That, 'c'= c,'c.kmeans' = c.kmeans, 'epsilon'= epsilon,'W' = W,'mu'= mu,'S'= S,'beta0'= beta0,'betahat'= betahat,'sigma2'= sigma2, 'lambda2'= lambda2,'tau2'= tau2,'sig2.dat'=sig2.dat))
  
}
