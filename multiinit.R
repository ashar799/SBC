multiinit = function(Y,c,beta, W, epsilon, ro, r, si,N,D, sig2.dat) {
  
  mu = matrix(data = NA, nrow = K, ncol = D)
  S = array(data = NA, dim =c(K,D,D))
  tau2 = matrix(data = NA, nrow = K, ncol = D)
  betahat = matrix(data = NA, nrow = K, ncol = D)
  lambda2 <- numeric(K)
  sigma2 <- rep(NA, K)
  beta0 <- rep(NA, K)
  
  
  source('priorPARAMETERS.R')
  disclass <- table(factor(c, levels = 1:K))
  activeclass <- which(disclass!=0)
  for ( j in 1:length(activeclass)){
    priorone <- priordraw(beta, W, epsilon, ro, r, si, N, D, sig2.dat)  
    mu[activeclass[j],] <- priorone$mu
    S[activeclass[j],1:D,1:D]  <- priorone$Sigma 
    beta0[activeclass[j]] <- priorone$beta0 
    sigma2[activeclass[j]] <- priorone$sigma2
    betahat[activeclass[j],1:D] <- priorone$betahat
    lambda2[activeclass[j]] <- priorone$lambda2 
    tau2[activeclass[j], 1:D] <- priorone$tau2
   }
  
  
  list('mu'=mu, 'beta0'=beta0, 'betahat'= betahat, 'sigma2' =sigma2, 'lambda2' = lambda2, 'tau2'= tau2,  'S' =S)  
}
