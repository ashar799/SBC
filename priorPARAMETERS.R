priordraw = function(beta, W, epsilon, ro, r, si,N,D, sig2.dat) {

  mut = matrix(data = NA, nrow = 1, ncol = D)
  St = array(data = NA, dim =c(D,D))

  tau2t = matrix(data = NA, nrow = 1, ncol = D)
  betahatt = matrix(data = NA, nrow = 1, ncol = D)
  sigma2t <- 0
  beta0t <- 0
  lambda2t <- 0
  
  
  res <- try(rWISHART(1, beta, solve((beta*W))), silent=TRUE)
  if (class(res) == "try-error"){
    St <- solve(W)
  } else{
    St <- rWISHART(1, beta, solve((beta*W)))
  }
  
  
  res2 <- try(as.matrix(rMVN(n=1, mean = epsilon, Q = ro*St)$x), silent = TRUE)
  if (class(res2) == "try-error"){
    mut <- epsilon
  } else{
    mut <- as.matrix(rMVN(n=1, mean = epsilon, Q = ro*St)$x)
  }
  
  
  
  
  
  lambda2t <- rgamma(1,shape = r, rate = si)
  for ( i in 1:D)  {
  tau2t[1, i] <- rgamma(1, shape = 1, rate = lambda2t)
  } 
  ## Approximating the Jeffery's prior by a Gamma distribution 
 
  sigma2t <- mean(rinvgamma(100, shape = 1, scale = 1))
  
  
  beta0t <- mean(rnorm(100, 0, sd = sig2.dat))
  ## 
  scaleof <- 0
  
  scaleof <- sqrt(abs(sigma2t/lambda2t))
  
  for ( i in 1 :D) {
    
    betahatt[1, i] <- urlaplace(1, location = 0, scale = scaleof) 
  }
  
  list('mu' = mut, 'Sigma' = St, 'beta0' = beta0t, 'sigma2' = sigma2t, 'betahat' = betahatt, 'lambda2'= lambda2t, 'tau2'= tau2t)
  

}