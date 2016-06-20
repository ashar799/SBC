posteriorhyperPLUS = function(c, Y, mu, S, epsilon, W, beta, ro ) {
  
  D = ncol(Y)
  numclust <- table(factor(c, levels = 1:K))
  activeclust <- which(numclust!=0)
  nactive <- length(activeclust)
  InvCov <- solve(cov(Y) + diag(D))
  meandata <- apply(Y, 2, mean )
  meandata <-  as.matrix(meandata)
 
  
  # Update the Epsilon paramter
  sum.precision <- matrix(0, nrow = D, ncol =D)
  sum.mean.precision <- matrix(0, nrow = D, ncol =1)
  for ( z in 1:nactive) {
    sum.precision <- sum.precision + ro * S[activeclust[z],1:D, 1:D]
    sum.mean.precision <-  sum.mean.precision + ro* S[activeclust[z],1:D, 1:D] %*% as.matrix(mu[activeclust[z],1:D]) 
  }
  precision.epsilon <- InvCov + sum.precision
  mean.epsilon <- solve(precision.epsilon) %*% ( InvCov %*% meandata + sum.mean.precision) 
  
  epsilon <- mvrnorm(n=1, mu = as.vector(mean.epsilon), Sigma = solve(precision.epsilon)) 
  
  # Update the ro paramter
  
  sum.ro <- 0
  for ( z in 1:nactive) {
    sum.ro <- sum.ro + t(as.matrix(mu[activeclust[z],1:D]- epsilon)) %*% S[activeclust[z],1:D, 1:D] %*% as.matrix(mu[activeclust[z],1:D]- epsilon)
  }
  ro <- rgamma(1, shape = (nactive/2 + 0.25 ), scale = (as.numeric(sum.ro) +0.5)^-1)
  


#### THE NEW STRUCTURE OF W INVOLVES THAT IT BE A DAIGONAL MATRIX WITH alpha_a variables
alpha_a <- c(rep(0,D))

for ( i in 1:D){
  alpha_a[i] <- rgamma(n =1, shape = 0.5* beta*nactive, rate = 0.5*beta* sum(S[1:nactive,i,i]))    
}
W <- diag(alpha_a)




  
  list('epsilon' = epsilon,'W' = W , 'ro' = ro) 
}
