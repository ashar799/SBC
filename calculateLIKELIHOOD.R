loglikelihood = function(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat) {
  
  disclass <- table(factor(c, levels = 1:K))
  activeclass <- which(disclass!=0)
  
  
  Ytemp <- Y
  Ytemp.scaled <- matrix(NA, nrow = N, ncol = D)
  
  
  
  for ( i in 1:length(activeclass)) {
    clust <- which(c == activeclass[i])
    
    if (length(clust)==1){
      Ytemp.scaled[clust,1:D] <- matrix(0, nrow =1, ncol =D)
      } else {
      Ytemp.scaled[clust,1:D] <- scale(Ytemp[clust,1:D], center = TRUE, scale = TRUE)
      }
  }
  
  loglikelihood <-c(0)
  GMMlikelihood <- c(0)
  AFTlikelihood <- c(0)
  
  for (j in 1:length(activeclass)) {
    loglikelihood[j] <- 0
    clust <- which(c == activeclass[j])
    luk1 <- c(0)
    luk2 <- c(0)
    for ( l in 1:length(clust)) {
        luk1[l] <- dMVN(x = as.vector(t(Y[clust[l],1:D])), mean = mu[activeclass[j],1:D],  Q = S[activeclass[j],1:D,1:D], log =TRUE)
        luk2[l] <-  dnorm(x = That[clust[l]], mean = beta0[activeclass[j]] + betahat[activeclass[j],1:D] %*% as.vector(t(Ytemp.scaled[clust[l],1:D])), sd = sqrt(sigma2[activeclass[j]]) , log =TRUE)
      } 
     
    GMMlikelihood[j] <- sum(luk1)
    AFTlikelihood[j] <- sum(luk2)
    loglikelihood[j] <-  GMMlikelihood[j] +  AFTlikelihood[j]
    
    
    }
    
 list('loglikelihood' = sum(loglikelihood),'GMMlikelihood' = sum(GMMlikelihood), 'AFTlikelihood' = sum(AFTlikelihood))
  
  
}