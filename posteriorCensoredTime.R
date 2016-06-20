updatetime = function(c, Y, Time,That, beta0, betahat, sigma2 ) {
  
  Ytemp <- matrix(NA, nrow = N, ncol = D)
  numclust <- table(factor(c, levels = 1:K))
  activeclass<- which(numclust!=0)
  
  
  ## Updating the Estimated Survival Time
  for ( i in 1:length(activeclass)) {
    
    clust <- which(c == activeclass[i])
    
    Ytemp[clust,1:D] <- scale(Y[clust,1:D], center = TRUE, scale = TRUE)
  }
  
      
      
  for ( h in 1:N){     
     if(Time[h,2] == 0){
     That[h]<- rtruncnorm(1, a = Time[h,1], b = Inf, mean = beta0[c[h]] + betahat[c[h],1:D ] %*% Ytemp[h,1:D] , sd = sqrt(sigma2[c[h]]) )
    }else {
      That[h] <- Time[h,1] }
  }
  
  list('time' = That) 
  
  
  
}
