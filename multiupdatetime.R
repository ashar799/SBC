multiupdatetime = function(c, Y1, Y2, Time,That, regy1, regy2) {
  
  Ytemp1 <- matrix(NA, nrow = N, ncol = D1)
  Ytemp2 <- matrix(NA, nrow = N, ncol = D2)
  
  numclust <- table(factor(c, levels = 1:K))
  activeclass<- which(numclust!=0)
  
  
  ## Updating the Estimated Survival Time
  for ( i in 1:length(activeclass)) {
    
    clust <- which(c == activeclass[i])
    
    Ytemp1[clust,1:D1] <- scale(Y1[clust,1:D1], center = TRUE, scale = TRUE)
    
    Ytemp2[clust,1:D2] <- scale(Y2[clust,1:D2], center = TRUE, scale = TRUE)
      
  }
  
  for ( h in 1:N){     
     if(Time[h,2] == 0){
     r1 <- rtruncnorm(1, a = Time[h,1], b = Inf, mean = regy1$beta0[c[h]] + regy1$betahat[c[h],1:D1] %*% Ytemp1[h,1:D1] , sd = sqrt(regy1$sigma2[c[h]]) )
     r2 <- rtruncnorm(1, a = Time[h,1], b = Inf, mean = regy2$beta0[c[h]] + regy2$betahat[c[h],1:D2] %*% Ytemp1[h,1:D2] , sd = sqrt(regy2$sigma2[c[h]]) )
     That[h] <-  (regy1$sigma2[c[h]] * r1 +  regy2$sigma2[c[h]] * r2) /(r1 + r2)  
     }else {
      That[h] <- Time[h,1] }
  }
  
  list('time' = That) 
  
  
  
}
