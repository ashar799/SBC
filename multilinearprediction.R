###### This function gives a linear predictor for the multi-view case
multipredictlinear = function(c, regy1, regy2 ) {
  
  
  
  numclust <- table(factor(c, levels = 1:K))
  activeclass<- which(numclust!=0)
  
  multiprediction <- c(0)
  
  Ytemp1 <- matrix(NA, nrow = N, ncol = D1)
  for ( i in 1:length(activeclass)) {
    
    clust <- which(c == activeclass[i])
    
    Ytemp1[clust,1:D1] <- scale(Y1[clust,1:D1], center = TRUE, scale = TRUE)
  }
  
  Ytemp2 <- matrix(NA, nrow = N, ncol = D2)
  for ( i in 1:length(activeclass)) {
    
    clust <- which(c == activeclass[i])
    
    Ytemp2[clust,1:D2] <- scale(Y2[clust,1:D2], center = TRUE, scale = TRUE)
  }
  
  for ( h in 1:N){  
    multiprediction[h]<- (regy1$sigma2[c[h]]^-1 *(regy1$beta0[c[h]] + regy1$betahat[c[h],1:D1] %*% Ytemp1[h,1:D1]) +  regy2$sigma2[c[h]]^-1 *(regy2$beta0[c[h]] + regy2$betahat[c[h],1:D2] %*% Ytemp2[h,1:D2]) ) / (regy1$sigma2[c[h]]^-1 + regy2$sigma2[c[h]]^-1)
  }
  
  
  list('predictlinear' =  multiprediction) 
  
  
}
