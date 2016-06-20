multiloglikelihood = function(c,Y1,Y2,D1,D2,That,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2 ) {
  
  gmmx1 <- gmmx1
  gmmx2 <- gmmx2
  regy1 <- regy1
  regy2 <- regy2
  
  disclass <- table(factor(c, levels = 1:K))
  activeclass <- which(disclass!=0)
  
  
  Ytemp1 <- Y1
  Ytemp2 <- Y2
  N = nrow(Y1)
  Ytemp1.scaled <- matrix(NA, nrow = N, ncol = D1)
  Ytemp2.scaled <- matrix(NA, nrow = N, ncol = D2)
 
  
  
  for ( i in 1:length(activeclass)) {
    clust <- which(c == activeclass[i])
    
    if (length(clust)==1){
      Ytemp1.scaled[clust,1:D1] <- matrix(0, nrow =1, ncol =D1)
      Ytemp2.scaled[clust,1:D2] <- matrix(0, nrow =1, ncol =D2)
      
    } else {
    Ytemp1.scaled[clust,1:D1] <- scale(Ytemp1[clust,1:D1], center = TRUE, scale = TRUE)
    Ytemp2.scaled[clust,1:D2] <- scale(Ytemp2[clust,1:D2], center = TRUE, scale = TRUE)
  }
  }
  
  
  loglikelihood <-c(0)
  for (j in 1:length(activeclass)) {
    
    loglikelihood[j] <- 0
    
#     clust <- which(c==activeclass[j])
#     if (length(clust)==1){
#       loglikelihood[j] <- 0
#     } else {
#             for ( l in 1:length(clust)) {
#                if (Time[clust[l],2]==1){
#                     loglikelihood[j] <- loglikelihood[j]   +  dMVN(x = as.vector(t(Ytemp1[clust[l],])), mean = gmmx1$mu[c[clust[l]],1:D1],  Q = gmmx1$S[c[clust[l]],1:D1,1:D1], log =TRUE) +  dMVN(x = as.vector(t(Ytemp2[clust[l],])), mean = gmmx2$mu[c[clust[l]],1:D2],  Q = gmmx2$S[c[clust[l]],1:D2,1:D2], log = TRUE)  +  dnorm(x = Time[clust[l],1], mean = regy1$beta0[c[clust[l]]] + regy1$betahat[c[clust[l]],] %*% as.vector(t(Ytemp1.scaled[clust[l],])), sd = sqrt(regy1$sigma2[c[clust[l]]]), log =TRUE) +  dnorm(x = Time[clust[l],1], mean = regy2$beta0[c[clust[l]]] + regy2$betahat[c[clust[l]],] %*% as.vector(t(Ytemp2.scaled[clust[l],])), sd = sqrt(regy2$sigma2[c[clust[l]]]), log = TRUE) 
#                } else{
#                     loglikelihood[j] <- loglikelihood[j]   +  dMVN(x = as.vector(t(Ytemp1[clust[l],])), mean = gmmx1$mu[c[clust[l]],1:D1],  Q = gmmx1$S[c[clust[l]],1:D1,1:D1], log =TRUE) +  dMVN(x = as.vector(t(Ytemp2[clust[l],])), mean = gmmx2$mu[c[clust[l]],1:D2],  Q = gmmx2$S[c[clust[l]],1:D2,1:D2], log =TRUE)  +  log(dtruncnorm(x = Time[clust[l],1], a = Time[clust[l],1], b = Inf , mean = regy1$beta0[c[clust[l]]] + regy1$betahat[c[clust[l]],] %*% as.vector(t(Ytemp1.scaled[clust[l],])), sd = sqrt(regy1$sigma2[c[clust[l]]])))  + log(dtruncnorm(x = Time[clust[l],1], a = Time[clust[l],1], b = Inf, mean = regy2$beta0[c[clust[l]]] + regy2$betahat[c[clust[l]],] %*% as.vector(t(Ytemp2.scaled[clust[l],])), sd = sqrt(regy2$sigma2[c[clust[l]]]))) 
#                 }
#            }
#     }
clust <- which(c == activeclass[j])
luk1 <- c(0)
luk2 <- c(0)
for ( l in 1:length(clust)) {
  if (Time[clust[l],2]==1){
    luk1[l] <- dMVN(x = as.vector(t(Ytemp1[clust[l],])), mean = gmmx1$mu[c[clust[l]],1:D1],  Q = gmmx1$S[c[clust[l]],1:D1,1:D1], log =TRUE) +  dMVN(x = as.vector(t(Ytemp2[clust[l],])), mean = gmmx2$mu[c[clust[l]],1:D2],  Q = gmmx2$S[c[clust[l]],1:D2,1:D2], log = TRUE) 
    luk2[l] <- dnorm(x = That[clust[l]], mean = regy1$beta0[c[clust[l]]] + regy1$betahat[c[clust[l]],] %*% as.vector(t(Ytemp1.scaled[clust[l],])), sd = sqrt(regy1$sigma2[c[clust[l]]]), log =TRUE) +  dnorm(x = That[clust[l]], mean = regy2$beta0[c[clust[l]]] + regy2$betahat[c[clust[l]],] %*% as.vector(t(Ytemp2.scaled[clust[l],])), sd = sqrt(regy2$sigma2[c[clust[l]]]), log = TRUE) 
  } else{
    luk1[l] <- dMVN(x = as.vector(t(Ytemp1[clust[l],])), mean = gmmx1$mu[c[clust[l]],1:D1],  Q = gmmx1$S[c[clust[l]],1:D1,1:D1], log =TRUE) +  dMVN(x = as.vector(t(Ytemp2[clust[l],])), mean = gmmx2$mu[c[clust[l]],1:D2],  Q = gmmx2$S[c[clust[l]],1:D2,1:D2], log =TRUE)  
    luk2[l] <- log(dtruncnorm(x = Time[clust[l],1], a = Time[clust[l],1], b = Inf, mean = regy1$beta0[c[clust[l]]] + regy1$betahat[c[clust[l]],] %*% as.vector(t(Ytemp1.scaled[clust[l],])), sd = sqrt(regy1$sigma2[c[clust[l]]])))  + log(dtruncnorm(x = Time[clust[l],1], a = Time[clust[l],1], b = Inf, mean = regy2$beta0[c[clust[l]]] + regy2$betahat[c[clust[l]],] %*% as.vector(t(Ytemp2.scaled[clust[l],])), sd = sqrt(regy2$sigma2[c[clust[l]]]))) 
  }
}
loglikelihood[j] <- sum(luk1) + sum(luk2)

}
  
  return(sum(loglikelihood))
  
  
}




# clust <- which(c == activeclass[j])
# luk1 <- c(0)
# luk2 <- c(0)
# for ( l in 1:length(clust)) {
# if (Time[clust[l],2]==1){
#   luk1[l] <- dMVN(x = as.vector(t(Ytemp1[clust[l],])), mean = gmmx1$mu[c[clust[l]],1:D1],  Q = gmmx1$S[c[clust[l]],1:D1,1:D1], log =TRUE) +  dMVN(x = as.vector(t(Ytemp2[clust[l],])), mean = gmmx2$mu[c[clust[l]],1:D2],  Q = gmmx2$S[c[clust[l]],1:D2,1:D2], log = TRUE) 
#   luk2[l] <- dnorm(x = That[clust[l]], mean = regy1$beta0[c[clust[l]]] + regy1$betahat[c[clust[l]],] %*% as.vector(t(Ytemp1.scaled[clust[l],])), sd = sqrt(regy1$sigma2[c[clust[l]]]), log =TRUE) +  dnorm(x = That[clust[l]], mean = regy2$beta0[c[clust[l]]] + regy2$betahat[c[clust[l]],] %*% as.vector(t(Ytemp2.scaled[clust[l],])), sd = sqrt(regy2$sigma2[c[clust[l]]]), log = TRUE) 
# } else{
#   luk1[l] <- dMVN(x = as.vector(t(Ytemp1[clust[l],])), mean = gmmx1$mu[c[clust[l]],1:D1],  Q = gmmx1$S[c[clust[l]],1:D1,1:D1], log =TRUE) +  dMVN(x = as.vector(t(Ytemp2[clust[l],])), mean = gmmx2$mu[c[clust[l]],1:D2],  Q = gmmx2$S[c[clust[l]],1:D2,1:D2], log =TRUE)  
#   luk2[l] <- log(dtruncnorm(x = Time[clust[l],1], a = Time[clust[l],1], b = Inf, mean = regy1$beta0[c[clust[l]]] + regy1$betahat[c[clust[l]],] %*% as.vector(t(Ytemp1.scaled[clust[l],])), sd = sqrt(regy1$sigma2[c[clust[l]]])))  + log(dtruncnorm(x = Time[clust[l],1], a = Time[clust[l],1], b = Inf, mean = regy2$beta0[c[clust[l]]] + regy2$betahat[c[clust[l]],] %*% as.vector(t(Ytemp2.scaled[clust[l],])), sd = sqrt(regy2$sigma2[c[clust[l]]]))) 
# }
# }
# sum(luk1) + sum(luk2)
