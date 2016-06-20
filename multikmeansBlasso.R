multikmeansBlasso = function(c,Y1,Y2,D1,D2,That,K, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2,surv.obj ) {
  
  gmmx1 <- gmmx1
  gmmx2 <- gmmx2
  
  regy1 <- regy1
  regy2 <- regy2
  
  
#   pc1 <- prcomp(Y1)
#   pc.pred1 <- predict(pc1,newdata = Y1)[,1]
#   pc2 <- prcomp(Y2)
#   pc.pred2 <- predict(pc2, newdata = Y2)[,2]
#   
  
  Yg <- cbind(Y1,Y2)
  pval =c(0)
  pval[1] =1
  
  G <- F
  k.data <- kmeans(Yg,F,nstart =10)
  c <- k.data$cluster

            
            
  Ygr <- cbind(Y1,Y2)          
  Dg <- D1 +D2
  mug = matrix(data = NA, nrow = K, ncol = Dg)
  betahatg = matrix(data = NA, nrow = K, ncol = Dg)
  tau2g = matrix(data = NA, nrow = K, ncol = Dg)
  sigma2 <- rep(NA, K)
  lambda2g <- numeric(K)
  beta0 <- rep(NA, K)
  
  
  
  source('priorPARAMETERS.R')
  prior.numclust <- table(factor(c, levels = 1:K))
  prior.activeclass<- which(prior.numclust!=0)
  
  ### The means are set to the cluster means
  
  for ( i in 1:length(prior.activeclass)){
    
    lclust <- which(c == prior.activeclass[i])
    
    
    mug[prior.activeclass[i],1:Dg] <-  apply(Ygr[lclust,],2,mean)
    
    gmmx1$S[prior.activeclass[i],1:D1,1:D1] <-  priordraw(gmmx1$beta, gmmx1$W, gmmx1$epsilon, gmmx1$ro, r, si,N,D1, sig2.dat)$Sigma
    
    gmmx2$S[prior.activeclass[i],1:D2,1:D2] <-  priordraw(gmmx2$beta, gmmx2$W, gmmx2$epsilon, gmmx2$ro, r, si,N,D2, sig2.dat)$Sigma
    
    lclust <- which(c == prior.activeclass[i])
    
    reg.blas <- 0
    
    sum <- c(0)
    
    coeff <- 0
    
    Ytemp <-  matrix(NA, nrow = length(lclust), ncol = Dg)
    
    Ytemp <- scale(Ygr[lclust,1:Dg], center = TRUE, scale = TRUE)
    
    ### Part where I use the MONOMVN PACKAGE
    
    Ttemp <- as.vector(That[lclust])
    
    ntemp <- length(lclust)
      
    reg.blas <- blasso(Ytemp, Ttemp, T = 200,thin =  10, RJ = TRUE, mprior = 0.0 ,normalize = TRUE, verb = 0)
      
    sum <- summary(reg.blas, burnin= 50)
      
    ## Selecting those features which are relevant
    
    coeff <- unlist(lapply(strsplit(sum$coef[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
      
    regy1$beta0[prior.activeclass[i]] <- coeff[1]
    
    regy2$beta0[prior.activeclass[i]] <- coeff[1]
    
      
    indexplusone <- Dg+1
      
    ind <- 2:indexplusone
      
    betahatg[prior.activeclass[i], ] <- coeff[ind]
      
    ta <- unlist(lapply(strsplit(sum$tau2i[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
      
    tau2g[prior.activeclass[i],] <- ta
      
    sigma2[prior.activeclass[i]] <- sum$s2[3]
      
    lambda2g[prior.activeclass[i]] <- sum$lambda2[3]
    
  }
  
  
  
  ## Deleting those values which are no longer relevant
  g <- table(factor(c, levels = 1:K))
  inactive <- which(g==0)
  
  for ( i in 1:length(inactive)){
    mug[inactive[i],1:Dg]  <- NA 
    gmmx1$S[inactive[i],1:D1,1:D1]  <- NA 
    gmmx2$S[inactive[i],1:D2,1:D2]  <- NA 
    regy1$beta0[inactive[i]] <- NA
    regy2$beta0[inactive[i]] <- NA
    regy2$sigma2[inactive[i]] <- NA
    regy1$sigma2[inactive[i]] <- NA
    betahatg[inactive[i],1:Dg] <- NA 
    lambda2g[inactive[i]] <- NA
    sigma2[inactive[i]] <- NA
    tau2g[inactive[i], 1:Dg] <- NA
  }
  
  indte <- D1+1
  
  gmmx1$mu <-  mug[,1:D1] 
  
  gmmx2$mu <-  mug[,indte:Dg]
  
  regy1$betahat <-  betahatg[,1:D1]
  
  regy2$betahat <- betahatg[,indte:Dg]
  
  regy1$tau2 <- tau2g[,1:D1]
  
  regy2$tau2 <-  tau2g[,indte:Dg]
  
  regy1$lambda2 <- lambda2g
  
  regy2$lambda2 <- lambda2g

  regy1$sigma2 <-  sigma2
 
  regy2$sigma2 <- sigma2
  
  
  ########################## THE HYPERPARAMETERS OF THE GMM are initialized to Empirical Bayes #################################  
   
  
  
  gmmx1$epsilon <- as.vector(apply(Y1,2,mean))
  gmmx1$W <- diag(diag(as.matrix(cov(Y1))))
  
  ##Updating the hyper parameter for the second data set
  
  gmmx2$epsilon <- as.vector(apply(Y2,2,mean))
  gmmx2$W <- diag(diag(as.matrix(cov(Y2))))
  
  
  list('c'=c,'gmmx1'=gmmx1,'gmmx2'= gmmx2, 'regy1'= regy1,'regy2'= regy2)  
  
}