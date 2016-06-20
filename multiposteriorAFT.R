posteriortimeparameterspenalized = function(c,Y, That, lambda2, tau2, sigma2, beta0, betahat, K, epsilon, W,  beta,ro,r, si, sig2.data,N, D ) {
  
  numclust <- table(factor(c, levels = 1:K))
  activeclass<- which(numclust!=0)
  
  
  for (j in 1:length(activeclass)) {
    
    reg.blas <- 0
    
    sum <- c(0)
    
    coeff <- 0
    
    
    ## A Temporary matrix that needs to store the standardized regressors
    clust <- which(c==activeclass[j])
    
    Ytemp <-  matrix(NA, nrow = length(clust), ncol = D)
    
    if (length(clust)==1){
      Ytemp <- matrix(0, nrow =1, ncol =D)
      
    } else {
      Ytemp <- scale(Y[clust,1:D], center = TRUE, scale = TRUE)
    }
    
    ### Part where I use the MONOMVN PACKAGE
   ## If the Cluster has just one member I draw the parameters from the prior 
    if (length(clust) > 1){
      Ttemp <- as.vector(That[clust])
      ntemp <- length(clust)
      ##reg.blas <- blasso(Ytemp, Ttemp, T =1000,thin = 10, RJ = TRUE, beta = as.vector(betahat[activeclass[j],]),lambda2 = lambda2[activeclass[j]],s2 = sigma2[activeclass[j]], mprior = 0.20 ,rd =c(r,si), ab = c(1,1),normalize = TRUE, verb = 0)
      reg.blas <- blasso(Ytemp, Ttemp, mprior = 0.20 ,rd =c(r,si), ab = c(1,1),normalize = TRUE, verb = 0 )
      
      sum <- summary(reg.blas, burnin= 100)
      
      ## Selecting those features which are relevant
      coeff <- unlist(lapply(strsplit(sum$coef[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
      
      
      
      beta0[activeclass[j]] <- coeff[1]
      
      indexplusone <- D+1
      ind <- 2:indexplusone
      betahat[activeclass[j], ] <- coeff[ind]
      
      ta <- unlist(lapply(strsplit(sum$tau2i[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
      tau2[activeclass[j],] <- ta
      
      sigma2[activeclass[j]] <- sum$s2[3]
      lambda2[activeclass[j]] <- sum$lambda2[3]
    } else {
      
      
      tempvector <- as.vector(That[clust])
      tempmean <- mean(tempvector)
      tmpscl <- scale(tempvector, center = TRUE, scale =FALSE)
      tempmatrix <- Ytemp
      tempnumber <- length(tempvector)
      
      
      tempD <- matrix( 0, nrow = D, ncol =D)
      
      if(any(is.na(tau2[activeclass[j],])) == TRUE)
      {
        tau2[activeclass[j],] <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)$tau2
      }
      
      betahat[activeclass[j],] <- as.vector(priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)$betahat)
      
      for ( i in 1:D ) {
        tempD[i,i] <- tau2[activeclass[j],i]
      }
      
      ## For updating the sparsity prior
      lambda2[activeclass[j]] <- rgamma(1, shape = r+D, rate = si + tr(tempD) )
      
      #For updating tau2
      
      for ( h in 1:D)  {
        tau2[activeclass[j], h] <- (rinv.gaussian(1,mu= sqrt(lambda2[activeclass[j]] * sigma2[activeclass[j]]/ (betahat[activeclass[j],h])^2), lambda = lambda2[activeclass[j]]))^-1
      } 
      
      #For updating sigma2
      ## For updating the sigma2 parameter we need temporary matrices
      
      tempprod <- NA
      
      tempscalesigma1 <- as.vector(tmpscl - Ytemp %*% betahat[activeclass[j], ])
      
      tempprod <- tempscalesigma1 %*% tempscalesigma1
      
      tempscalesigma2 <- NA
      
      tempscalesigma2 <- t(betahat[activeclass[j], ] %*% solve(tempD) %*% betahat[activeclass[j], ] )
      
      
      sigma2[activeclass[j]] <- rinvgamma(1, shape = 1+ 0.5 * (tempnumber +D -1), scale = 1 + (0.5* (tempprod + tempscalesigma2 )) )
      ## This is because the error of the model may make it computationally infeasible
      
      
      ## For updating Betahat we need some matrices
      tempD <- matrix( 0, nrow = D, ncol =D)
      for ( i in 1:D ) {
        tempD[i,i] <- tau2[activeclass[j],i]
      }
      
      tempA <-   matrix(NA, nrow = D, ncol = D)
      
      tempA <- t(Ytemp) %*% Ytemp + solve(tempD)
      
      
      betahat[activeclass[j],] <- mvrnorm(1, mu = solve(tempA) %*% t(tempmatrix) %*% tmpscl, Sigma=  sigma2[activeclass[j]] * solve(tempA))
      
      
      beta0[activeclass[j]] <- rnorm(1, mean = tempmean, sd= sqrt(sigma2[activeclass[j]]/tempnumber))
      
      
      
    }
    
  }
  
  
  
  
  
  
  list('beta0' = beta0,'sigma2' = sigma2, 'betahat' = betahat, 'lambda2' = lambda2, 'tau2' =  tau2 )
}
