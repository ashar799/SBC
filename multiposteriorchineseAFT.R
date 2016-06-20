multiposteriorchineseAFT = function(c,Y1,Y2,D1,D2,That, K, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2) {
  
  numclust <- table(factor(c, levels = 1:K))
  activeclass<- which(numclust!=0)
  
  
  
  
  gmmx1 <- gmmx1
  gmmx2 <- gmmx2
  
  regy1 <- regy1
  regy2 <- regy2

  
  ctemp <- c
  Ytemp1 <- Y1
  Ytemp2 <- Y2
  
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
  
  
  
  
  ## This can't be parallelized !!!!!
  for(l in 1:N)  {
    
    temp <- ctemp[l]
    cminus <- ctemp
    cminus[l] <- NA
    ## The table function helps converting the data point specific indicator variables to class specific indicator variables
    g <- table(factor(cminus, levels = 1:K))
    active <- which(g!=0)
    
    
    ### Length of the number of active clusters
    kminus <- length(active)
    ## Two Auxilary Variables
    ## The name of the auxilary variables are taken to be one and two more than the maximum value in the already active cluster set
    active <- append(active, max(active)+1)
    active <- append(active, max(active)+1)
    
    
    
    ## If the observation was singelton (i.e no other point was associated with it then we assign to kminus +1 parameter)
    if(length(which(cminus==temp))==0)  
    {
      ## The kminus+1 parameter gets the value of the temporary variable
      ctemp[l] <- active[kminus+1]
      gmmx1$mu[active[kminus+1],1:D1] <- gmmx1$mu[temp,1:D1]
      gmmx1$S[active[kminus]+1,1:D1,1:D1] <- gmmx1$S[temp,1:D1,1:D1]
      regy1$beta0[active[kminus+1]] <- regy1$beta0[temp]
      regy1$betahat[active[kminus+1], 1:D1] <- regy1$betahat[temp, 1:D1]
      regy1$sigma2[active[kminus+1]] <- regy1$sigma2[temp]
      regy1$lambda2[active[kminus+1]] <- regy1$lambda2[temp]
      regy1$tau2[active[kminus+1], 1:D1] <- regy1$tau2[temp, 1:D1]
      
      
      gmmx2$mu[active[kminus+1],1:D2] <- gmmx2$mu[temp,1:D2]
      gmmx2$S2[active[kminus]+1,1:D2,1:D2] <- gmmx2$S2[temp,1:D2,1:D2]
      regy2$beta0[active[kminus+1]] <- regy2$beta0[temp]
      regy2$betahat[active[kminus+1], 1:D2] <- regy2$betahat[temp, 1:D2]
      regy2$sigma2[active[kminus+1]] <- regy2$sigma2[temp]
      regy2$lambda2[active[kminus+1]] <- regy2$lambda2[temp]
      regy2$tau2[active[kminus+1], 1:D2] <- regy2$tau2[temp, 1:D2]
    
      
      
      ## Also the second auxilary variable should be drawn from the prior distribution  
      
      source('priordraw.R')
      priorone1 <- NA
      priorone2 <- NA
  ##### priorone1 <- priordraw(beta, gmmx1$W, gmmx1$epsilon, ro, r, si,N,D1, sig2.dat)
      
      repeat
      {
        priorone1 <- priordraw(beta, gmmx1$W, gmmx1$epsilon, gmmx1$ro, r, si,N,D1, sig2.dat)
        res <- try(chol(priorone1$Sigma), silent = TRUE)
        if (class(res) != "try-error"){
        break 
        }
        }
      gmmx1$mu[active[kminus+2],1:D1]  <- priorone1$mu  
      gmmx1$S[active[kminus+2],1:D1,1:D1]  <- priorone1$Sigma 
      regy1$beta0[active[kminus+2]] <- priorone1$beta0 
      regy1$sigma2[active[kminus+2]] <- priorone1$sigma2
      regy1$betahat[active[kminus+2],1:D1] <- priorone1$betahat 
      regy1$lambda2[active[kminus+2]] <- priorone1$lambda2 
      regy1$tau2[active[kminus+2], 1:D1] <- priorone1$tau2
#       priorone2 <- priordraw(beta, gmmx2$W, gmmx2$epsilon, ro, r, si,N,D2, sig2.dat)
 
     repeat
       {
         priorone2 <- priordraw(beta, gmmx2$W, gmmx2$epsilon, gmmx2$ro, r, si,N,D2, sig2.dat)
         res <- try(chol(priorone2$Sigma), silent = TRUE)
         if (class(res) != "try-error"){
         break 
        }
      }  
 
      gmmx2$mu[active[kminus+2],1:D2]  <- priorone2$mu  
      gmmx2$S[active[kminus+2],1:D2,1:D2]  <- priorone2$Sigma 
      regy2$beta0[active[kminus+2]] <- priorone2$beta0 
      regy2$sigma2[active[kminus+2]] <- priorone2$sigma2
      regy2$betahat[active[kminus+2],1:D2] <- priorone2$betahat 
      regy2$lambda2[active[kminus+2]] <- priorone2$lambda2 
      regy2$tau2[active[kminus+2], 1:D2] <- priorone2$tau2
      
      
      ## The number of points in the cluster
      
  
    } else {
      
      priorone1 <- NA
      priorone2 <- NA
      ### Draw the values of two auxilary parameters from Prior Distribution
      source('priordraw.R')
      #priorone1 <- priordraw(beta, gmmx1$W, gmmx1$epsilon, ro, r, si,N,D1, sig2.dat)
      repeat {
      priorone1 <- priordraw(beta, gmmx1$W, gmmx1$epsilon, gmmx1$ro, r, si,N,D1, sig2.dat)
      res <- try(chol(priorone1$Sigma), silent = TRUE)
      if (class(res) != "try-error"){
      break 
        }
      }
      gmmx1$mu[active[kminus+1],1:D1]  <- priorone1$mu  
      gmmx1$S[active[kminus+1],1:D1,1:D1]  <- priorone1$Sigma 
      regy1$beta0[active[kminus+1]] <- priorone1$beta0 
      regy1$sigma2[active[kminus+1]] <- priorone1$sigma2
      regy1$betahat[active[kminus+1],1:D1] <- priorone1$betahat 
      regy1$lambda2[active[kminus+1]] <- priorone1$lambda2 
      regy1$tau2[active[kminus+1], 1:D1] <- priorone1$tau2
      # priorone2 <- priordraw(beta, gmmx2$W, gmmx2$epsilon, ro, r, si,N,D2, sig2.dat)
       repeat {
      priorone2 <- priordraw(beta, gmmx2$W, gmmx2$epsilon, gmmx2$ro, r, si,N, D2, sig2.dat)
      res <- try(chol(priorone2$Sigma), silent = TRUE)
      if (class(res) != "try-error"){
      break 
       }
      }  

      gmmx2$mu[active[kminus+1],1:D2]  <- priorone2$mu  
      gmmx2$S[active[kminus+1],1:D2,1:D2]  <- priorone2$Sigma 
      regy2$beta0[active[kminus+1]] <- priorone2$beta0 
      regy2$sigma2[active[kminus+1]] <- priorone2$sigma2
      regy2$betahat[active[kminus+1],1:D2] <- priorone2$betahat 
      regy2$lambda2[active[kminus+1]] <- priorone2$lambda2 
      regy2$tau2[active[kminus+1], 1:D2] <- priorone2$tau2
      
      source('priordraw.R')
      #priorone1 <- priordraw(beta, gmmx1$W, gmmx1$epsilon, ro, r, si,N,D1, sig2.dat)
      repeat {
        priorone1 <- priordraw(beta, gmmx1$W, gmmx1$epsilon, gmmx1$ro, r, si,N,D1, sig2.dat)
        res <- try(chol(priorone1$Sigma),silent = TRUE)
        if (class(res) != "try-error"){
          break 
        }
      }
      gmmx1$mu[active[kminus+2],1:D1]  <- priorone1$mu  
      gmmx1$S[active[kminus+2],1:D1,1:D1]  <- priorone1$Sigma 
      regy1$beta0[active[kminus+2]] <- priorone1$beta0 
      regy1$sigma2[active[kminus+2]] <- priorone1$sigma2
      regy1$betahat[active[kminus+2],1:D1] <- priorone1$betahat 
      regy1$lambda2[active[kminus+2]] <- priorone1$lambda2 
      regy1$tau2[active[kminus+2], 1:D1] <- priorone1$tau2

      ##priorone2 <- priordraw(beta, gmmx2$W, gmmx2$epsilon, ro, r, si,N,D2, sig2.dat)
      repeat {
        priorone2 <- priordraw(beta, gmmx2$W, gmmx2$epsilon, gmmx2$ro, r, si,N,D2, sig2.dat)
        res <- try(chol(priorone2$Sigma), silent = TRUE)
        if (class(res) != "try-error"){
          break 
        }
      }  
      gmmx2$mu[active[kminus+2],1:D2]  <- priorone2$mu  
      gmmx2$S[active[kminus+2],1:D2,1:D2]  <- priorone2$Sigma 
      regy2$beta0[active[kminus+2]] <- priorone2$beta0 
      regy2$sigma2[active[kminus+2]] <- priorone2$sigma2
      regy2$betahat[active[kminus+2],1:D2] <- priorone2$betahat 
      regy2$lambda2[active[kminus+2]] <- priorone2$lambda2 
      regy2$tau2[active[kminus+2], 1:D2] <- priorone2$tau2
      
   
    }
    
    
    
    #######################################################
    
    
    posterior <- matrix(NA, nrow = length(active), ncol = 1)
    
  ## Calculating the probabalities for drawing the value of c_i from the active classes
      for (j in 1:kminus) {
       
      posterior[j] <- log(g[active[j]] /(N-1+alpha)) + dMVN(x = as.vector(t(Ytemp1[l,])), mean = gmmx1$mu[active[j],1:D1],  Q = gmmx1$S[active[j],1:D1,1:D1], log = TRUE) +  dMVN(x = as.vector(t(Ytemp2[l,])), mean = gmmx2$mu[active[j],1:D2],  Q = gmmx2$S[active[j],1:D2,1:D2], log =TRUE) +    dnorm(x = That[l], mean = regy1$beta0[active[j]] + regy1$betahat[active[j],] %*% as.vector(t(Ytemp1.scaled[l,])), sd = sqrt(regy1$sigma2[active[j]]),log = TRUE ) +  dnorm(x = That[l], mean = regy2$beta0[active[j]] + regy2$betahat[active[j],] %*% as.vector(t(Ytemp2.scaled[l,])), sd = sqrt(regy2$sigma2[active[j]]), log = TRUE )  
      }
    
    
      posterior[kminus+1] <- log((0.5 * alpha) /(N-1+alpha)) + dMVN(x = as.vector(t(Ytemp1[l,])), mean = gmmx1$mu[active[kminus+1],1:D1],  Q = gmmx1$S[active[kminus+1],1:D1,1:D1], log = TRUE) +  dMVN(x = as.vector(t(Ytemp2[l,])), mean = gmmx2$mu[active[kminus+1],1:D2],  Q = gmmx2$S[active[kminus+1],1:D2,1:D2], log = TRUE)  +  dnorm(x = That[l], mean = regy1$beta0[active[kminus+1]] + regy1$betahat[active[kminus+1],] %*% as.vector(t(Ytemp1.scaled[l,])), sd = sqrt(regy1$sigma2[active[kminus+1]]), log =TRUE ) +  dnorm(x = That[l], mean = regy2$beta0[active[kminus+1]] + regy2$betahat[active[kminus+1],] %*% as.vector(t(Ytemp2.scaled[l,])), sd = sqrt(regy2$sigma2[active[kminus+1]]), log =TRUE)   
      posterior[kminus+2] <- log((0.5 * alpha) /(N-1+alpha)) + dMVN(x = as.vector(t(Ytemp1[l,])), mean = gmmx1$mu[active[kminus+2],1:D1],  Q = gmmx1$S[active[kminus+2],1:D1,1:D1], log = TRUE) +  dMVN(x = as.vector(t(Ytemp2[l,])), mean = gmmx2$mu[active[kminus+2],1:D2],  Q = gmmx2$S[active[kminus+2],1:D2,1:D2], log = TRUE)  +  dnorm(x = That[l], mean = regy1$beta0[active[kminus+2]] + regy1$betahat[active[kminus+2],] %*% as.vector(t(Ytemp1.scaled[l,])), sd = sqrt(regy1$sigma2[active[kminus+2]]), log =TRUE ) +  dnorm(x = That[l], mean = regy2$beta0[active[kminus+2]] + regy2$betahat[active[kminus+2],] %*% as.vector(t(Ytemp2.scaled[l,])), sd = sqrt(regy2$sigma2[active[kminus+2]]), log =TRUE)   
    
    
    
#     posterior[kminus+1] <- (0.5 * alpha) /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]) *  dnorm(x = That[l], mean = beta0[active[kminus+1]] + betahat[active[kminus+1],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[kminus+1]]) )
#     
#     posterior[kminus+2] <- (0.5 * alpha) /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+2],1:D], Q = S[active[kminus+2],1:D,1:D]) *  dnorm(x = That[l], mean = beta0[active[kminus+2]] + betahat[active[kminus+2],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[kminus+2]]) )
#     
    
    ## Calculating the normalization constant for probabilities
    post <- exp(posterior) 
    
    ctemp[l] <- sample(active, 1, prob= post, replace = TRUE)
    }
    
 
  
  c <- ctemp
  


  ######## Delete those observations that are not associcated with no data point #################3
  g <- table(factor(c, levels = 1:K))
  inactive <- which(g==0)
  for ( i in 1:length(inactive)){
    gmmx1$mu[inactive[i],]  <- NA 
    gmmx1$S[inactive[i],1:D1,1:D1]  <- NA 
    gmmx2$mu[inactive[i],]  <- NA 
    gmmx2$S[inactive[i],1:D2,1:D2]  <- NA 
    regy1$beta0[inactive[i]] <- NA
    regy2$beta0[inactive[i]] <- NA
    regy2$sigma2[inactive[i]] <- NA
    regy1$sigma2[inactive[i]] <- NA
    regy1$betahat[inactive[i],] <- NA
    regy2$betahat[inactive[i],] <- NA
    regy1$lambda2[inactive[i]] <- NA
    regy2$lambda2[inactive[i]] <- NA
    regy1$tau2[inactive[i],] <- NA
    regy2$tau2[inactive[i],] <- NA
    
}

  
  
  
  list('c'=c,'gmmx1'=gmmx1,'gmmx2'= gmmx2,'regy1'= regy1,'regy2'= regy2)
  
}
