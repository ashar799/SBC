posteriorchineseAFT = function(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat) {
  
  
  Ytemp <- matrix(NA, nrow = N, ncol = D)
  ctemp <- c
  
  ## This can't be parallelized !!!!!
  for(l in 1:N)  {
    
    temp <- ctemp[l]
    cminus <- ctemp
    cminus[l] <- NA
    ## The table function helps converting the data point specific indicator variables to class specific indicator variables
    g <- table(factor(cminus, levels = 1:K))
    active <- which(g!=0)
    
    
    
    kminus <- length(active)
    ## Two Auxilary Variables
    ## The name of the auxilary variables are taken to be one and two more than the maximum value in the already active cluster set
    active <- append(active, max(active)+1)
    active <- append(active, max(active)+1)
    
    
    
    ## If the observation was singelton (i.e no other point was associated with it then we assign to kminus +1 parameter)
    if(length(which(cminus==temp))==0 || length(which(cminus==temp))==1 )  
    {
      ## The kminus+1 parameter gets the value of the temporary variable
      ctemp[l] <- active[kminus+1]
      mu[active[kminus+1],1:D] <- mu[temp,1:D]
      S[active[kminus]+1,1:D,1:D] <- S[temp,1:D,1:D]
      beta0[active[kminus+1]] <- beta0[temp]
      betahat[active[kminus+1], 1:D] <- betahat[temp, 1:D]
      sigma2[active[kminus+1]] <- sigma2[temp]
      lambda2[active[kminus+1]] <- lambda2[temp]
      tau2[active[kminus+1], 1:D] <- tau2[temp, 1:D]
      
      ## Also the second auxilary variable should be drawn from the prior distribution  
      
      priorone <- NA
      priorone <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
      mu[active[kminus+2],1:D]  <- priorone$mu  
      S[active[kminus+2],1:D,1:D]  <- priorone$Sigma  
      beta0[active[kminus+2]] <- priorone$beta0 
      sigma2[active[kminus+2]] <- priorone$sigma2
      betahat[active[kminus+2],1:D] <- priorone$betahat 
      lambda2[active[kminus+2]] <- priorone$lambda2 
      tau2[active[kminus+2], 1:D] <- priorone$tau2
      
      ## As we have to deal with centred matrices and if this point is alone in its cluster then
      for ( k in 1:D){
        Ytemp[l,k] <- 0
      }
      
    }else{
      
      
      
      ## We have to deal with centred matrices
      clust <- which(ctemp == temp)
      tempmatrix <- Y[clust,1:D]
      sd.tempmatrix <- apply(tempmatrix, 2, function(x) sd(x))
      mean.tempmatrix <- apply(tempmatrix, 2, mean)
      
      for ( k in 1:D){
       if (sd.tempmatrix[k] == 0){
         sd.tempmatrix[k] = 1
       }
      }
      
      for ( k in 1:D){
        Ytemp[l,k] <- (Y[l,k] - mean.tempmatrix[k])/(sd.tempmatrix[k])
      }
      
      priortwo <- NA
      priortwo <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
      mu[active[kminus+1],1:D]  <- priortwo$mu  
      S[active[kminus+1],1:D,1:D]  <- priortwo$Sigma[1:D,1:D]  
      beta0[active[kminus+1]] <- priortwo$beta0 
      sigma2[active[kminus+1]] <- priortwo$sigma2
      betahat[active[kminus+1],1:D] <- priortwo$betahat 
      lambda2[active[kminus+1]] <- priortwo$lambda2 
      tau2[active[kminus+1], 1:D] <- priortwo$tau2
      
      
      priorthree <- NA
      priorthree <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
      mu[active[kminus+2],1:D]  <- priorthree$mu  
      S[active[kminus+2],1:D,1:D]  <- priorthree$Sigma[1:D,1:D]  
      beta0[active[kminus+2]] <- priorthree$beta0 
      sigma2[active[kminus+2]] <- priorthree$sigma2
      betahat[active[kminus+2],1:D] <- priorthree$betahat 
      lambda2[active[kminus+2]] <- priorthree$lambda2 
      tau2[active[kminus+2], 1:D] <- priorthree$tau2
    }
    
    
    
    #######################################################
    
    
    posterior <- matrix(NA, nrow = length(active), ncol = 1)
    
    
    
    
    ## Calculating the probabalities for drawing the value of c_i from the active classes
    for (j in 1:kminus) {
      res <- try(dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[j],1:D], Q = S[active[j],1:D,1:D]), silent=TRUE)
      if (class(res) == "try-error"){
        posterior[j] <- 0
      } else{
        posterior[j] <- g[active[j]] /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[j],1:D], Q = S[active[j],1:D,1:D]) *  dnorm(x = That[l], mean = beta0[active[j]] + betahat[active[j],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[j]]) )
      }
      
    }
    
    res <- try(dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]), silent=TRUE)
    if (class(res) == "try-error"){
       posterior[kminus+1] <- 0
    } else{
      posterior[kminus+1] <- (0.5 * alpha) /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]) *  dnorm(x = That[l], mean = beta0[active[kminus+1]] + betahat[active[kminus+1],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[kminus+1]]) )
    }
    
    res2 <- try(dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+2],1:D], Q= S[active[kminus+2],1:D,1:D]), silent=TRUE)
    if (class(res2) == "try-error"){
      posterior[kminus+2] <- 0
    } else{
      posterior[kminus+2] <- (0.5 * alpha) /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+2],1:D], Q = S[active[kminus+2],1:D,1:D]) *  dnorm(x = That[l], mean = beta0[active[kminus+2]] + betahat[active[kminus+2],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[kminus+2]]) )
    }
    
    
#     posterior[kminus+1] <- (0.5 * alpha) /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]) *  dnorm(x = That[l], mean = beta0[active[kminus+1]] + betahat[active[kminus+1],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[kminus+1]]) )
#     
#     posterior[kminus+2] <- (0.5 * alpha) /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+2],1:D], Q = S[active[kminus+2],1:D,1:D]) *  dnorm(x = That[l], mean = beta0[active[kminus+2]] + betahat[active[kminus+2],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[kminus+2]]) )
#     
    
    ## Calculating the normalization constant for probabilities
    normalization <- sum(posterior) 
    
    if (normalization < 1e-200 || normalization ==Inf){
      ctemp[l] <- sample(active, 1, prob = rep(1,length(active)), replace = TRUE)
    } else {  
      ctemp[l] <- sample(active, 1, prob= posterior, replace = TRUE)
    }
    
  }
  
  c <- ctemp
  ## Delete those observations that are not associcated with no data point
  g <- table(factor(c, levels = 1:K))
  inactive <- which(g==0)
  
  for ( i in 1:length(inactive)){
    mu[inactive[i],1:D]  <- NA 
    S[inactive[i],1:D,1:D]  <- NA  
    beta0[inactive[i]] <- NA 
    sigma2[inactive[i]] <- NA
    betahat[inactive[i],1:D] <- NA 
    lambda2[inactive[i]] <- NA
    tau2[inactive[i], 1:D] <- NA
  }
  
  
  
  list('indicator' = c,'mean' = mu,'precision' = S,  'beta0' = beta0,'betahat2'= betahat, 'sigma2'=sigma2, 'lambda2'= lambda2, 'tau2' = tau2)
  
}
