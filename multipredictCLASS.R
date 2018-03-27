#### this file just takes the new molecular data
#### THIS function predicts the Class of the new Data Points
#### It is Based on the PredictChineseAFT function

multipredictCLASS = function(Y1.test, Y2.test){
  
  
  N.new <<- nrow(Y1.test)
  c.new.list <- list(0)
  ## The number of posterior samples
  Nps <<- as.integer(iter/ iter.thin)
  That.new <- time.new 
  
  print("GOING THROUGH MCMC Samples")
  pb <- txtProgressBar(min = 1, max = Nps , style = 3)
  
  
  gmmx1.tmp <- list(0)
  gmmx2.tmp <- list(0)
  regy1.tmp <- list(0)
  regy2.tmp <- list(0)
  
  Ytemp1 <- Y1.test
  Ytemp2 <- Y2.test
  Ytemp1.scaled <- matrix(NA, nrow = N, ncol = D1)
  Ytemp2.scaled <- matrix(NA, nrow = N, ncol = D2)
  
  modelweights <- c(0)
  
  
  for (count in 1:Nps){
    
    ## Assign the parameters to the posterior sample
    ctemp <- c.list[[count]]
    gmmx1.tmp <- est.gmmx1[[count]]
    gmmx2.tmp <- est.gmmx2[[count]]
    regy1.tmp <- est.regy1[[count]]
    regy2.tmp <- est.regy2[[count]]
    g <- table(factor(ctemp, levels = 1:K))
    activeclass <- which(g!=0)
    ## The table function helps converting the data point specific indicator variables to class specific indicator variables
    kminus <- length(activeclass)
    ## Two Auxilary Variables
    ## The name of the auxilary variables are taken to be one and two more than the maximum value in the already active cluster set
    activeclass <- append(activeclass, max(activeclass)+1)
    activeclass <- append(activeclass, max(activeclass)+1)
    active <- activeclass 
    ### Assigning values to parameters 
    
    priorone1 <- NA
    priorone2 <- NA
    ### Draw the values of two auxilary parameters from Prior Distribution
    source('priorPARAMETERS.R')
    #priorone1 <- priordraw(beta, gmmx1$W, gmmx1$epsilon, ro, r, si,N,D1, sig2.dat)
    repeat {
      priorone1 <- priordraw(gmmx1.tmp$beta, gmmx1.tmp$W, gmmx1.tmp$epsilon, gmmx1.tmp$ro, r, si,N,D1, sig2.dat)
      res <- try(chol(priorone1$Sigma), silent = TRUE)
      if (class(res) != "try-error"){
        break 
      }
    }
    gmmx1.tmp$mu[active[kminus+1],1:D1]  <- priorone1$mu  
    gmmx1.tmp$S[active[kminus+1],1:D1,1:D1]  <- priorone1$Sigma 
    regy1.tmp$beta0[active[kminus+1]] <- priorone1$beta0 
    regy1.tmp$sigma2[active[kminus+1]] <- priorone1$sigma2
    regy1.tmp$betahat[active[kminus+1],1:D1] <- priorone1$betahat 
    regy1.tmp$lambda2[active[kminus+1]] <- priorone1$lambda2 
    regy1.tmp$tau2[active[kminus+1], 1:D1] <- priorone1$tau2
    
    repeat {
      priorone2 <- priordraw(gmmx2.tmp$beta, gmmx2.tmp$W, gmmx2.tmp$epsilon, gmmx2.tmp$ro, r, si,N, D2, sig2.dat)
      res <- try(chol(priorone2$Sigma), silent = TRUE)
      if (class(res) != "try-error"){
        break 
      }
    }  
    
    gmmx2.tmp$mu[active[kminus+1],1:D2]  <- priorone2$mu  
    gmmx2.tmp$S[active[kminus+1],1:D2,1:D2]  <- priorone2$Sigma 
    regy2.tmp$beta0[active[kminus+1]] <- priorone2$beta0 
    regy2.tmp$sigma2[active[kminus+1]] <- priorone2$sigma2
    regy2.tmp$betahat[active[kminus+1],1:D2] <- priorone2$betahat 
    regy2.tmp$lambda2[active[kminus+1]] <- priorone2$lambda2 
    regy2.tmp$tau2[active[kminus+1], 1:D2] <- priorone2$tau2
    
    source('priorPARAMETERS.R')
    #priorone1 <- priordraw(beta, gmmx1$W, gmmx1$epsilon, ro, r, si,N,D1, sig2.dat)
    repeat {
      priorone1 <- priordraw(beta, gmmx1$W, gmmx1$epsilon, gmmx1$ro, r, si,N,D1, sig2.dat)
      res <- try(chol(priorone1$Sigma),silent = TRUE)
      if (class(res) != "try-error"){
        break 
      }
    }
    gmmx1.tmp$mu[active[kminus+2],1:D1]  <- priorone1$mu  
    gmmx1.tmp$S[active[kminus+2],1:D1,1:D1]  <- priorone1$Sigma 
    regy1.tmp$beta0[active[kminus+2]] <- priorone1$beta0 
    regy1.tmp$sigma2[active[kminus+2]] <- priorone1$sigma2
    regy1.tmp$betahat[active[kminus+2],1:D1] <- priorone1$betahat 
    regy1.tmp$lambda2[active[kminus+2]] <- priorone1$lambda2 
    regy1.tmp$tau2[active[kminus+2], 1:D1] <- priorone1$tau2
    
    ##priorone2 <- priordraw(beta, gmmx2$W, gmmx2$epsilon, ro, r, si,N,D2, sig2.dat)
    repeat {
      priorone2 <- priordraw(beta, gmmx2$W, gmmx2$epsilon, gmmx2$ro, r, si,N,D2, sig2.dat)
      res <- try(chol(priorone2$Sigma), silent = TRUE)
      if (class(res) != "try-error"){
        break 
      }
    }  
    gmmx2.tmp$mu[active[kminus+2],1:D2]  <- priorone2$mu  
    gmmx2.tmp$S[active[kminus+2],1:D2,1:D2]  <- priorone2$Sigma 
    regy2.tmp$beta0[active[kminus+2]] <- priorone2$beta0 
    regy2.tmp$sigma2[active[kminus+2]] <- priorone2$sigma2
    regy2.tmp$betahat[active[kminus+2],1:D2] <- priorone2$betahat 
    regy2.tmp$lambda2[active[kminus+2]] <- priorone2$lambda2 
    regy2.tmp$tau2[active[kminus+2], 1:D2] <- priorone2$tau2
    
    #######################################################
    ctemp.new = c(0)
    
    weights.final <- c(0)
    ## This can't be parallelized !!!!!
    for(l in 1:N.new)  {
      
      
      posterior <- matrix(NA, nrow = length(active), ncol = 1)
      Y.new.sc1 <- matrix(0, nrow = N.new, ncol =D1)
      Y.new.sc2 <- matrix(0, nrow = N.new, ncol =D2)
      ## Calculating the probabalities for drawing the value of c_i from the active classes
      for (j in 1:kminus) {
        
        clust <- which(ctemp == active[j])
        
        posterior[j] <- log(g[active[j]] /(N-1+alpha)) + dMVN(x = as.vector(t(Ytemp1[l,])), mean = gmmx1.tmp$mu[active[j],1:D1],  Q = gmmx1.tmp$S[active[j],1:D1,1:D1], log = TRUE) +  dMVN(x = as.vector(t(Ytemp2[l,])), mean = gmmx2.tmp$mu[active[j],1:D2],  Q = gmmx2.tmp$S[active[j],1:D2,1:D2], log =TRUE) 
      }
      
      posterior[kminus+1] <- log((0.5 * alpha) /(N-1+alpha)) + dMVN(x = as.vector(t(Ytemp1[l,])), mean = gmmx1.tmp$mu[active[kminus+1],1:D1],  Q = gmmx1.tmp$S[active[kminus+1],1:D1,1:D1], log = TRUE) +  dMVN(x = as.vector(t(Ytemp2[l,])), mean = gmmx2.tmp$mu[active[kminus+1],1:D2],  Q = gmmx2.tmp$S[active[kminus+1],1:D2,1:D2], log = TRUE) 
      posterior[kminus+2] <- log((0.5 * alpha) /(N-1+alpha)) + dMVN(x = as.vector(t(Ytemp1[l,])), mean = gmmx1.tmp$mu[active[kminus+2],1:D1],  Q = gmmx1.tmp$S[active[kminus+2],1:D1,1:D1], log = TRUE) +  dMVN(x = as.vector(t(Ytemp2[l,])), mean = gmmx2.tmp$mu[active[kminus+2],1:D2],  Q = gmmx2.tmp$S[active[kminus+2],1:D2,1:D2], log = TRUE)  
      
      ## Calculating the normalization constant for probabilities
      post <- exp(posterior) 
      
      if (sum(post) > 0){
        ctemp.new[l] <- sample(active, 1, prob= post, replace = TRUE)
      } else {
        ctemp.new[l] <- sample(active, 1)
      }
    
      weights.final[l] <- dMVN(x = as.vector(t(Ytemp1[l,1:D1])), mean = gmmx1.tmp$mu[ctemp.new[l] ,1:D1],  Q = gmmx1.tmp$S[ctemp.new[l] ,1:D1,1:D1], log = TRUE) +  dMVN(x = as.vector(t(Ytemp2[l,1:D2])), mean = gmmx2.tmp$mu[ctemp.new[l],1:D2],  Q = gmmx2.tmp$S[ctemp.new[l],1:D2,1:D2], log =TRUE) 
      }
    modelweights[count] <- sum(weights.final)
    c.new.list[[count]] <- ctemp.new 
    Sys.sleep(0.1)
    setTxtProgressBar(pb, count)
    
  }
  
  c.matrix.new <- matrix(NA, nrow = N.new, ncol = Nps)
  for( h in 1:Nps){
    c.matrix.new[,h] <- c.new.list[[h]]
  }
  
  
  c.matrix.new <<- c.matrix.new
  
  ### Build A consensus clustering based on the posterior matrix for both training and testing labels
  psm2 <- comp.psm(t(c.matrix.new))
  mpear2 <- maxpear(psm2)
  adjustedRandIndex(c.true.new,mpear2$cl)
  
  
  
  ### Generally the MPEAR output needs post-processing
  ### If we build a cluster specific sbc approach
  
  c.sbc.new <<- mpear2$cl
  
  #### To calculate the posterior probabilities
  
  test.modelweights <<- modelweights
  
}
