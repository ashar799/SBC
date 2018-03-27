### This is the multi view extension of the single data source case
### This function takes the posterior parameters AND predicts the time for the new points
### The fundamental assumption is that EACH NEW TEST POINT IS CONDITIONALLY INDEPENDENT on the OTHER POINTS
### We predict value of one point at a time
### The final output is Time for the new samples, ONE AT A TIME


multipredictchineseAFTtime = function(Y1.test, Y2.test){
  
  
  c.new.list <- list(0)
  ## The number of posterior samples
  
  post.time  = matrix(NA,nrow = nrow(Y1.test), ncol = Nps)
  cind <- c(0) 
  
  
  N.new <- nrow(Y1.test)
  
  gmmx1.tmp <- list(0)
  gmmx2.tmp <- list(0)
  regy1.tmp <- list(0)
  regy2.tmp <- list(0)
  
  Ytemp1 <- Y1.test
  Ytemp2 <- Y2.test
  
  
  
  print("GOING THROUGH MCMC Samples")
  pb <- txtProgressBar(min = 1, max = Nps , style = 3)
  
  for (count in 1:Nps){
    
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
      priorone1 <- priordraw(gmmx1$beta, gmmx1$W, gmmx1$epsilon, gmmx1$ro, r, si,N,D1, sig2.dat)
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
      priorone2 <- priordraw(gmmx2$beta, gmmx2$W, gmmx2$epsilon, gmmx2$ro, r, si,N,D2, sig2.dat)
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
    
    #################################################
    Y1.new.scaled.list <- list(0)
    Y2.new.scaled.list <- list(0)
    
    ###### Some quantities used to store probabilities  
    posteriortime <- matrix(0, nrow = length(active), ncol = N.new)
    posteriortimeweight <- matrix(0, nrow = length(active), ncol = N.new)
    weights <- matrix(0, nrow = length(active), ncol = N.new)
    
    
    ####### This can't be parallelized !!!!! #####################################
    
    for(l in 1:N.new)  {
      
      ## Calculating the Expectations and also the normalization constant for the Expectation
      for (j in 1:kminus) {
        
        clust <- which(ctemp == active[j])
        obj.t1 <- scale(Y1[clust,1:D1], center = TRUE, scale = TRUE)
        obj.t2 <- scale(Y2[clust,1:D2], center = TRUE, scale = TRUE)
        
        Y1.new.scaled.list[[j]] <- scale(Ytemp1[,1:D1], center = attr(obj.t1,"scaled:center"), scale = (attr(obj.t1,"scaled:scale")))
        Y2.new.scaled.list[[j]]<- scale(Ytemp2[,1:D2], center = attr(obj.t2,"scaled:center"), scale = (attr(obj.t2,"scaled:scale")))
        
      }
      
      for (j in (kminus+1):(kminus+2)) {
        obj.t1 <- scale(Y1[,1:D1], center = TRUE, scale = TRUE)
        obj.t2 <- scale(Y2[,1:D2], center = TRUE, scale = TRUE)
        Y1.new.scaled.list[[j]] <- scale(Ytemp1, center = attr(obj.t1,"scaled:center"), scale = (attr(obj.t1,"scaled:scale")))
        Y2.new.scaled.list[[j]] <- scale(Ytemp2, center = attr(obj.t2,"scaled:center"), scale = (attr(obj.t2,"scaled:scale")))
      }
      
    }
    
    
    for(l in 1:N.new) {
      
      for (j in 1:kminus){
        posteriortime[j,l] <-  (regy1.tmp$sigma2[active[j]]^-1 *(regy1.tmp$beta0[active[j]] + regy1.tmp$betahat[active[j],1:D1] %*% Y1.new.scaled.list[[j]][l,1:D1]) +  regy2.tmp$sigma2[active[j]]^-1 *(regy2.tmp$beta0[active[j]] + regy2$betahat[active[j],1:D2] %*% Y1.new.scaled.list[[j]][l,1:D1]) ) / (regy1.tmp$sigma2[active[j]]^-1 + regy2.tmp$sigma2[active[j]]^-1)
        
        posteriortimeweight[j,l] <- log(g[active[j]])  +  dMVN(as.vector(t(Ytemp1[l,1:D1])), mean = gmmx1.tmp$mu[active[j],1:D1], Q = gmmx1.tmp$S[active[j],1:D1,1:D1], log = TRUE)  + dMVN(as.vector(t(Ytemp2[l,1:D2])), mean = gmmx2.tmp$mu[active[j],1:D2], Q = gmmx2.tmp$S[active[j],1:D2,1:D2], log =TRUE) 
      }
      
      
      res <- try(dMVN(x = as.vector(t(Ytemp1[l,])), mean = gmmx1.tmp$mu[active[kminus+1],1:D1],  Q = gmmx1.tmp$S[active[kminus+1],1:D1,1:D1]) *  dMVN(x = as.vector(t(Ytemp2[l,])), mean = gmmx2.tmp$mu[active[kminus+1],1:D2],  Q = gmmx2.tmp$S[active[kminus+1],1:D2,1:D2]), silent =TRUE )
      if (class(res) == "try-error"){
        posteriortime[kminus+1,l] <- 0
        posteriortimeweight[kminus+1,l] <- -Inf
      } else{
        posteriortime[kminus+1,l] <-  ( regy1.tmp$sigma2[active[kminus+1]]^-1 *(regy1.tmp$beta0[active[kminus+1]] + regy1.tmp$betahat[active[kminus+1],1:D1] %*% Y1.new.scaled.list[[kminus+1]][l,1:D1] ) +  regy2.tmp$sigma2[active[kminus+1]]^-1 *(regy2.tmp$beta0[active[kminus+1]] + regy2.tmp$betahat[active[kminus+1],1:D2] %*% Y2.new.scaled.list[[kminus+1]][l,1:D1]) ) / (regy1.tmp$sigma2[active[kminus+1]]^-1 + regy2.tmp$sigma2[active[kminus+1]]^-1)
        posteriortimeweight[kminus+1,l] <-   log(alpha) +  dMVN(x = as.vector(t(Ytemp1[l,])), mean = gmmx1.tmp$mu[active[kminus+1],1:D1],  Q = gmmx1.tmp$S[active[kminus+1],1:D1,1:D1], log =TRUE) +  dMVN(x = as.vector(t(Ytemp2[l,])), mean = gmmx2.tmp$mu[active[kminus+1],1:D2],  Q = gmmx2.tmp$S[active[kminus+1],1:D2,1:D2], log= TRUE)
      }
      
      res2 <- try(dMVN(x = as.vector(t(Ytemp1[l,])), mean = gmmx1.tmp$mu[active[kminus+2],1:D1],  Q = gmmx1.tmp$S[active[kminus+2],1:D1,1:D1]) *  dMVN(x = as.vector(t(Ytemp2[l,])), mean = gmmx2.tmp$mu[active[kminus+2],1:D2],  Q = gmmx2.tmp$S[active[kminus+2],1:D2,1:D2]), silent =TRUE )
      if (class(res) == "try-error"){
        posteriortime[kminus+2,l] <- 0
        posteriornormtime[kminus+2,l] <- -Inf
      } else{
        posteriortime[kminus+2,l] <-  (regy1.tmp$sigma2[active[kminus+2]]^-1 *(regy1.tmp$beta0[active[kminus+2]] + regy1.tmp$betahat[active[kminus+2],1:D1] %*% Y1.new.scaled.list[[kminus+2]][l,1:D1] ) +  regy2.tmp$sigma2[active[kminus+2]]^-1 *(regy2.tmp$beta0[active[kminus+2]] + regy2.tmp$betahat[active[kminus+2],1:D2] %*% Y2.new.scaled.list[[kminus+2]][l,1:D1]) ) / (regy1.tmp$sigma2[active[kminus+2]]^-1 + regy2.tmp$sigma2[active[kminus+2]]^-1)
        posteriortimeweight[kminus+2,l] <-   log(alpha) +  dMVN(x = as.vector(t(Ytemp1[l,])), mean = gmmx1.tmp$mu[active[kminus+2],1:D1],  Q = gmmx1.tmp$S[active[kminus+2],1:D1,1:D1], log =TRUE) +  dMVN(x = as.vector(t(Ytemp2[l,])), mean = gmmx2.tmp$mu[active[kminus+2],1:D2],  Q = gmmx2.tmp$S[active[kminus+2],1:D2,1:D2], log= TRUE)
      }
      
      weights[,l] <- exp(posteriortimeweight[,l])/sum(exp(posteriortimeweight[,l]))
    }
    
    for ( l in 1:N.new){
      post.time[l,count] <- as.numeric(t(posteriortime[,l]) %*% weights[,l])
    } 
    
    cind[count] <- as.numeric(survConcordance(Surv(exp(time.new),censoring.new) ~ exp(-post.time[,count]))[1]) 
    # print(cind[count])
    
    #   Sys.sleep(0.1)
      setTxtProgressBar(pb, count)
    #     
    
  }
  
  #### To calculate average values over MCMC samples
  post.time.avg <<- apply(post.time[,1:Nps],1,mean)
  predCIndex.sbc <<- cind
  
  
}
