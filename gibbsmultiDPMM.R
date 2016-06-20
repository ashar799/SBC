## Gibbs Iterations for the multi view DPMM 

gibbsmultiDPMM = function(){
  
  source('priorPARAMETERS.R')
  param <- NA
  paramtime1 <- NA
  paramtime2 <- NA
  cognate <- NA
  hypercognate1 <- NA
  hypercognate2 <- NA
  loglike<- rep(0, iter)  
  est.regy1 <- list(0)
  est.regy2 <- list(0)
  est.gmmx1 <- list(0)
  est.gmmx2 <- list(0)
  c.list <- list(0)
  That.list <- list(0)
  alpha.list <- list(0)
  
  
  randy <- c(0)
  likli <- c(0)
  print("GIBB'S SAMPLING")
  pb <- txtProgressBar(min = 1, max = iter , style = 3)
  count = 1
  #################### GIBBS ITERATION ###################################################
  
  for (o in 1:iter) {
    
    
    ################## PARAMETERS OF THE DP Mixture Model ######################################################
    ## Updating the parameters based on the observations 
    source('posteriorGMM.R')
    param <- posteriorGMMparametrs(c,Y1,gmmx1$mu,gmmx1$S, alpha, K, gmmx1$epsilon, gmmx1$W, gmmx1$beta, gmmx1$ro,N,D1 )
    gmmx1$mu <- param$mean
    gmmx1$S <- param$precision
    param2 <- posteriorGMMparametrs(c,Y2,gmmx2$mu,gmmx2$S, alpha,K, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro,N,D2 )
    gmmx2$mu <- param2$mean
    gmmx2$S <- param2$precision
    
    
    source('multiposteriorAFT.R')
    paramtime2 <- posteriortimeparameterspenalized(c,Y2, That, regy2$lambda2, regy2$tau2, regy2$sigma2, regy2$beta0, regy2$betahat, K, gmmx2$epsilon, gmmx2$W, beta, ro, r, si, sig2.data,N, D2)
    regy2$beta0 <- paramtime2$beta0
    regy2$betahat <- paramtime2$betahat
    regy2$sigma2 <- paramtime2$sigma2
    regy2$lambda2 <- paramtime2$lambda2
    regy2$tau2 <- paramtime2$tau2
    
    
    
    paramtime1 <- posteriortimeparameterspenalized(c,Y1, That, regy1$lambda2, regy1$tau2, regy1$sigma2, regy1$beta0, regy1$betahat, K, gmmx1$epsilon, gmmx1$W, beta, ro, r, si, sig2.data,N, D1)
    regy1$beta0 <- paramtime1$beta0
    regy1$betahat <- paramtime1$betahat
    regy1$sigma2 <- paramtime1$sigma2
    regy1$lambda2 <- paramtime1$lambda2
    regy1$tau2 <- paramtime1$tau2
    
    
    
    
    
    
    ########################## THE HYPERPARAMETERS OF THE GMM #################################  
    source('posteriorhyperGMM.R')  
    # Updating the hyper paramters for the first data set
    hypercognate <- posteriorhyperPLUS(c, Y1, gmmx1$mu, gmmx1$S, gmmx1$epsilon, gmmx1$W, gmmx1$beta, gmmx1$ro )
    gmmx1$epsilon <- hypercognate$epsilon
    tmpW <- hypercognate$W
    gmmx1$W <- matrix(as.matrix(tmpW),nrow = D1, ncol =D1)
    gmmx1$ro <- hypercognate$ro
    
    
    
    ##Updating the hyper parameter for the second data set
    hypercognate2 <- posteriorhyperPLUS(c, Y2, gmmx2$mu, gmmx2$S, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro )
    gmmx2$epsilon <- hypercognate2$epsilon
    tmpW2 <- hypercognate2$W
    gmmx2$W <- matrix(as.matrix(tmpW2),nrow = D2, ncol =D2)
    gmmx2$ro <- hypercognate2$ro
    
    
    ### Updating Beta parameter for the first view #################
    #   source('posteriorbeta.R')
    #   if( o%%10 == 0){
    #     res <- try(posteriorbeta(c, gmmx1$beta, D1, gmmx1$S, gmmx1$W))
    #     if (class(res) == "try-error"){
    #       gmmx1$beta = gmmx1$beta
    #     } else{
    #       gmmx1$beta <- posteriorbeta(gmmx1$beta, D1, gmmx1$S, gmmx1$W)
    #       
    #     }
    #   } 
    #   ### Updating Beta parameter for the second view #################
    #   source('posteriorbeta.R')
    #   if( o%%10 == 0){
    #     res <- try(posteriorbeta(c, gmmx2$beta, D2, gmmx2$S, gmmx2$W))
    #     if (class(res) == "try-error"){
    #       gmmx2$beta = gmmx2$beta
    #     } else{
    #       gmmx2$beta <- posteriorbeta(gmmx2$beta, D2, gmmx2$S, gmmx2$W)
    #       
    #     }
    #   } 
    #   
    
    
    ################# INDICATOR VARIABLE ##################################################################
    ## Updating the indicator variables and the parameters
    source('multiposteriorCLASS.R') 
    cognate <- multiposteriorchineseAFT(c,Y1,Y2,D1,D2,That, K, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
    c <- cognate$c
    gmmx1 <- cognate$gmmx1
    gmmx2 <- cognate$gmmx2
    regy1 <- cognate$regy1
    regy2 <- cognate$regy2
    
    
    
    ########################### The Concentration Parameter #################################################################
    
    
    source('posterioralpha.R') 
    # Updating the concentration parameter
    alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
    
    
    
    ######################## The Censored Times ###########################################################
    source('multiupdatetime.R')
    # Updating the Time Variable
    ti <- NA
    ti <- multiupdatetime(c, Y1, Y2, Time,That, regy1, regy2)
    That <- ti$time
    
    
    
    if(o%% iter.thin == 0 ){
      est.regy1[[count]] <- regy1
      est.regy2[[count]] <- regy2
      est.gmmx1[[count]] <- gmmx1
      est.gmmx2[[count]] <- gmmx2
      c.list[[count]] <- c
      That.list[[count]] <- That
      alpha.list[[count]] <-  alpha
      count <- count +1
    }
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, o)
  } 
  
  assign("est.gmmx1", est.gmmx1, envir = .GlobalEnv)
  assign("est.gmmx2", est.gmmx2, envir = .GlobalEnv)
  assign("est.regy1", est.regy1, envir = .GlobalEnv)
  assign("est.regy2", est.regy2, envir = .GlobalEnv)
  assign("c.list", c.list, envir = .GlobalEnv)
  assign("alpha.list", alpha.list, envir = .GlobalEnv)
  assign("That.list", alpha.list, envir = .GlobalEnv)
  
  
  
}
