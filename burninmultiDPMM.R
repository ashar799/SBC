## Burnin Iterations for the multi view DPMM 

burninmultiDPMM = function(){
  
source('priorPARAMETERS.R')
source('multilikelihood.R')
param <- NA
paramtime1 <- NA
paramtime2 <- NA
cognate <- NA
hypercognate1 <- NA
hypercognate2 <- NA
loglike<- rep(0, iter)  




burnin.likli <- c(0)
gmm.likli <- c(0)
aft.likli <- c(0)
randy <- c(0)

#################### BURNIN PHASE ###################################################
print("BURNIN...PHASE")
for (o in 1:iter.burnin) {
  
  
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
  paramtime2 <- posteriortimeparameterspenalized(c,Y2, That, regy2$lambda2, regy2$tau2, regy2$sigma2, regy2$beta0, regy2$betahat, K, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro, r, si, sig2.data,N, D2)
  regy2$beta0 <- paramtime2$beta0
  regy2$betahat <- paramtime2$betahat
  regy2$sigma2 <- paramtime2$sigma2
  regy2$lambda2 <- paramtime2$lambda2
  regy2$tau2 <- paramtime2$tau2
  
  
  paramtime1 <- posteriortimeparameterspenalized(c,Y1, That, regy1$lambda2, regy1$tau2, regy1$sigma2, regy1$beta0, regy1$betahat, K, gmmx1$epsilon, gmmx1$W, gmmx1$beta, gmmx1$ro, r, si, sig2.data,N, D1)
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
  hypercognate2 <- posteriorhyperPLUS(c, Y2, gmmx2$mu, gmmx2$S, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro)
  gmmx2$epsilon <- hypercognate2$epsilon
  tmpW2 <- hypercognate2$W
  gmmx2$W <- matrix(as.matrix(tmpW2),nrow = D2, ncol =D2)
  gmmx2$ro <- hypercognate2$ro
  
  

  
  
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
  
  
  
  ####################### The Censored Times ###########################################################
  source('multiupdatetime.R')
  # Updating the Time Variable
  ti <- NA
  ti <- multiupdatetime(c, Y1, Y2, Time,That, regy1, regy2)
  That <- ti$time
  
  
  
  
  ##################### Print SOME Statistics #####################################################
  randy[o] <- adjustedRandIndex(c.kmeans,as.factor(c))
  print(randy[o])
  cg <- multiloglikelihood(c,Y1,Y2,D1,D2,That,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  burnin.likli[o] <- cg$loglikelihood
  gmm.likli[o] <- cg$GMMlikelihood
  aft.likli[o] <- cg$AFTlikelihood 
  
  print(burnin.likli[o])
  print(gmm.likli[o])
  print(aft.likli[o])
  print(o/iter.burnin)
  
} 

assign("alpha", alpha, envir = .GlobalEnv)
assign("gmmx1", gmmx1, envir = .GlobalEnv)
assign("gmmx2", gmmx2, envir = .GlobalEnv)
assign("regy1", regy1, envir = .GlobalEnv)
assign("regy2", regy2, envir = .GlobalEnv)
assign("c", c, envir = .GlobalEnv)
assign("burnin.likli", burnin.likli, envir = .GlobalEnv)
assign("gmm.likli", gmm.likli, envir = .GlobalEnv)
assign("aft.likli", gmm.likli, envir = .GlobalEnv)


plot(gmm.likli, main = 'GMM Burnin Iterations')
plot(aft.likli, main = 'AFT Burnin Iterations')
plot(burnin.likli, main = 'Overall Burnin Iterations')
}

