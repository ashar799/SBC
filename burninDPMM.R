## Gibb's sampling 
burninDPMM = function(){
  
source('posteriorCLASS.R')
source('posteriorGMM.R')
source('posteriorAFT.R')
source('posteriorCensoredTime.R')
source('priorPARAMETERS.R')
source('calculateLIKELIHOOD.R')
source('posteriorhyperGMM.R')
source('posterioralpha.R')
source('posteriorbeta.R')



iter.burnin = iter.burnin
cognate <- NA
param <- NA
paramtime <- NA
loglike<- rep(0, iter)  
timeparam <- NA
time.predicted <- c(0)
cindex <- c(0)


likli <- c(0)
gmm.likli <- c(0)
aft.likli <- c(0)


cog <- loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)
likli[1] <- cog$loglikelihood 
gmm.likli[1] <- cog$GMMlikelihood
aft.likli[1] <- cog$AFTlikelihood

rmse <- c(0)
randy <- c(0)



o =1
#################### BURNIN PHASE ###################################################
o.iter = o
print("BURNIN...PHASE.. LOGLIKELIHOOD")

pb <- txtProgressBar(min = o.iter, max = iter.burnin , style = 3)
for (o in o.iter:iter.burnin) {
  
  
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  param <- posteriorGMMparametrs(c,Y,mu,S, alpha,K, epsilon, W, beta, ro,N,D )
  mu <- param$mean
  S <- param$precision
  paramtime <- posteriortimeparameters(c, That, lambda2,tau2,sigma2,beta0, betahat, Y, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.data)
  beta0 <- paramtime$beta0
  betahat <- paramtime$betahat
  sigma2 <- paramtime$sigma2
  lambda2 <- paramtime$lambda2
  tau2 <- paramtime$tau2
  
  ########################## THE HYPERPARAMETERS OF THE GMM #################################  
  #  Updating the hyper paramters
  hypercognate <- posteriorhyperPLUS (c, Y, mu, S, epsilon, W, beta, ro )
  epsilon <- hypercognate$epsilon
  tmpW <- hypercognate$W
  W <- matrix(as.matrix(tmpW),nrow = D, ncol =D)
  ro <- hypercognate$ro
  
#   if( o%%10 == 0){
#     res <- try(posteriorbeta(c, beta, D, S, W))
#     if (class(res) == "try-error"){
#       beta = beta
#     } else{
#       beta <- posteriorbeta(c, beta, D, S, W)
#       
#     }
#   } 
#   
  ################# INDICATOR VARIABLE ##################################################################
  ## Updating the indicator variables and the parameters
  cognate <- posteriorchineseAFT(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)
  c <- cognate$indicator
  mu <- cognate$mean
  S <- cognate$precision
  beta0 <- cognate$beta0
  betahat <- cognate$betahat
  sigma2 <- cognate$sigma2
  lambda2 <- cognate$lambda2
  tau2 <- cognate$tau2
  
  ########################### The Concentration Parameter #################################################################
  # Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  ####################### The Censored Times ###########################################################
  # Updating the Time Variable
  ti <- NA
  ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
  That <- ti$time
  
  
  ##################### Print SOME Statistics #####################################################
  randy[o] <- adjustedRandIndex(c.kmeans,as.factor(c))
  print(randy[o])
  cog <- loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)
  likli[o+1] <- cog$loglikelihood 
  gmm.likli[o+1] <- cog$GMMlikelihood
  aft.likli[o+1] <- cog$AFTlikelihood
  
  print(likli[o+1])
  print(gmm.likli[o+1])
  print(aft.likli[o+1])
  
  print(o/iter.burnin)
  
  
  ##### Print the status bar
  Sys.sleep(0.1)
  #setTxtProgressBar(pb, o)


 ####### If At all the W get's NA
  if (sum(is.na(diag(W))+ 0) > 0){
    W <- diag(diag(cov(Y)))
    }
} 

assign("alpha", alpha, envir = .GlobalEnv)
assign("ro", ro, envir = .GlobalEnv)
assign("That", That, envir = .GlobalEnv)
assign("c", c, envir = .GlobalEnv)
assign("epsilon", epsilon, envir = .GlobalEnv)
assign("W", W, envir = .GlobalEnv)
assign("mu", mu, envir = .GlobalEnv)
assign("S", S, envir = .GlobalEnv)
assign("beta0", beta0, envir = .GlobalEnv)
assign("betahat", betahat, envir = .GlobalEnv)
assign("sigma2", sigma2, envir = .GlobalEnv)
assign("lambda2", lambda2, envir = .GlobalEnv)
assign("tau2", tau2, envir = .GlobalEnv)
assign("randy.burnin", randy, envir = .GlobalEnv)
assign("rmse.burnin", rmse, envir = .GlobalEnv)
assign("likli.burnin", likli, envir = .GlobalEnv)
assign("gmm.burnin", gmm.likli, envir = .GlobalEnv)
assign("aft.burnin", aft.likli, envir = .GlobalEnv)

plot(likli, main = 'Burnin Iterations')
plot(gmm.likli, main = 'GMM Burnin Iterations')
plot(aft.likli, main = 'AFT Burnin Iterations')



}

