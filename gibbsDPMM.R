gibbsDPMM = function(){
  

############## GIBBS SAMPLING WITH THINNING ######################################################
nmrse <- c(0)
mu.list <- list(0)
S.list <- list(0)
beta0.list <- list(0)
betahat.list <- list(0) 
sigma2.list <- list(0)
lambda2.list <- list(0)
tau2.list <- list(0)
c.list <- list(0)
That.list <- list(0)
likli.gibbs <- c(0)
W.list <- list(0)


source('posteriorCLASS.R')
source('posteriorGMM.R')
source('posteriorAFT.R')
source('posteriorCensoredTime.R')
source('priorPARAMETERS.R')
source('calculateLIKELIHOOD.R')
source('posteriorhyperGMM.R')
source('posterioralpha.R')
source('posteriorbeta.R')


print("GIBB'S SAMPLING")
pb <- txtProgressBar(min = 1, max = iter , style = 3)
count = 1
for (o in 1:iter) {
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
  
  
  ######################## The Censored Times ###########################################################
  # Updating the Time Variable
  ti <- NA
  ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
  That <- ti$time
  
  
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  
  if(o%% iter.thin == 0 ){
    mu.list[[count]] <- mu
    S.list[[count]] <- S
    W.list[[count]] <- W
    beta0.list[[count]] <- beta0
    betahat.list[[count]] <- betahat  
    sigma2.list[[count]] <- sigma2
    lambda2.list[[count]] <- lambda2
    tau2.list[[count]] <- tau2
    c.list[[count]] <- c
    That.list[[count]] <- That
    
#     nmrse[count] <- calcrmse(time.real,time.predicted)$rmse
    count <- count +1
  }
  
# likli.gibbs[o] <- loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)

  
  
  #   print(o/iter) 
  #   print(loglike[o])
  #   print(cindex)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, o)
} 


assign("c.list", c.list, envir = .GlobalEnv)
assign("mu.list", mu.list, envir = .GlobalEnv)
assign("S.list", S.list, envir = .GlobalEnv)
assign("beta0.list", beta0.list, envir = .GlobalEnv)
assign("betahat.list", betahat.list, envir = .GlobalEnv)
assign("sigma2.list", sigma2.list, envir = .GlobalEnv)
assign("lambda2.list", lambda2.list, envir = .GlobalEnv)
assign("tau2.list", tau2.list, envir = .GlobalEnv)
assign("That.list", That.list, envir = .GlobalEnv)
assign("sigma2.list", sigma2.list, envir = .GlobalEnv)
assign("W.list", W.list, envir = .GlobalEnv)




}
