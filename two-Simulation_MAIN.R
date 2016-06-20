## Testing the DPplusAFT model with Gibb's Sampling
### Script to test if the model works with higher dimensions
## With irrelevant features and with Time Censoring 

rm(list = ls())
setwd("~/Dropbox/Code/DPmixturemodel/SBC")

########################################### Generating Molecular data sets ##########################################################3
## Number of points
N = 200

## Number of Clusters
F = 2

## Percentage of Points belonging to different data
p.dist = c(0.3,0.7)


##### The following Code Generates the Molecular Data and the Survival Times
{
## Source the file
source('i-simulate.R')
D1 = 50
con<- generatedata(N,D1,F,p.dist)
Y1 <-  con$Y
betalist1 <- con$beta
timepur1 <- con$timepur
sdlist1 <- list(0)
for ( j in 1:F){
  sdlist1[[j]] <- sd(as.vector(timepur1[[j]]))
}

## Generate the second data sets
source('i-simulate.R')
D2 = 30
con2 <- generatedata(N,D2,F,p.dist)
Y2 <-  con2$Y
betalist2 <- con2$beta
timepur2 <- con2$timepur
sdlist2 <- list(0)
for ( j in 1:F){
  sdlist2[[j]] <- sd(as.vector(timepur2[[j]]))
}

## Generating Time Data

## Percentage of Noise/Overlap in Time Data
prob.noise = 0.05

## Total Percentage of censoring
prob.censoring = 0.10


## The pure time is generated using weighted addition of the individual curves
time.pur.final <- list(0)
for ( i in 1:F){
  time.pur.final[[i]] <- (sdlist1[[i]])^-2 * as.vector(timepur1[[i]]) + (sdlist2[[i]])^-2 * as.vector(timepur2[[i]])
}

## Checking the distributions with plotting them
dat <- data.frame(xx = c(as.vector(timepur1[[1]]),as.vector(timepur2[[1]]),as.vector(time.pur.final[[1]])),yy = rep(letters[1:3],each = as.integer(N*(p.dist[1]))))

ggplot(dat,aes(x=xx)) + 
  geom_density(data=subset(dat,yy == 'a'),fill = "blue", alpha = 0.2) +
  geom_density(data=subset(dat,yy == 'b'),fill = "green", alpha = 0.2) +
  geom_density(data=subset(dat,yy == 'c'),fill = "yellow", alpha = 0.2) 
    


## Simulating Time Data which is ONE dimensional from a mixture distribution of two clusters
time.cluster <- MixSim(MaxOmega = prob.noise, K = F, p = 1, int =c(0,3))


## time.noise.list containts the cluster specific noises and means ( cluster variance is inverse sum of the individiual covariances)
time.noise.list <- list(0)
for ( i in 1:F){
  time.noise.list[[i]] <- rnorm(as.integer(N * p.dist[i]), mean = time.cluster$Mu[i], sd = sqrt((sdlist1[[i]]^-2 + sdlist2[[i]]^-2 )^-1))
}

time.list <- list(0)
for ( i in 1:F){
  time.list[[i]] <- time.pur.final[[i]] + time.noise.list[[i]]
}

#################################################################
############# Individual Times Per Data Set ##################################
#################################################################


## It would be nice to compare the individual times with the joint times
time.noise.list1 <- list(0)
for ( i in 1:F){
  time.noise.list1[[i]] <- rnorm(as.integer(N * p.dist[i]), mean = time.cluster$Mu[i], sd = sdlist1[[i]])
}

time.noise.list2 <- list(0)
for ( i in 1:F){
  time.noise.list2[[i]] <- rnorm(as.integer(N * p.dist[i]), mean = time.cluster$Mu[i], sd = sdlist2[[i]])
}

time.list1 <- list(0)
for ( i in 1:F){
  time.list1[[i]] <- timepur1[[i]] + time.noise.list1[[i]]
}

time.list2 <- list(0)
for ( i in 1:F){
  time.list2[[i]] <- timepur2[[i]] + time.noise.list2[[i]]
}


#############################################################################
## True Labels for the points
c.true <- c(0)
for ( i in 1:F){
  c.true <- rbind(as.matrix(c.true) , as.matrix(c(rep(i, as.integer(N * p.dist[i])))))  
}
c.true <- as.factor(c.true[-1,])



## Check PCs
pc <- prcomp(Y1)
pc.pred <- predict(pc,newdata = Y1)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.true)






####################################### MAKING TIME from cluster data ########################################################
## Real time without Censoring
time.real <- c(0)
for (i in 1:F){
  time.real <- c(time.real, time.list[[i]])
} 
time.real <- time.real[-1]
time.real <- as.vector(unlist(time.real))


####################################### Adding CENSORING INFORMATION  ################################################################
## Adding the censoring information to the TIME
prob.censoring = 0.1
censoring <- rbinom(n = N, size =1, prob = 1- prob.censoring)
right.censoring.time <- min(time.real)  

time <- time.real

index.time <- which(censoring==0)
for ( q in 1:length(index.time)){
  time[index.time[q]] <- right.censoring.time
  
}


## Boxplots for Vizualization of the time Data without censoring
boxplot(time.list)
### A Quick ManWhittney U / Kruskal test  test to check if the time's of the two cluster are significantly different
## wilcox.test(as.vector(time.list[[1]]), as.vector(time.list[[2]]), alternative = "two.sided")
kruskal.test(time, c.true)

}

############################# PARAMETERS for GIBB's SAMPLING ######################################
iter = 50
iter.burnin = 100
iter.thin  =5

################################# GIBBS SAMPLING  ###################################################
#### The following code initializes the parameters of the two data sets
{
Time <- cbind(time, censoring) 

K = as.integer(N/2)

## HYPER PRIORS
## Hyper parameters of the DP
shape.alpha <- 2
rate.alpha <- 1
## Hyperparameters for the GMM
beta  = (D1 +D2)
ro = 0.5
## Initialize the c using chinese restaurant process
source('rchinese.R')
alpha  = rgamma(1, shape = shape.alpha, rate = rate.alpha )
c <-  rchinese(N,alpha)
f <- table(factor(c, levels = 1:max(c)))
#Sparsity controlling hyperparameter of the BAYESIAN LASSO MODEL
r =1
si = 1.78


### LETS MAKE A LIST "gmmx" to store parameters/hyperprameters for X and "regy" to store paameters for Regression Y
## For the First Data Set
gmmx1 <- list(0)
gmmx1$epsilon <-  as.vector(apply(Y1,2,mean))
gmmx1$W <- cov(Y1)
gmmx1$mu <- matrix(data = NA, nrow = K, ncol = D1)
gmmx1$S <-  array(data = NA, dim =c(K,D1,D1))

regy1 <- list(0)
regy1$lambda2 <- numeric(K)
regy1$tau2 = matrix(data = NA, nrow = K, ncol = D1)
regy1$betahat = matrix(data = NA, nrow = K, ncol = D1)
regy1$sigma2 <- rep(NA, K)
regy1$beta0 <- rep(NA, K)

## For the second data set
gmmx2 <- list(0)
gmmx2$epsilon <-  as.vector(apply(Y2,2,mean))
gmmx2$W <- cov(Y2)
gmmx2$mu <- matrix(data = NA, nrow = K, ncol = D2)
gmmx2$S <-  array(data = NA, dim =c(K,D2,D2))

regy2 <- list(0)
regy2$lambda2 <- numeric(K)
regy2$tau2 = matrix(data = NA, nrow = K, ncol = D2)
regy2$betahat = matrix(data = NA, nrow = K, ncol = D2)
regy2$sigma2 <- rep(NA, K)
regy2$beta0 <- rep(NA, K)


###### To initialize the parameters for all the data sets
That <-  time
####### We can use a simple Linear Model to get some estimates of the variance##########
Yg <- cbind(Y1,Y2)
Dg <- (D1 + D2)

## Fitting a linear model to the whole model
Ysc <- scale(Yg[1:N,1:Dg], center = TRUE, scale =TRUE)
lm.data <- lm(time ~ Ysc)
sig2.dat <-  var(lm.data$residuals)


## Set Some Initial Values for the Cluster Parameters
source('multi-initialize.R')
## For the first data set
cont1 <- multiinit(Y1,c,beta, gmmx1$W, gmmx1$epsilon, ro, r, si,N,D1, sig2.dat)
gmmx1$mu <- cont1$mu
gmmx1$S <- cont1$S
regy1$lambda2 <- cont1$lambda2
regy1$tau2 <- cont1$tau2
regy1$betahat <- cont1$betahat
regy1$sigma2 <- cont1$sigma2
regy1$beta0 <- cont1$beta0

## For the second data set
cont2 <- multiinit(Y2,c,beta, gmmx2$W, gmmx2$epsilon, ro, r, si,N,D2, sig2.dat)
gmmx2$mu <- cont2$mu
gmmx2$S <- cont2$S
regy2$lambda2 <- cont2$lambda2
regy2$tau2 <- cont2$tau2
regy2$betahat <- cont2$betahat
regy2$sigma2 <- cont2$sigma2
regy2$beta0 <- cont2$beta0



## Initialization part for the parmaters of AFT Model with k-means and Bayesian Lasso and Normal Bayesian Regression
source('multikmeansBlasso.R')
km <- multikmeansBlasso(c,Y1,Y2,D1,D2,That, F,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2 )
c <- km$c
gmmx1 <- km$gmmx1
gmmx2 <- km$gmmx2 
regy1 <- km$regy1
regy2 <- km$regy2


## Adjusted Initial Rand INDEX measure
randindexi <- adjustedRandIndex(c.true,as.factor(c))

## Initial Likelihood
source('multilikelihood.R')
likint <- multiloglikelihood( c,Y1,Y2,D1,D2,That,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)


}


############################################################################
### BURNIN AND GIBBS SAMPLING ##############################################
## Gibb's sampling #########################################################

cognate <- NA
param <- NA
paramtime <- NA
loglike<- rep(0, iter)  
timeparam <- NA
time.predicted <- c(0)
cindex <- c(0)
randy <- c(0)
likli <- c(0)

#################### BURNIN PHASE ###################################################
print("BURNIN...PHASE")
for (o in 1:iter.burnin) {
  
  
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  source('posteriorGMM.R')
  param <- posteriorGMMparametrs(c,Y1,gmmx1$mu,gmmx1$S, alpha, K, gmmx1$epsilon, gmmx1$W, beta, ro,N,D1 )
  gmmx1$mu <- param$mean
  gmmx1$S <- param$precision
  param2 <- posteriorGMMparametrs(c,Y2,gmmx2$mu,gmmx2$S, alpha,K, gmmx2$epsilon, gmmx2$W, beta, ro,N,D2 )
  gmmx2$mu <- param2$mean
  gmmx2$S <- param2$precision
  
  
  source('multiposteriorAFT.R')
  paramtime1 <- posteriortimeparameterspenalized(c,Y1, That, regy1$lambda2, regy1$tau2, regy1$sigma2, regy1$beta0, regy1$betahat, K, gmmx1$epsilon, gmmx1$W, beta, ro, r, si, sig2.data,N, D1)
  regy1$beta0 <- paramtime1$beta0
  regy1$betahat <- paramtime1$betahat
  regy1$sigma2 <- paramtime1$sigma2
  regy1$lambda2 <- paramtime1$lambda2
  regy1$tau2 <- paramtime1$tau2
  
  paramtime2 <- posteriortimeparameterspenalized(c,Y2, That, regy2$lambda2, regy2$tau2, regy2$sigma2, regy2$beta0, regy2$betahat, K, gmmx2$epsilon, gmmx2$W, beta, ro, r, si, sig2.data,N, D2)
  regy2$beta0 <- paramtime2$beta0
  regy2$betahat <- paramtime2$betahat
  regy2$sigma2 <- paramtime2$sigma2
  regy2$lambda2 <- paramtime2$lambda2
  regy2$tau2 <- paramtime2$tau2
  
  
  
  ########################## THE HYPERPARAMETERS OF THE GMM #################################  
  source('posteriorhyperGMM.R')  
   
  # Updating the hyper paramters for the first data set
  hypercognate <- posteriorhyperPLUS(c, Y1, gmmx1$mu, gmmx1$S, gmmx1$epsilon, gmmx1$W, beta, ro,D1 )
  gmmx1$epsilon <- hypercognate$epsilon
  tmpW <- hypercognate$W
  gmmx1$W <- matrix(as.matrix(tmpW),nrow = D1, ncol =D1)
  
  ##Updating the hyper parameter for the second data set
  hypercognate2 <- posteriorhyperPLUS(c, Y2, gmmx2$mu, gmmx2$S, gmmx2$epsilon, gmmx2$W, beta, ro,D2 )
  gmmx2$epsilon <- hypercognate2$epsilon
  tmpW2 <- hypercognate2$W
  gmmx2$W <- matrix(as.matrix(tmpW2),nrow = D2, ncol =D2)
      
    
  ################# INDICATOR VARIABLE ##################################################################
  ## Updating the indicator variables and the parameters
  source('multiposteriorCLASS.R')  
  cognate <- multiposteriorchineseAFT(c,Y1,Y2,D1,D2,That, F,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  c <- cognate$c
  gmmx1 <- cognate$gmmx1
  gmmx2 <- cognate$gmmx2
  regy1 <- cognate$regy1
  regy2 <- cognate$regy2
  
 
  
  ########################### The Concentration Parameter #################################################################
  source('posterioralpha.R') 
  # Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  ##################### Print SOME Statistics #####################################################
  randy[o] <- adjustedRandIndex(c.true,as.factor(c))
  print(randy[o])
  likli[o] <- multiloglikelihood( c,Y1,Y2,D1,D2,That, F,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  print(likli[o])
  print(o/iter.burnin)

} 







#####################   GIBBS SAMPLING WITH THINNING ######################################################

est.regy1 <- list(0)
est.regy2 <- list(0)
est.gmmx1 <- list(0)
est.gmmx2 <- list(0)
c.list <- list(0)





print("GIBB'S SAMPLING")
count = 1
for (o in 1:iter) {
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  source('posteriorGMMparametrs.R')
  param <- posteriorGMMparametrs(c,Y1,gmmx1$mu,gmmx1$S, alpha, K, gmmx1$epsilon, gmmx1$W, beta, ro,N,D1 )
  gmmx1$mu <- param$mean
  gmmx1$S <- param$precision
  param2 <- posteriorGMMparametrs(c,Y2,gmmx2$mu,gmmx2$S, alpha,K, gmmx2$epsilon, gmmx2$W, beta, ro,N,D2 )
  gmmx2$mu <- param2$mean
  gmmx2$S <- param2$precision
  
  
  source('posteriortimeparameterspenalized.R')
  paramtime1 <- posteriortimeparameterspenalized(c,Y1, That, regy1$lambda2, regy1$tau2, regy1$sigma2, regy1$beta0, regy1$betahat, K, gmmx1$epsilon, gmmx1$W, beta, ro, r, si, sig2.data,N, D1)
  regy1$beta0 <- paramtime1$beta0
  regy1$betahat <- paramtime1$betahat
  regy1$sigma2 <- paramtime1$sigma2
  regy1$lambda2 <- paramtime1$lambda2
  regy1$tau2 <- paramtime1$tau2
  
  paramtime2 <- posteriortimeparameterspenalized(c,Y2, That, regy2$lambda2, regy2$tau2, regy2$sigma2, regy2$beta0, regy2$betahat, K, gmmx2$epsilon, gmmx2$W, beta, ro, r, si, sig2.data,N, D2)
  regy2$beta0 <- paramtime2$beta0
  regy2$betahat <- paramtime2$betahat
  regy2$sigma2 <- paramtime2$sigma2
  regy2$lambda2 <- paramtime2$lambda2
  regy2$tau2 <- paramtime2$tau2
  
  
  
  ########################## THE HYPERPARAMETERS OF THE GMM #################################  
  source('posteriorhyper.R')  
  
  # Updating the hyper paramters for the first data set
  hypercognate <- posteriorhyper (c, Y1, gmmx1$mu, gmmx1$S, gmmx1$epsilon, gmmx1$W, beta, ro,D1 )
  gmmx1$epsilon <- hypercognate$epsilon
  tmpW <- hypercognate$W
  gmmx1$W <- matrix(as.matrix(tmpW),nrow = D1, ncol =D1)
  
  ##Updating the hyper parameter for the second data set
  hypercognate2 <- posteriorhyper (c, Y2, gmmx2$mu, gmmx2$S, gmmx2$epsilon, gmmx2$W, beta, ro,D2 )
  gmmx2$epsilon <- hypercognate2$epsilon
  tmpW2 <- hypercognate2$W
  gmmx2$W <- matrix(as.matrix(tmpW2),nrow = D2, ncol =D2)
  
  
  
  
  ################# INDICATOR VARIABLE ##################################################################
  ## Updating the indicator variables and the parameters
  source('multiposteriorchineseAFT.R')  
  cognate <- multiposteriorchineseAFT(c,Y1,Y2,D1,D2,That, F,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  c <- cognate$c
  gmmx1 <- cognate$gmmx1
  gmmx2 <- cognate$gmmx2
  regy1 <- cognate$regy1
  regy2 <- cognate$regy2
  
  
  
  ########################### The Concentration Parameter #################################################################
  
  
  source('posterioralpha.R') 
  # Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  ##################### Print SOME Statistics #####################################################
  randy[o] <- adjustedRandIndex(c.true,as.factor(c))
  print(randy[o])
  likli[o] <- multiloglikelihood( c,Y1,Y2,D1,D2,That, F,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  print(likli[o])
  print(o/iter.burnin)
  
  
  
  if(o%% iter.thin == 0 ){
   est.regy1[[count]] <- regy1
   est.regy2[[count]] <- regy2
   est.gmmx1[[count]] <- gmmx1
   est.gmmx2[[count]] <- gmmx2
   count <- count +1
  }
  
  
  
  
  print(o/iter) 
  #   print(loglike[o])
  #   print(cindex)
} 

########## ANLAYSING THE OUTPUT #######################################################

count <- count -1

## Selecting that clustering which gives the maximum POSTERIOR PROBABILITY
ri <- 1:count

index.good <- which.max(unlist(lapply(ri,function(x) adjustedRandIndex(c.true,as.factor(c.list[[x]])))))

## FINAL VALUES
c.final <- c.list[[index.good]]
mu.final <- mu.list[[index.good]] 
beta0.final <- beta0.list[[index.good]]
betahat.final <-   betahat.list[[index.good]]
sigma2.final <- sigma2.list[[index.good]]  
lambda2.final <- lambda2.list[[index.good]] 
That.final <- That.list[[index.good]]


## FINAL VALUES
## Adjusted Rand INDEX measure
randindex.final <- adjustedRandIndex(c.true,as.factor(c.final))


source('predicttime.R')
time.predicted.final <- predicttime(c.final,Y, That.final,Time, beta0.final, betahat.final, sigma2.final)$predicttime


## Prelimnary estimates of the RAND and C-INDEX index
source('calcindex.R')
cindex.final <- calcindex(c.final,Time,time.predicted.final)$cindex

## Calcuating the RandIndex and the C-index
ra <- c(0)
c1 <- c(0)
c2 <- c(0)


for ( i in 1:count){
  ra[i] <- adjustedRandIndex(c.true,as.factor(c.list[[i]]))
  time.pr <- predicttime(c.list[[i]],Y, That.list[[i]],Time, beta0.list[[i]], betahat.list[[i]], sigma2.list[[i]])$predicttime
  c1[i] <-  calcindex(c.list[[i]],Time,time.pr)$cindex[1]
  c2[i] <-  calcindex(c.list[[i]],Time,time.pr)$cindex[2]
}

pdf('/home/bit/ashar/Dropbox/WarsawTalk/Boxplots.pdf')
boxplot(ra,c1,c2, names = c("RandInd","C-Ind1","C-Ind2"), main = "Rand and Concordance Index for Simulation")
leg.text <- c("samples = 200", "dims = 50", "cluster =2", "relevant.dims = 4")
legend("topright", leg.text)
dev.off()




surv.ob <- Surv(Time[,1],Time[,2])
surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)


pdf('/home/bit/ashar/Dropbox/WarsawTalk/KaplanMeier.pdf')
plot(surv.fit, col = c("blue", "green"))
title("Kaplan-Meier Curves\nfor the Simulation")
leg.text <- c("LogRank Test p-value of 7.79e-09")
legend("topright", leg.text)
dev.off()
