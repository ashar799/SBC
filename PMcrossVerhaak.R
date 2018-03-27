##################################################################################################
#### THE SBC model in Cross Validation setting ####################################################
##### This file Runs a SINGLE RUN FOR the cross validation 5 times cross validation on the Verhaak Cancer Data Set ####
setwd("~/Dropbox/Code/DPmixturemodel/SBC")
rm(list =ls())
source('import.R')

## Set the initial conditions
### The Cross-vaidation repeat number
u =1
### The fold number
v =1

### Load the data
############ Load Training set ###################################
load('/home/bit/ashar/ExpressionSets/ONE_VIEW/Verhark/OriginalVerhaakData.RData')
Y.full <- exprs.norm
mode(Y.full) <- "numeric"

###### Load Verhaark Gene Signature #################################
signature.prelim <- read.xlsx('/home/bit/ashar/ExpressionSets/ONE_VIEW/Verhark/ClaNC840_centroids.xls', sheetIndex =1)
signature.vk <- signature.prelim[3:nrow(signature.prelim),1]
verhaak.signature <- signature.vk[signature.vk %in% colnames(Y.full)]

##### Load Pheno Data for the training data set #####################
load('/home/bit/ashar/ExpressionSets/ONE_VIEW/Verhark/phenoVerhaark.RData')


######## Getting survival times and status #####################
time.pre <- as.numeric((pheno[,3]))
censoring.pre <- pheno[,2]
c.verhaak <- pheno[,4]
levels(c.verhaak)[1:4] <- c(1:4)

############################# PARAMETERS for GIBB's SAMPLING ####
iter = 200
iter.burnin = 200
iter.thin  = 5
Nps = as.integer(iter/ iter.thin)

### Some parameters for the DP mixture model ###################
k = 2
F =k


##Define the key functions to be used
source('rchinese.R')
source('priorPARAMETERS.R')
source('posteriorCensoredTime.R')
source('loadSBCvijversignature.R')
source('posteriorCLASS.R')
source('posteriorGMM.R')
source('posteriorAFT.R')
source('posteriorCensoredTime.R')
source('priorPARAMETERS.R')
source('calculateLIKELIHOOD.R')
source('posteriorhyperGMM.R')
source('posterioralpha.R')

#############################################################################################
########################## BEGIN THE CROSS-VALIDATION FOLD ###################################
##############################################################################################
set.seed(42*u)
#### Define the folds
folds <- createFolds(c.verhaak, k = 5, list = TRUE, returnTrain = FALSE)
test.index <- folds[[v]]
####################

####### If we define the splits ourselves ####
############### Defining Our data without prefiltering #########################
Y.pre.train <- Y.full[-test.index, ]
Y.pre.test <-  Y.full[test.index,]
c.true <- c.verhaak[-test.index]
c.true.new <- c.verhaak[test.index]
time <- time.pre[-test.index]
censoring <- censoring.pre[-test.index]
time.new <- time.pre[test.index]
censoring.new <- censoring.pre[test.index]


### setting up the survival objects
smod <-  Surv(exp(time), censoring)
smod.new <- Surv(exp(time.new), censoring.new)


#### Use Verhaak's classification to come up with separability and c-Index##
#### Using Verhaak singature and classification training and Testing Data Sets ##########################
Y.verhaak.train <- Y.pre.train[,colnames(Y.pre.train) %in% verhaak.signature]
Y.verhaak.test <- Y.pre.test[,colnames(Y.pre.test) %in% verhaak.signature]

###
recov.logrank.verhaak<- unlist(survdiff(smod ~ c.true))$chisq
pred.logrank.verhaak  <-  unlist(survdiff(smod.new ~ c.true.new))$chisq

####### Recovering C-Indexes ##############
### Using L1 penalized Cox-PH ############
linear.vv.recovery <- c(0)
linear.vv.prediction <- c(0)
for ( q in 1:4){
  ind <- which(c.true == q)
  ind.new <- which(c.true.new == q)
  reg.cox <- cv.glmnet(x = Y.verhaak.train[ind,], y = smod[ind], family = "cox", maxit = 10000000)
  linear.vv.recovery[ind] <- predict(object =reg.cox, newx = Y.verhaak.train[ind,], s= "lambda.min")
  linear.vv.prediction[ind.new] <- predict(object =reg.cox, newx = Y.verhaak.test[ind.new,], s= "lambda.min")
}
recovCIndex.vv.pcox <- as.numeric(survConcordance(smod ~ linear.vv.recovery)[1])
predCIndex.vv.pcox <- as.numeric(survConcordance(smod.new ~ linear.vv.prediction)[1])


### Use k-NN to predict using the VIJVER classification
label.train <- as.factor(c.true)
### One has to to tune the k-NN classifier for k ###
fitControl <- trainControl(method = "repeatedcv", number = 5,repeats = 5)
### Tune the parameter k 
knnFit <- train(x = Y.verhaak.train, y = label.train, method =  "knn", trControl = fitControl, tuneLength = 5)
knnPredict <- predict(knnFit,newdata = Y.verhaak.test)
label.test <- knnPredict
linear.vv.kk.prediction <- c(0)
for ( q in 1:2){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.cox <- cv.glmnet(x = Y.verhaak.train[ind,], y = smod[ind], family = "cox", maxit = 10000000)
  linear.vv.kk.prediction[ind.new] <- predict(object =reg.cox, newx = Y.verhaak.test[ind.new,], s= "lambda.min")
}
predCIndex.vv.kk.pcox<- as.numeric(survConcordance(smod.new ~ linear.vv.kk.prediction)[1])
pred.logrank.vv.kk <- unlist(survdiff(smod.new ~ label.test))$chisq




#### Use SBC signature ##
source('loadSBCverhaaksignature.R')
rel <- loadSBCverhaaksignature(Y.pre.train, time, censoring)
signature.sbc <- rel$signature.sbc


###### We get the signature and then we can define the folds #####
########## Creating Training and Test Data ##########################
Y <- Y.pre.train[,colnames(Y.pre.train) %in% signature.sbc ]
Y.new <- Y.pre.test[,colnames(Y.pre.test) %in% signature.sbc ] 

### Setting the important parameters
D <- ncol(Y)
N <- nrow(Y)
N.new <- nrow(Y.new)
K <- as.integer(N/5)

####### Use HC clustering + k-NN prediction  #########################################
###### Use SBC signature        ##################################################
#### Use k-means clustering 
gr.km <- kmeans(Y, F, nstart =10)
label.train <- as.factor(gr.km$cluster)
### One has to to tune the k-NN classifier for k ###
fitControl <- trainControl(method = "repeatedcv", number = 5,repeats = 5)
### Tune the parameter k 
knnFit <- train(x = Y, y = label.train, method =  "knn", trControl = fitControl, tuneLength = 5)
knnPredict <- predict(knnFit,newdata = Y.new )
label.test <- knnPredict

recov.logrank.kk <- unlist(survdiff(smod ~ label.train))$chisq
pred.logrank.kk <- unlist(survdiff(smod.new ~ label.test))$chisq



###### penCox on top of the clustering ###################################################################
######## Penalized Cox PH with k-means clustering###########################################
linear.kk.recovery <- c(0)
linear.kk.prediction <- c(0)
for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.cox <- cv.glmnet(x = Y[ind,], y = smod[ind], family = "cox", maxit = 10000000)
  linear.kk.recovery[ind] <- predict(object =reg.cox, newx = Y[ind,], s= "lambda.min")
  linear.kk.prediction[ind.new] <- predict(object =reg.cox, newx = Y.new[ind.new,], s= "lambda.min")
}
recovCIndex.kk.pcox <- as.numeric(survConcordance(smod ~ linear.kk.recovery)[1])
predCIndex.kk.pcox  <- as.numeric(survConcordance(smod.new ~ linear.kk.prediction)[1])

# ######################### Initialize the Parameters ##############################
Time <- cbind(time, censoring)
## HYPER PRIORS
## Hyper parameters of the DP
shape.alpha <- 2
rate.alpha <- 1
## Hyperparameters for the GMM
beta  = D+1
ro = 0.5



alpha  = rgamma(1, shape = shape.alpha, rate = rate.alpha )
c <-  rchinese(N,alpha)
f <- table(factor(c, levels = 1:max(c)))

## Empirical Bayes Estimate of the Hyperparameters
epsilon = as.vector(apply(Y,2,mean))
W = diag(diag(cov(Y)))


## Initialization of the parameters for Gaussian Mixture
mu = matrix(data = NA, nrow = K, ncol = D)
S = array(data = NA, dim =c(K,D,D))


#Sparsity controlling hyperparameter of the BAYESIAN LASSO MODEL
r =1
si = 1.78


##Actual parameters
lambda2 <- numeric(K)
tau2 = matrix(data = NA, nrow = K, ncol = D)
betahat = matrix(data = NA, nrow = K, ncol = D)
sigma2 <- rep(NA, K)
beta0 <- rep(NA, K)
That <-  c(0)

## Fitting a linear model to the whole model
Ysc <- scale(Y[1:N,1:D], center = TRUE, scale =TRUE)
lm.data <- lm(time ~ Ysc)
sig2.dat <-  var(lm.data$residuals)


## Set Some Initial Values for the Cluster Parameters


disclass <- table(factor(c, levels = 1:K))
activeclass <- which(disclass!=0)
for ( j in 1:length(activeclass)){
  
  priorone <- priordraw(beta, W, epsilon, ro, r, si, N, D, sig2.dat)
  mu[activeclass[j],] <- (priorone$mu)
  S[activeclass[j],1:D,1:D]  <- priorone$Sigma
  beta0[activeclass[j]] <- priorone$beta0
  sigma2[activeclass[j]] <- priorone$sigma2
  betahat[activeclass[j],1:D] <- priorone$betahat
  lambda2[activeclass[j]] <- priorone$lambda2
  tau2[activeclass[j], 1:D] <- priorone$tau2
}

# The Time has to be initialized
ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
That <- ti$time


################# K-Means BLASSO INITIALIZATION ############################################
G <- F

k.data <- kmeans(Y,F,nstart =10)

c <- k.data$cluster

c.kmeans <- c
#### Under special cases
###### c <- c.true


prior.numclust <- table(factor(c, levels = 1:K))
prior.activeclass <- which(prior.numclust!=0)

### The means are set using the k-means
for ( i in 1:length(prior.activeclass)){
  mu[prior.activeclass[i],1:D] <-  k.data$centers[i,1:D]
  
  S[prior.activeclass[i],1:D,1:D] <-  priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)$Sigma
  
  lclust <- which(c == prior.activeclass[i])
  
  reg.blas <- 0
  
  sum <- c(0)
  
  coeff <- 0
  
  Ytemp <-  matrix(NA, nrow = length(lclust), ncol = D)
  
  Ytemp <- scale(Y[lclust,1:D], center = TRUE, scale = TRUE)
  
  
  ### Part where I use the MONOMVN PACKAGE
  
  Ttemp <- as.vector(That[lclust])
  
  ntemp <- length(lclust)
  
  reg.blas <- blasso(Ytemp, Ttemp, T = 300,thin =  50, RJ = TRUE, mprior = 0.0 ,normalize = TRUE, verb = 0)
  
  sum <- summary(reg.blas, burnin= 100)
  
  ## Selecting those features which are relevant
  
  coeff <- unlist(lapply(strsplit(sum$coef[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
  
  beta0[prior.activeclass[1]] <- coeff[1]
  
  indexplusone <- D+1
  
  ind <- 2:indexplusone
  
  betahat[prior.activeclass[i], ] <- coeff[ind]
  
  ta <- unlist(lapply(strsplit(sum$tau2i[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
  
  tau2[prior.activeclass[i],] <- ta
  
  sigma2[prior.activeclass[i]] <- sum$s2[3]
  
  lambda2[prior.activeclass[i]] <- sum$lambda2[3]
  
}

## Deleting those values which are no longer relevant
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

############# Train the Model #########################################

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
  
  ####### If At all the W get's NA
  if (sum(is.na(diag(W))+ 0) > 0){
    W <- diag(diag(cov(Y)))
  }
}


# ### Gibb's Sampling
# ############## GIBBS SAMPLING WITH THINNING ######################################################
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
  
  likli.gibbs[o] <- loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)
  
  
  setTxtProgressBar(pb, o)
  Sys.sleep(0.1)
  
  
}





# ########## Analyze the fit ##########################################
# ### Good feature selection from heatmap plus cindex plus randindex
########## ANLAYSING THE MCMC samples AND CALCULATING METRICES #######################################################
source('linearprediction.R')
count <- Nps
final.rand <- c(0)

############ The Matrices that will store the results #################################################
##### Class Assignments ########################
c.matrix <- matrix(NA, nrow = N, ncol = count)
recovCIndex.sbc <- c(0)
recovCIndex.sbc.paft <- c(0)
################ Begin Analysig the MCMC samples #######################################################

require(Hmisc)

for (h in 1:Nps){
  ### See C-Index (concordance index) returned by the model
  surv.aft <- Surv(exp(time),censoring)
  ### Predict Time from the model
  tem.tim <- predicttime(c.list[[h]], Y, That.list[[h]], Time, beta0.list[[h]], betahat.list[[h]], sigma2.list[[h]])$predicttime
  library(Hmisc)
  recovCIndex.sbc[h] <-  survConcordance(surv.aft ~ exp(-tem.tim))[[1]]
  
  ## try to fit a Cluster-specific AFT model on top of it ##
  
  ###### penAFT ###################################################################
  ######## Penalized AFT with k-means clustering ######################################################
  sbc.aft <- c(0)
  for ( q in 1:F){
    ind <- which((c.list[[h]]) == q)
    L= length(ind)
    
    time.tmp <- time[ind]
    censoring.tmp <- censoring[ind]
    Y.tmp <- Y[ind,]
    
    reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
    coeff.pred <- coef(object =reg, newx = Y.tmp, s= "lambda.min")
    sbc.aft[ind] <- predict(object = reg, newx = Y.tmp, s = "lambda.min")
  }
  recovCIndex.sbc.paft[[h]] <- as.numeric(survConcordance(smod ~ exp(-sbc.aft))[1])
  
  c.matrix[,h] <- c.list[[h]]
  
}




###############################################
###### Calculating POINT ESTIMATES ############
###############################################
psm <- comp.psm(t(c.matrix))
mpear <- maxpear(psm)

### If we build a cluster specific sbc approach
c.final <- mpear$cl
c.sbc <- mpear$cl

active <- as.numeric(rownames(table(c.final)))


############ Time Covariate Slopes FOR Relevant Clusters ############
list.betahat <- list(0)

count = Nps

for ( i in 1:count){
  list.betahat[[i]] <- (betahat.list[[i]][active,] != 0) +0
}

Q <- length(active)
matrix.betahat <- array(data = NA, dim =c(Q,count,D))

for ( z in 1:Q){
  for ( x  in 1:count){
    matrix.betahat[z,x,] <- list.betahat[[x]][z,]
  }
}

final.betahat <- apply(matrix.betahat,c(1,3),mean)
### Probability of betahat of genes FOR ONE SIMULATION
##colnames(final.betahat) =  c(rep("relevant",rel.D),rep("irrelevant",irrel.D))
heatmapdata <- as.data.frame(final.betahat)
heatmap.2(t(as.matrix(heatmapdata)),dendrogram="none", margins=c(6,10), main = "Posterior prob. \n for Selection \n ", cexCol = 0.85, cexRow = 0.7, Rowv = FALSE)

recov.logrank.sbc <- unlist(survdiff(smod ~ c.sbc))$chisq


######## Predict CLASS MEMBERSHIP on  New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
N.new <- nrow(Y.new)
c.new.list <- list(0)
## The number of posterior samples

print("GOING THROUGH MCMC Samples")
pb <- txtProgressBar(min = 1, max = Nps , style = 3)

ctemp.new <- c(0)
modelweights <- c(0)

for (count in 1:Nps){
  
  ## Assign the parameters to the posterior sample
  ctemp <- c.list[[count]]
  mu <- mu.list[[count]]
  S <- S.list[[count]]
  
  g <- table(factor(ctemp, levels = 1:K))
  activeclass <- which(g!=0)
  ## The table function helps converting the data point specific indicator variables to class specific indicator variables
  kminus <- length(activeclass)
  # active <- activeclass
  #Two Auxilary Variables
  #The name of the auxilary variables are taken to be one and two more than the maximum value in the already active cluster set
  activeclass <- append(activeclass, max(activeclass)+1)
  activeclass <- append(activeclass, max(activeclass)+1)
  active <- activeclass
  
  
  
  ### Assigning values to parameters
  
  priortwo <- NA
  priortwo <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
  mu[active[kminus+1],1:D]  <- priortwo$mu
  S[active[kminus+1],1:D,1:D]  <- priortwo$Sigma[1:D,1:D]
  
  
  
  
  
  priorthree <- NA
  priorthree <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
  mu[active[kminus+2],1:D]  <- priorthree$mu
  S[active[kminus+2],1:D,1:D]  <- priorthree$Sigma[1:D,1:D]
  
  
  
  ###### Some quantities used to store probabilities
  posteriorweight <- matrix(0, nrow = length(active), ncol = N.new)
  weights <- matrix(0, nrow = length(active), ncol = N.new)
  weights.final <- c(0)
  ctemp.new <- c(0)
  #   ## This can't be parallelized !!!!!
  for(l in 1:N.new)  {
    
    ## Calculating the Expectations and also the normalization constant for the Expectation
    for (j in 1:kminus) {
      posteriorweight[j,l] <- log(g[active[j]]/ (N-1+alpha))  +  dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[j],1:D], Q = S[active[j],1:D,1:D], log =TRUE)
    }
    
    res <- try(dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]), silent=TRUE)
    if (class(res) == "try-error"){
      posteriorweight[kminus+1,l] <- -Inf
    } else{
      posteriorweight[kminus+1,l] <-  log(alpha/ (N-1+alpha)) +  dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D], log = TRUE)
    }
    
    res2 <- try(dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+2],1:D], Q= S[active[kminus+2],1:D,1:D]), silent=TRUE)
    if (class(res) == "try-error"){
      posteriorweight[kminus+2,l] <- -Inf
      
    } else{
      posteriorweight[kminus+2,l] <-  log(alpha/ (N-1+alpha))  +  dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+2],1:D], Q= S[active[kminus+2],1:D,1:D], log = TRUE)
    }
    
    weights[,l] <- exp(posteriorweight[,l])/sum(exp(posteriorweight[,l]))
    
    
    if (sum(exp(posteriorweight[,l])) < 1e-200){
      ctemp.new[l] <- sample(active, 1, prob = rep(1,length(active)), replace = TRUE)
    } else {
      ctemp.new[l] <- sample(active, 1, prob= weights[,l], replace = TRUE)
      
    }
    
    weights.final[l] <- dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[ctemp.new[l],1:D], Q = S[ctemp.new[l],1:D,1:D], log =TRUE)
    
  }
  
  
  modelweights[count] <- sum(weights.final)
  
  c.new.list[[count]] <- ctemp.new
  Sys.sleep(0.1)
  setTxtProgressBar(pb, count)
  
}

#
#

## Converting the list to a matrix

c.matrix.new <- matrix(NA, nrow = N.new, ncol = Nps)
for( h in 1:Nps){
  c.matrix.new[,h] <- c.new.list[[h]]
}

####### Use Ad-hoc methods to calculate the actual separability and C-Index ##################
### This function takes the posterior parameters AND predicts CLUSTER MEMEBERSHIP for the new points
### The posterior parameters being the RECOVERED CLUSTER MEMBERSHIP
#### THEN IT FITS A LOGISTIC/MULTINOMIAL regression FOR LINEAR SEPARABILITY
#### WE CAN ALSO EXPLORE NON-LINEAR Classifiers as k-SVM (RBF kernel)
### Check to see if this gives differences to state of the art

### POSSIBLE APPROACHES TO GET POINT ESTIMATES #####



### POSSIBILITY 1
### Use MPEAR APPROACH to get POINT ESTIMATES #############################################################################
psm2 <- comp.psm(t(c.matrix.new))
mpear2 <- maxpear(psm2)
### Generally the MPEAR output needs post-processing
c.sbc.new.mpear <- mpear2$cl




### POSSIBILITY 2
### Use Logistic regression to get the labels
reg <- cv.glmnet(x = Y, y = as.factor(c.sbc), family = "multinomial")
reg.new <- predict(object = reg, newx = Y.new, s = "lambda.min", type="class") 
c.sbc.new.log <- as.numeric(reg.new)


### POSSIBILITY 3
### Use k-Nearest neighbour
label.train <- as.factor(c.sbc)
### One has to to tune the k-NN classifier for k ###
fitControl <- trainControl(method = "repeatedcv", number = 5,repeats = 5)
### Tune the parameter k 
knnFit <- train(x = Y, y = label.train, method =  "knn", trControl = fitControl, tuneLength = 5)
c.sbc.new.knn <- predict(knnFit,newdata = Y.new )
#### Choose that configuration which has the highest difference in survival curves ####################

### POSSIBILITY 4
### Use k-svm 

library(e1071)
library(kernlab)
df <- data.frame(cbind(c.sbc,Y))
sig <- sigest(c.sbc ~.,data = df)
range.sig <-  seq(from = sig[1],to = sig[3], length.out = 25)
obj.base <- tune(svm, train.x = Y, train.y = factor(c.sbc), ranges = list(gamma = 2^(-10:14), cost = 2^(-10:14)), tunecontrol = tune.control(sampling = "cross"))

sigma.final <- obj.base$best.parameters$gamma
cost.final <- obj.base$best.parameters$cost

# Training the Classifier
ksvm.verk = ksvm(Y, factor(c.sbc), kernel = "rbfdot",kpar = list(sigma =sigma.final),C= cost.final, prob.model =TRUE)
## Predicting
labels.ksvm <- predict(ksvm.verk ,Y.new, type = "response")
pred.logrank.rbf <- survdiff(smod.new ~ labels.ksvm)
c.sbc.new.rbf <- labels.ksvm


### POSSIBILITY 5
### Use one of the actual cluster assignments

c.matrix.new.sbc <- matrix(0, nrow = N.new, ncol = Nps + 4)
c.matrix.new.sbc[,1:Nps] <- c.matrix.new[,1:Nps]
c.matrix.new.sbc[,(Nps+1)] <- c.sbc.new.mpear
c.matrix.new.sbc[,(Nps+2)] <- c.sbc.new.log
c.matrix.new.sbc[,(Nps+3)] <- c.sbc.new.knn 
c.matrix.new.sbc[,(Nps+4)] <- c.sbc.new.rbf


Nps.mod <- Nps +4
lr <- c(0)
for (j in 1:Nps.mod){
  lr[j] <-  1 - pchisq(unlist(survdiff(smod.new ~ c.matrix.new.sbc[,j]))$chisq,df = nlevels(as.factor(c.matrix.new.sbc[,i])) -1 )
}

c.sbc.new <- c.matrix.new.sbc[,which.min(lr)]


######## Predict C_INDEX on  New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
N.new <- nrow(Y.new)
c.new.list <- list(0)
## The number of posterior samples
#That.new <- time.new
post.time  = matrix(NA,nrow = nrow(Y.new), ncol = Nps)
print("GOING THROUGH MCMC Samples")
pb <- txtProgressBar(min = 1, max = Nps , style = 3)

cind <- c(0)

# modelweights <- c(0)

predCIndex.sbc <- c(0)


for (count in 1:Nps){
  
  ## Assign the parameters to the posterior sample
  ctemp <- c.list[[count]]
  mu <- mu.list[[count]]
  S <- S.list[[count]]
  beta0 <- beta0.list[[count]]
  betahat  <- betahat.list[[count]]
  sigma2  <- sigma2.list[[count]]
  g <- table(factor(ctemp, levels = 1:K))
  activeclass <- which(g!=0)
  ## The table function helps converting the data point specific indicator variables to class specific indicator variables
  kminus <- length(activeclass)
  # active <- activeclass
  #Two Auxilary Variables
  #The name of the auxilary variables are taken to be one and two more than the maximum value in the already active cluster set
  activeclass <- append(activeclass, max(activeclass)+1)
  activeclass <- append(activeclass, max(activeclass)+1)
  active <- activeclass
  
  
  
  ### Assigning values to parameters
  
  priortwo <- NA
  priortwo <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
  mu[active[kminus+1],1:D]  <- priortwo$mu
  S[active[kminus+1],1:D,1:D]  <- priortwo$Sigma[1:D,1:D]
  beta0[active[kminus+1]] <- priortwo$beta0
  sigma2[active[kminus+1]] <- priortwo$sigma2
  betahat[active[kminus+1],1:D] <- priortwo$betahat
  
  
  
  
  priorthree <- NA
  priorthree <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
  mu[active[kminus+2],1:D]  <- priorthree$mu
  S[active[kminus+2],1:D,1:D]  <- priorthree$Sigma[1:D,1:D]
  beta0[active[kminus+2]] <- priorthree$beta0
  sigma2[active[kminus+2]] <- priorthree$sigma2
  betahat[active[kminus+2],1:D] <- priorthree$betahat
  
  
  #######################################################
  ########### Scaling the data ########################
  Y.new.scaled.list <- list(0)
  for (j in 1:kminus) {
    clust <- which(ctemp == active[j])
    
    if(length(clust) > 1){
      obj.t <- scale(Y[clust,1:D], center = TRUE, scale = TRUE)
      Y.new.scaled.list[[j]] <- scale(Y.new, center = attr(obj.t,"scaled:center"), scale = (attr(obj.t,"scaled:scale")))
    } else {
      obj.t <- scale(Y[,1:D], center = TRUE, scale = TRUE)
      Y.new.scaled.list[[j]] <- scale(Y.new, center = attr(obj.t,"scaled:center"), scale = (attr(obj.t,"scaled:scale")))
      
    }
    
    
  }
  
  for (j in (kminus+1):(kminus+2)) {
    obj.t <- scale(Y[,1:D], center = TRUE, scale = TRUE)
    
    Y.new.scaled.list[[j]] <- scale(Y.new, center = attr(obj.t,"scaled:center"), scale = (attr(obj.t,"scaled:scale")))
  }
  
  
  ###### Some quantities used to store probabilities
  posteriortime <- matrix(0, nrow = length(active), ncol = N.new)
  posteriortimeweight <- matrix(0, nrow = length(active), ncol = N.new)
  weights <- matrix(0, nrow = length(active), ncol = N.new)
  
  
  
  ## This can't be parallelized !!!!!
  for(l in 1:N.new)  {
    
    ## Calculating the Expectations and also the normalization constant for the Expectation
    for (j in 1:kminus) {
      
      posteriortime[j,l] <-  beta0[active[j]] + betahat[active[j],1:D] %*% as.vector(t(Y.new.scaled.list[[j]][l,]))
      
      posteriortimeweight[j,l] <- log(g[active[j]]/ (N-1+alpha))  +  dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[j],1:D], Q = S[active[j],1:D,1:D], log =TRUE)
    }
    
    res <- try(dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]), silent=TRUE)
    if (class(res) == "try-error"){
      posteriortime[kminus+1,l] <- 0
      posteriortimeweight[kminus+1,l] <- -Inf
    } else{
      posteriortime[kminus+1,l] <-  beta0[active[kminus+1]] + betahat[active[kminus+1],1:D] %*% as.vector(t(Y.new.scaled.list[[j]][l,]))
      posteriortimeweight[kminus+1,l] <-  log(alpha/ (N-1+alpha)) +  dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D], log = TRUE)
    }
    
    res2 <- try(dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+2],1:D], Q= S[active[kminus+2],1:D,1:D]), silent=TRUE)
    if (class(res) == "try-error"){
      posteriortime[kminus+2,l] <- 0
      posteriortimeweight[kminus+2,l] <- -Inf
      
    } else{
      posteriortime[kminus+2,l] <- beta0[active[kminus+2]] + betahat[active[kminus+2],1:D] %*% as.vector(t(Y.new.scaled.list[[j]][l,]))
      posteriortimeweight[kminus+2,l] <-  log(alpha/ (N-1+alpha))  +  dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+2],1:D], Q= S[active[kminus+2],1:D,1:D], log = TRUE)
    }
    
    
    
    weights[,l] <- exp(posteriortimeweight[,l])/sum(exp(posteriortimeweight[,l]))
  }
  
  
  
  for ( l in 1:N.new){
    post.time[l,count] <- as.numeric(t(posteriortime[,l]) %*% weights[,l])
    
  }
  
  # modelweights[count] <- sum(exp((1/N.new) *apply(posteriortimeweight,1,sum)))
  
  predCIndex.sbc[count] <- as.numeric(survConcordance(smod.new ~ exp(-post.time[,count]))[1])
  #
  #    print(cind[count])
  
  ## Calculting the Model Weight
  
  Sys.sleep(0.1)
  setTxtProgressBar(pb, count)
  
  
}



#### Use PC method of the reviewer #####
pX <- prcomp(Y.pre.train)
pc.recov <- predict(pX,newdata = Y.pre.train)
pc.pred <- predict(pX,newdata = Y.pre.test)
cox.fit <- coxph(smod ~ ., data=as.data.frame(pc.recov[,1:20]))
linear.PC.recovery <- predict(cox.fit, newdata=as.data.frame(pc.recov[,1:20]))
linear.PC.prediction <- predict(cox.fit,  newdata= as.data.frame(pc.pred[,1:20]))

#### Getting the C-Indices ####
recovCIndex.PC <- as.numeric(survConcordance(smod ~ linear.PC.recovery)[1])
predCIndex.PC <- as.numeric(survConcordance(smod.new ~ linear.PC.prediction)[1])


filename <- paste('/home/bit/ashar/ExpressionSets/CROSS_VALIDATION/Verhaak/','repeat',u,'split',v,'.RData',sep = "")
save(list =ls(), file = filename)


#### Add these two extra methods
load("/home/bit/ashar/ExpressionSets/CROSS_VALIDATION/Verhaak/repeat1split1.RData")
#### Use Un-clustered Penalized L1 Cox Model
reg.full.cox <- cv.glmnet(x = Y.pre.train, y = smod, family = "cox", maxit = 10000000)
linear.full.recovery <- predict(object =reg.full.cox, newx = Y.pre.train, s= "lambda.min")
linear.full.prediction <- predict(object =reg.full.cox, newx = Y.pre.test, s= "lambda.min")
recovCIndex.NA.pcox <- as.numeric(survConcordance(smod ~ linear.full.recovery)[1])
predCIndex.NA.pcox <- as.numeric(survConcordance(smod.new ~ linear.full.prediction)[1])

#### Use Un-clustered Penalized L1 Cox Model
reg.sbc.cox <- cv.glmnet(x = Y, y = smod, family = "cox", maxit = 10000000)
linear.sbc.recovery<- predict(object =reg.sbc.cox, newx = Y, s= "lambda.min")
linear.sbc.prediction <- predict(object =reg.sbc.cox, newx = Y.new, s= "lambda.min")
recovCIndex.NAS.pcox <- as.numeric(survConcordance(smod ~ linear.sbc.recovery)[1])
predCIndex.NAS.pcox <- as.numeric(survConcordance(smod.new ~ linear.sbc.prediction)[1])

save(list =ls(), file = filename)




