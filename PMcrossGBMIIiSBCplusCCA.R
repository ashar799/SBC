##### This file runs Poor Man's cross-validation results on the i-SBC model ######
##### The i-SBC code with CCA pre-processing  is tested ##########################
setwd("~/Dropbox/Code/DPmixturemodel/SBC")
rm(list =ls())
source('import.R')

## Set the initial conditions
### The Cross-vaidation repeat number
u =1
### The fold number
v =1

##############################################################################################
#### Load the two data sets ##################################################################
##### Load Pheno Data and Load Verhaak Data Sets  ############################################
load('/home/bit/ashar/ExpressionSets/ONE_VIEW/Verhark/phenoVerhaark.RData')
load("/home/bit/ashar/ExpressionSets/TWO_VIEW/TCGA_GBM/IntegratedDataSurv.rda")
load("/home/bit/ashar/ExpressionSets/TWO_VIEW/TCGA_GBM/mRNA.rda")
load("/home/bit/ashar/ExpressionSets/TWO_VIEW/TCGA_GBM/miRNA.rda")

### Selecting the 189 patients which have common type of genomic data ###
patient.tcga <- rownames(mRNA)[(rownames(mRNA) %in% rownames(miRNA)) & (rownames(mRNA) %in% pheno$Patient)]


###############################################################
Y.mrna.tcga <- mRNA[patient.tcga,]
Y.mirna.tcga <- miRNA[patient.tcga,]

######## Getting survival times and status #####################
######## Lets Get the Pheno Data #####################################################################
pheno.ord <- pheno[match(patient.tcga,pheno$Patient),]
pheno.ord$Patient == rownames(Y.mrna.tcga)
pheno.ord$Patient == rownames(Y.mirna.tcga)

time.pre <- as.numeric(pheno.ord[,3])
censoring.pre <- pheno.ord[,2]
c.verhaak <- pheno.ord[,4]
levels(c.verhaak)[1:4] <- c(1:4)

############################# PARAMETERS for GIBB's SAMPLING ####
iter = 200
iter.burnin = 200
iter.thin  = 5
Nps = as.integer(iter/ iter.thin)


#############################################################################################
########################## BEGIN THE CROSS-VALIDATION FOLD ###################################
##############################################################################################
set.seed(42*u)
#### Define the folds
folds <- createFolds(c.verhaak, k = 5, list = TRUE, returnTrain = FALSE)
test.index <- folds[[v]]


###############################################################################
####### If we define the splits ourselves #####################################
############### Defining Our data without prefiltering #########################
Y1.pre.train <- Y.mrna.tcga[-test.index, ]
Y2.pre.train <- Y.mirna.tcga[-test.index,]
Y1.pre.test <-  Y.mrna.tcga[test.index,]
Y2.pre.test <-  Y.mirna.tcga[test.index,]

c.true <- c.verhaak[-test.index]
c.true.new <- c.verhaak[test.index]

time <- time.pre[-test.index]
censoring <- censoring.pre[-test.index]
time.new <- time.pre[test.index]
censoring.new <- censoring.pre[test.index]

#### Use SBC signature ##
source('loadiSBCGBMIIsignature.R')
rel <- loadiSBCGBMIIsignature(Y1.pre.train, Y2.pre.train, time, censoring)
signature1.sbc <- rel$signature1.sbc
signature2.sbc <- rel$signature2.sbc

###### We get the signature and then we can define the folds #####
########## Creating Training and Test Data ##########################
Y1 <- Y1.pre.train[,colnames(Y1.pre.train) %in% signature1.sbc ]
Y2 <- Y2.pre.train[,colnames(Y2.pre.train) %in% signature2.sbc ]

Y1.new <- Y1.pre.test[,colnames(Y1.pre.test) %in% signature1.sbc ]
Y2.new <- Y2.pre.test[,colnames(Y2.pre.test) %in% signature2.sbc]

###### If using CCA features on top
#### Use SBC signature ##
source('loadiSBCGBMIIsignatureCCA.R')
rel <- loadiSBCGBMIIsignatureCCA(Y1, Y2, Y1.new, Y2.new)
Y1 <- rel$Y1
Y2 <- rel$Y2
Y1.new <- rel$Y1.new
Y2.new <- rel$Y2.new


###################################################
D1 <- ncol(Y1)
D2 <- ncol(Y2)

#### For other methods ###########################
Y <- as.matrix(cbind(Y1,Y2))
Y.new <- as.matrix(cbind(Y1.new,Y2.new))
D <- D1 + D2

############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 200
iter.thin  = 2
k = 2
F =k
N <- nrow(Y1)
K <-  as.integer(N/5)
Time <- cbind(time,censoring)


########################### Sourcing key files #################
source('rchinese.R')
source('multiinit.R')
source('multilikelihood.R')
source('priorPARAMETERS.R')
source('multilikelihood.R')
source('multikmeansBlasso.R')
source('posteriorGMM.R')
source('multiposteriorAFT.R')
source('posteriorhyperGMM.R') 
source('posterioralpha.R') 
source('multiposteriorCLASS.R')
source('multiupdatetime.R')


########################################################################################################
############################ COMPARISON only with k-means + K-Nearest Neighbours #######################
########################################################################################################

### Comparison with base methods
smod <-  Surv(exp(time), censoring)
smod.new <- Surv(exp(time.new), censoring.new)


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
  reg.cox <- cv.glmnet(x = Y[ind,1:D], y = Surv(time[ind],censoring[ind]), family = "cox", maxit = 10000000)
  linear.kk.recovery[ind] <- predict(object =reg.cox, newx = Y[ind,], s= "lambda.min")
  linear.kk.prediction[ind.new] <- predict(object =reg.cox, newx = Y.new[ind.new,], s= "lambda.min")
}
recovCIndex.kk.pcox <- as.numeric(survConcordance(smod ~ linear.kk.recovery)[1])
predCIndex.kk.pcox  <- as.numeric(survConcordance(smod.new ~ linear.kk.prediction)[1])


############### Use Principal Components without clustering ###############

#### Use PC method of the reviewer #####
Y.pre.train <- as.matrix(cbind(Y1.pre.train, Y2.pre.train))
Y.pre.test <-  as.matrix(cbind(Y1.pre.test, Y2.pre.test))

pX <- prcomp(Y.pre.train)
pc.recov <- predict(pX,newdata = Y.pre.train)
pc.pred <- predict(pX,newdata = Y.pre.test)
cox.fit <- coxph(smod ~ ., data=as.data.frame(pc.recov[,1:20]))
linear.PC.recovery <- predict(cox.fit, newdata=as.data.frame(pc.recov[,1:20]))
linear.PC.prediction <- predict(cox.fit,  newdata= as.data.frame(pc.pred[,1:20]))

#### Getting the C-Indices ####
recovCIndex.PC <- as.numeric(survConcordance(smod ~ linear.PC.recovery)[1])
predCIndex.PC <- as.numeric(survConcordance(smod.new ~ linear.PC.prediction)[1])



##################################################
####### Initialize ###############################
##################################################
##GIBBS SAMPLING INITIALIZATION  ################
#################################################
## HYPER PRIORS


## Hyper parameters of the DP
shape.alpha <- 2
rate.alpha <- 1


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
gmmx1$W <- diag(diag(cov(Y1)))
gmmx1$mu <- matrix(data = NA, nrow = K, ncol = D1)
gmmx1$S <-  array(data = NA, dim =c(K,D1,D1))
gmmx1$ro <- 0.5
gmmx1$beta <- D1 +1




regy1 <- list(0)
regy1$lambda2 <- numeric(K)
regy1$tau2 = matrix(data = NA, nrow = K, ncol = D1)
regy1$betahat = matrix(data = NA, nrow = K, ncol = D1)
regy1$sigma2 <- rep(NA, K)
regy1$beta0 <- rep(NA, K)

## For the second data set
gmmx2 <- list(0)
gmmx2$epsilon <-  as.vector(apply(Y2,2,mean))
gmmx2$W <- diag(diag(cov(Y2)))
gmmx2$mu <- matrix(data = NA, nrow = K, ncol = D2)
gmmx2$S <-  array(data = NA, dim =c(K,D2,D2))
gmmx2$ro <- 0.5
gmmx2$beta <- D2 +1



regy2 <- list(0)
regy2$lambda2 <- numeric(K)
regy2$tau2 = matrix(data = NA, nrow = K, ncol = D2)
regy2$betahat = matrix(data = NA, nrow = K, ncol = D2)
regy2$sigma2 <- rep(NA, K)
regy2$beta0 <- rep(NA, K)


###### To initialize the parameters for all the data sets
That <-  time


## Fitting a linear model to the whole model
Ysc <- scale(Y[1:N,1:D], center = TRUE, scale =TRUE)
lm.data <- lm(time ~ Ysc)
sig2.dat <-  var(lm.data$residuals)


## Set Some Initial Values for the Cluster Parameters


## For the first data set
cont1 <- multiinit(Y1,c, gmmx1$beta, gmmx1$W, gmmx1$epsilon, gmmx1$ro, r, si,N,D1, sig2.dat)
gmmx1$mu <- cont1$mu
gmmx1$S <- cont1$S
regy1$lambda2 <- cont1$lambda2
regy1$tau2 <- cont1$tau2
regy1$betahat <- cont1$betahat
regy1$sigma2 <- cont1$sigma2
regy1$beta0 <- cont1$beta0

## For the second data set
cont2 <- multiinit(Y2,c, gmmx1$beta, gmmx2$W, gmmx2$epsilon, gmmx2$ro, r, si,N,D2, sig2.dat)
gmmx2$mu <- cont2$mu
gmmx2$S <- cont2$S
regy2$lambda2 <- cont2$lambda2
regy2$tau2 <- cont2$tau2
regy2$betahat <- cont2$betahat
regy2$sigma2 <- cont2$sigma2
regy2$beta0 <- cont2$beta0


## Initialization part for the parmaters of AFT Model with k-means and Bayesian Lasso and Normal Bayesian Regression
km <- multikmeansBlasso(c,Y1,Y2,D1,D2,That,K, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2,surv.obj )
c <- km$c
c.kmeans <- c

gmmx1 <- km$gmmx1
gmmx2 <- km$gmmx2 
regy1 <- km$regy1
regy2 <- km$regy2


## Initial Likelihood
likli.int <- multiloglikelihood( c,Y1,Y2,D1,D2,That,K, beta, ro, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)

###############################################################################################################
######################## Train the Model ######################################################################
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
  
  param <- posteriorGMMparametrs(c,Y1,gmmx1$mu,gmmx1$S, alpha, K, gmmx1$epsilon, gmmx1$W, gmmx1$beta, gmmx1$ro,N,D1 )
  gmmx1$mu <- param$mean
  gmmx1$S <- param$precision
  param2 <- posteriorGMMparametrs(c,Y2,gmmx2$mu,gmmx2$S, alpha,K, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro,N,D2 )
  gmmx2$mu <- param2$mean
  gmmx2$S <- param2$precision
  
  
  
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
  cognate <- multiposteriorchineseAFT(c,Y1,Y2,D1,D2,That, K, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  c <- cognate$c
  gmmx1 <- cognate$gmmx1
  gmmx2 <- cognate$gmmx2
  regy1 <- cognate$regy1
  regy2 <- cognate$regy2
  
  
  
  ########################### The Concentration Parameter #################################################################
  
  # Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  
  ####################### The Censored Times ###########################################################
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
  param <- posteriorGMMparametrs(c,Y1,gmmx1$mu,gmmx1$S, alpha, K, gmmx1$epsilon, gmmx1$W, gmmx1$beta, gmmx1$ro,N,D1 )
  gmmx1$mu <- param$mean
  gmmx1$S <- param$precision
  param2 <- posteriorGMMparametrs(c,Y2,gmmx2$mu,gmmx2$S, alpha,K, gmmx2$epsilon, gmmx2$W, gmmx2$beta, gmmx2$ro,N,D2 )
  gmmx2$mu <- param2$mean
  gmmx2$S <- param2$precision
  
  
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
  cognate <- multiposteriorchineseAFT(c,Y1,Y2,D1,D2,That, K, r, si,sig2.dat,gmmx1, gmmx2, regy1, regy2)
  c <- cognate$c
  gmmx1 <- cognate$gmmx1
  gmmx2 <- cognate$gmmx2
  regy1 <- cognate$regy1
  regy2 <- cognate$regy2
  
  
  
  ########################### The Concentration Parameter #################################################################
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



#######################################################################################################
############################ Analyzing the MCMC output #################################################
#######################################################################################################
Y <- as.matrix(cbind(Y1,Y2))
Y.new <- as.matrix(cbind(Y1.new, Y2.new))
Y.test <- Y.new



#######################################################################
############# TRAININIG DATA ###########################################
#######################################################################
########################################################################
#### This function calculates some important metrices for the TRAINING DATA Data
#### C-Index
#### Point Estimate of Clsuter Assignments based on m-pear
########## ANLAYSING THE MCMC samples AND CALCULATING METRICES #######################################################

Nps = as.integer(iter/ iter.thin)
count <- Nps


############ The Matrices that will store the results #################################################
cindex.final1 <- c(0)
cindex.final2 <- c(0)   
recovCIndex.isbc <- c(0)  
recovCIndex.isbc.paft <- c(0)
################ Begin Analysig the MCMC samples #######################################################



for (h in 1:Nps){
  ### Adjusted Rand Indices
  surv.aft <- Surv(exp(time),censoring)
  
  ### Predict Time from the model
  source('linearprediction.R')
  tem.tim1 <- as.vector(unlist(predicttime(c.list[[h]], Y1, That, Time, est.regy1[[h]]$beta0, est.regy1[[h]]$betahat, est.regy1[[h]]$sigma2)))
  tem.tim2 <- as.vector(unlist(predicttime(c.list[[h]], Y2, That, Time,  est.regy2[[h]]$beta0, est.regy2[[h]]$betahat, est.regy2[[h]]$sigma2)))
  
  cindex.final1[h] <-  survConcordance(surv.aft ~ exp(-tem.tim1))[[1]]
  cindex.final2[h] <-  survConcordance(surv.aft ~ exp(-tem.tim2))[[1]]
  
  source('multilinearprediction.R')
  tem.tim <- as.vector(unlist(multipredictlinear(c.list[[h]],  est.regy1[[h]], est.regy2[[h]] )))
  
  recovCIndex.isbc[h] <-  survConcordance(surv.aft ~ exp(-tem.tim))[[1]]
  
  
  
  ###### penAFT ###################################################################
  ######## Penalized AFT with k-means clustering ######################################################
  isbc.aft <- c(0)
  for ( q in 1:F){
    ind <- which((c.list[[h]]) == q)
    L= length(ind)
    
    time.tmp <- time[ind]
    censoring.tmp <- censoring[ind]
    Y.tmp <- Y[ind,]
    
    reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
    coeff.pred <- coef(object =reg, newx = Y.tmp, s= "lambda.min")
    isbc.aft[ind] <- predict(object = reg, newx = Y.tmp, s = "lambda.min") 
  }
  recovCIndex.isbc.paft[h] <- as.numeric(survConcordance(smod ~ exp(-isbc.aft))[1])
}



###############################################
###### Calculating POINT ESTIMATES ############
###############################################
##### Class Assignments ########################
c.matrix <- matrix(NA, nrow = N, ncol = count)
for ( i in 1:count){
  c.matrix[,i] <- c.list[[i]]
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


############ Time Covariate Slopes FOR Relevant Clusters and Heatmap Plots ############
list.betahat1 <- list(0)

for ( i in 1:count){
  list.betahat1[[i]] <- (est.regy1[[i]]$betahat[active,] != 0) +0
}

Q <- length(active)
matrix.betahat1 <- array(data = NA, dim =c(Q,count,D1))

for ( z in 1:Q){
  for ( x  in 1:count){
    matrix.betahat1[z,x,] <- list.betahat1[[x]][z,]
  }
}

final.betahat1 <- apply(matrix.betahat1,c(1,3),mean)
### Probability of betahat of genes FOR ONE SIMULATION
##colnames(final.betahat) =  c(rep("relevant",rel.D),rep("irrelevant",irrel.D))
heatmapdata1 <- as.data.frame(final.betahat1)
#heatmap.2(t(as.matrix(heatmapdata1)),dendrogram="none", col =cm.colors(180), margins=c(6,10), main = "Posterior prob. \n for Selection for Data set 1 ", cexCol = 0.85, cexRow = 0.7, Rowv = FALSE)


########################## For the second data set ####################################################
list.betahat2 <- list(0)

for ( i in 1:count){
  list.betahat2[[i]] <- (est.regy2[[i]]$betahat[active,] != 0) +0
}

Q <- length(active)
matrix.betahat2 <- array(data = NA, dim =c(Q,count,D2))

for ( z in 1:Q){
  for ( x  in 1:count){
    matrix.betahat2[z,x,] <- list.betahat2[[x]][z,]
  }
}

final.betahat2 <- apply(matrix.betahat2,c(1,3),mean)
heatmapdata2 <- as.data.frame(final.betahat2)
#heatmap.2(t(as.matrix(heatmapdata2)),dendrogram="none", col =cm.colors(180), margins=c(6,10), main = "Posterior prob. \n for Selection for Data set 2 ", cexCol = 0.85, cexRow = 0.7, Rowv = FALSE)

############# Calculating the recovery logrank statistic
recov.logrank.sbc <- unlist(survdiff(smod ~ c.sbc))$chisq


##### Predicting the Class of the new data points ##################
Y1.test <- Y1.new
Y2.test <- Y2.new

N.new <- nrow(Y1.test)
c.new.list <- list(0)
## The number of posterior samples
Nps <- as.integer(iter/ iter.thin)
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
reg <- cv.glmnet(x = Y, y = as.factor(c.sbc), family = "binomial")
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

predCIndex.sbc <- c(0)

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
  
  predCIndex.sbc[count] <- as.numeric(survConcordance(smod.new ~ exp(-post.time[,count]))[1])
  
  
  #   Sys.sleep(0.1)
  setTxtProgressBar(pb, count)
  #     
  
}


filename <- paste('/home/bit/ashar/ExpressionSets/CROSS_VALIDATION/GBMIICCA/','repeat',u,'split',v,'.RData',sep = "")
save(list =ls(), file = filename)

#### Add these two extra methods
load("/home/bit/ashar/ExpressionSets/CROSS_VALIDATION/GBMIICCA/repeat1split1.RData")
#### Use Un-clustered Penalized L1 Cox Model
reg.sbc.cox <- cv.glmnet(x = Y, y = smod, family = "cox", maxit = 10000000)
linear.sbc.recovery<- predict(object =reg.sbc.cox, newx = Y, s= "lambda.min")
linear.sbc.prediction <- predict(object =reg.sbc.cox, newx = Y.new, s= "lambda.min")
recovCIndex.NAS.pcox <- as.numeric(survConcordance(smod ~ linear.sbc.recovery)[1])
predCIndex.NAS.pcox <- as.numeric(survConcordance(smod.new ~ linear.sbc.prediction)[1])

save(list =ls(), file = filename)



