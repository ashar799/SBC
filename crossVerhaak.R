##################################################################################################
#### THE SBC model in Cross Validation setting ####################################################
##### This file Runs a cross validation 5 times cross validation on the Verhaak data set ####
rm(list =ls())
source('import.R')

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
time.pre <- as.numeric(exp(pheno[,3]))
censoring.pre <- pheno[,2]
c.verhaak <- pheno[,4]
levels(c.verhaak)[1:4] <- c(1:4)


#### Prepare for the parallel computation ###
#### Prepare for the parallel computation ###
library(parallel)
number.of.cores <- detectCores()
cluster.parallel <- makeCluster(number.of.cores -1, type="FORK")
library(doParallel)
registerDoParallel(cluster.parallel)
library(foreach)
getDoParWorkers()


############################# PARAMETERS for GIBB's SAMPLING ####
iter = 200
iter.burnin = 200
iter.thin  = 5
Nps = as.integer(iter/ iter.thin)

### Some parameters for the DP mixture model ###################
k = 4
F =k


#### Define key variables which store the results ####
# recov.logrank.verhaak <- c(0)
# pred.logrank.verhaak <- c(0)
# 
# recovCIndex.vv.pcox <- c(0)
# predCIndex.vv.pcox <- c(0)
# 
# predCIndex.vv.kk.pcox <- c(0)
# pred.logrank.vv.kk <- c(0)
# 
# 
# recov.logrank.kk <- c(0)
# pred.logrank.kk <- c(0)
# 
# 
# recovCIndex.kk.pcox <- c(0)
# predCIndex.kk.pcox <- c(0)
# 
# recov.logrank.sbc <- c(0)
# pred.logrank.sbc <- c(0)
# 
# recovCIndex.PC <- c(0)
# predCIndex.PC <- c(0)
# 
# 
# recovCIndex.sbc.matrix <- matrix(0, nrow =  Nps, ncol = 5)
# recovCIndex.sbc.paft.matrix <- matrix(0, nrow =  Nps, ncol = 5)
# 
# predCIndex.sbc.matrix <- matrix(0, nrow =  Nps, ncol = 5)
# predCIndex.sbc.pcox <- c(0)

# recovCIndex.PC <- c(0)
# predCIndex.PC <- c(0)

#signature.sbc.list <- list(0)


#############################################################################################
########################## BEGIN THE CROSS-VALIDATION FOLD ###################################
##############################################################################################

#############################################################################################
########################## BEGIN THE CROSS-VALIDATION FOLD ###################################
##############################################################################################

#### Define the folds
set.seed(42)
folds <- createFolds(c.verhaak, k = 5, list = TRUE, returnTrain = FALSE)



#### The actual cross-validation loop contains around 200 lines of code
### Begins here #####



results <- foreach(i = 1:5, .export = ls(),.packages = list.p) %dopar%{
  
  test.index <- folds[[i]]
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
  smod <-  Surv(time, censoring)
  smod.new <- Surv(time.new, censoring.new)
  
  
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
  
  
  ######################### Initialize the Parameters ##############################
  source('initialize.R')
  initialize()
  
  ########### Train the Model #########################################
  source('burninDPMM.R')
  burninDPMM()
  
  source('gibbsDPMM.R')
  gibbsDPMM()
  
  ########## Analyze the fit ##########################################
  ### Good feature selection from heatmap plus cindex plus randindex
  source('MCMCanalyze.R')
  MCMCanalyze()
  recov.logrank.sbc <- unlist(survdiff(smod ~ c.sbc))$chisq
  recovCIndex.sbc
  recovCIndex.sbc.paft
  
  ######## Predict CLASS MEMBERSHIP on  New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
  source('predictCLASS.R')
  predictCLASS(Y.new)
  
  ####### Use Ad-hoc methods to calculate the actual separability and C-Index ##################
  source('predictADHOCverhaak.R')
  pred.logrank.sbc <- unlist(survdiff(smod.new ~ c.sbc.new))$chisq
  
  
  ######## Predict C_INDEX on  New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
  source('predictTIME.R')
  predictchineseAFTtime(Y.new)
  predCIndex.sbc
  
  
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
  
  
  list('recov.logrank.verhaak'= recov.logrank.verhaak,
         'pred.logrank.verhaak' = pred.logrank.verhaak,
         'recovCIndex.vv.pcox' = recovCIndex.vv.pcox, 
         'predCIndex.vv.pcox'  = predCIndex.vv.pcox,
         'predCIndex.vv.kk.pcox' = predCIndex.vv.kk.pcox,
         'pred.logrank.vv.kk'  = pred.logrank.vv.kk,
         'recov.logrank.kk' = recov.logrank.kk,
         'pred.logrank.kk' = pred.logrank.kk,
         'recovCIndex.kk.pcox' = recovCIndex.kk.pcox,
         'predCIndex.kk.pcox' = predCIndex.kk.pcox,
         'recov.logrank.sbc' = recov.logrank.sbc,
         'pred.logrank.sbc' = pred.logrank.sbc,
         'recovCIndex.PC' = recovCIndex.PC,
         'predCIndex.PC' = predCIndex.PC,
         'recovCIndex.sbc' = recovCIndex.sbc,
         'recovCIndex.sbc.paft' = recovCIndex.sbc.paft,
         'predCIndex.sbc' = predCIndex.sbc,
          'signature.sbc' = signature.sbc)
  
  
  
}

stopImplicitCluster()

