########## This file compares prediction on test cases with different methods
####### So we use k-means + SVM to predict labels on new data points ###########

TESTCOMPARISON = function(){
  

############ Predicting New Class Labels using SVM #################################  
# gr.km <- kmeans(Y, F, nstart =10)
# label.train <- gr.km$cluster
# svms <- sapply(2^(-10:14), function(cost) cross(ksvm(Y, factor(label.train), C=cost, kernel="vanilladot", kpar=list(), cross=5)))
# mysvm <- ksvm(Y, factor(label.train), C=2^(-10:14)[which.min(svms)], kernel="vanilladot", kpar=list(), cross=10) # accuracy ~97%
# pred.svm <-  predict(mysvm, Y.new)
# predRandIndex.svm <<- adjustedRandIndex(c.true.new, pred.svm)


############ Predicting the New Class Labels using kNN #############################
gr.km <- kmeans(Y, F, nstart =10)
label.train <- gr.km$cluster
knear <- knn(train = Y, test = Y.new, cl = label.train, k = F)
predRandIndex.knear <<- adjustedRandIndex(c.true.new, knear)

######## Predicting New C-Indices based on a Penalized Cox or AFT model#################### 
######## Penalized Cox PH ###########################################

linear.pred.cox <- c(0)
### see if we can use glmnet
reg.pcox <- cv.glmnet(x = Y, y = Surv(exp(time), censoring), family = "cox")
linear.pred.cox <- predict(object =reg.pcox, newx = Y.new, s= "lambda.min")
smod <-  Surv(exp(time.new), censoring.new)
predCIndex.cox <<- as.numeric(survConcordance(smod ~ linear.pred.cox)[1])

##### Penalized AFT Model #############################################
linear.pred.paft <- c(0)
### see if we can use glmnet
reg.paft <- cv.glmnet(x = Y, y = time, family = "gaussian")
linear.pred.paft <- predict(object = reg.paft, newx = Y.new, s= "lambda.min")
smod <-  Surv(exp(time.new), censoring.new)
predCIndex.aft <<- as.numeric(survConcordance(smod ~ exp(-linear.pred.paft))[1])

  
###### K-Means + KNN 

gr.km <- kmeans(Y, F, nstart =10)
label.train <- gr.km$cluster
label.test <- knn(train = Y, test = Y.new, cl = label.train, k = F)


#### penAFT ###################################################################

linear.kkpaft.prediction <- c(0)
for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.aft <- cv.glmnet(x = Y[ind,], y = time[ind], family = "gaussian")
  linear.kkpaft.prediction[ind.new] <- predict(object =reg.aft, newx = Y.new[ind.new,], s= "lambda.min")
}
predCIndex.kkpaft <<- as.numeric(survConcordance(smod ~ exp(-linear.kkpaft.prediction))[1])


###### penCox ###################################################################

linear.kkpcox.prediction <- c(0)
for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.aft <- cv.glmnet(x = Y[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox")
  linear.kkpcox.prediction[ind.new] <- predict(object =reg.aft, newx = Y.new[ind.new,], s= "lambda.min")
}
predCIndex.kkpcox <<- as.numeric(survConcordance(smod ~ linear.kkpcox.prediction)[1])


#### PenFLXMIX  to make predictions on Clusters ##################################################################
data <- data.frame(time, Y)
data.new <- data.frame(time.new, Y.new)
## The cross validation folds for choosing lambda
fo <- sample(rep(seq(10), length = nrow(data)))
gr.flx <- flexmix(time ~ ., data = data, k = F, cluster = gr.km$cluster, model = FLXMRglmnet(foldid = fo, adaptive= FALSE), control = list(iter.max = 500))
cl.flx <- clusters(gr.flx, newdata = data.new)
rand.flx.prediction <<-  adjustedRandIndex(c.true.new,as.factor(cl.flx))



# ####### PREMIUM ############################
# direc <- as.character('/home/bit/ashar/ownCloud/Research/DPMMSIMULATIONS/OneView/D20Noise20perOverlap01per/premium')
# setwd(dir)
# data.premium <- list(NULL)
# data.premium$inputData <- as.data.frame(cbind(Y,time))
# data.premium$covNames <- paste0(rep("Variable",D),as.character(c(1:D)))
# names(data.premium$inputData) <- c(data.premium$covNames,"outcome")
# data.premium$xModel <- "Normal"
# data.premium$yModel <- "Normal"
# data.premium$ncovariates <- D
# data.premium$outputData <- as.data.frame(Y.new)
# names(data.premium$outputData) <- paste0(rep("Variable",D),as.character(c(1:D)))
# 
# 
# ####### Run the Profile Regression
# 
# runInfoObj <-profRegr(yModel=data.premium$yModel, xModel= data.premium$xModel, nSweeps= 200, nBurn= 100, data= data.premium$inputData, output= "outcome", nFilter= 5 ,covNames= data.premium$covNames,nClusInit= k,reportBurnIn=FALSE, predict = data.premium$outputData)
# dissimObj<- calcDissimilarityMatrix(runInfoObj)
# clusObj<-calcOptimalClustering(dissimObj)
# 
# 
# ######## Rand Indices ##############################
# rand.recovery.premium <<- adjustedRandIndex(c.true,clusObj$clustering)
# rand.predicted.premium <<-  adjustedRandIndex(c.true.new,clusObj$clusteringPred)
# 
# ########## Making predictions #####################
# riskProfileObj <- calcAvgRiskAndProfile(clusObj)
# 
# output_predictions <- calcPredictions(riskProfileObj)
# 
# cindex.predicted.premium <<-  as.numeric(survConcordance(Surv(exp(time.new),censoring.new) ~ exp(-output_predictions$predictedY))[1]) 

  
  
  
}