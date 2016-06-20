########## This file compares prediction on test cases with different methods
####### So we use k-means + SVM to predict labels on new data points ###########
####### Then fit penalized Cox- Models with the predicted Labels


predictionGroundTruth = function(){
  
  
  ############ Predicting New Class Labels #################################  
  Y <- cbind(Y1,Y2)
  Y.new <- cbind(Y1.test,Y2.test)
  
  gr.km <- kmeans(Y, F, nstart =10)
  label.train <- gr.km$cluster
  svms <- sapply(2^(-10:14), function(cost) cross(ksvm(Y, factor(label.train), C=cost, kernel="vanilladot", kpar=list(), cross=5)))
  mysvm <- ksvm(Y, factor(label.train), C=2^(-10:14)[which.min(svms)], kernel="vanilladot", kpar=list(), cross=10) # accuracy ~97%
  pred.svm <-  predict(mysvm, Y.new)
  predRandIndex.svm <<- adjustedRandIndex(c.true.new, pred.svm)
  
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
  
  
  
  
  
  
  
  
  
}