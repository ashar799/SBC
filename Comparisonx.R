############ Ground Truth on TRAINING DATA ###################################
############## PREDICTIONS ON TESTING DATA  ##################################

### K-means + CoxPH
### K-means + AFT

### K-means + Penalized CoxPH
### K-means + Penalized AFT

### Mixture of Factor Analyzers
#### Sparse KMeans
#### Sparse Hierarchical clustering 


Comparisonx = function(){
  
  smod <-  Surv(exp(time), censoring)
  smod.new <-  Surv(exp(time.new), censoring.new)

  ### Fitting A Penalized Cox Proportional Hazard's Model
  reg.pcox <- cv.glmnet(x = Y, y = smod, family = "cox")
  lp <- predict(object =reg.pcox, newx = Y, s= "lambda.min")
  recovCIndex.na.cox <<- as.numeric(survConcordance(smod ~lp)[1])
  
  ######## Prediction with penalized Cox Proportional Hazards Model ###########################################
  linear.pred.cox <- c(0)
  ### see if we can use glmnet
  reg.pcox <- cv.glmnet(x = Y, y = Surv(exp(time), censoring), family = "cox")
  linear.pred.cox <- predict(object =reg.pcox, newx = Y.new, s= "lambda.min")
  predCIndex.na.cox <<- as.numeric(survConcordance(smod.new ~ linear.pred.cox)[1])
  
  
  #### Fitting a penalized AFT Model ####
  reg <- cv.glmnet(x = Y, y = time, family = "gaussian")
  linear.aft <- predict(object = reg, newx = Y, s = "lambda.min") 
  cindex.paft <- as.numeric(survConcordance(smod ~ exp(-linear.aft))[1])
  recovCIndex.na.aft <<- as.numeric(cindex.paft)
  ##### Prediction using the penlized AFT #############################################
  linear.pred.paft <- c(0)
  ### see if we can use glmnet
  reg.paft <- cv.glmnet(x = Y, y = time, family = "gaussian")
  linear.pred.paft <- predict(object = reg.paft, newx = Y.new, s= "lambda.min")
  predCIndex.na.aft <<- as.numeric(survConcordance(smod.new ~ exp(-linear.pred.paft))[1])
  
  #############################################
  ########### K-means #########################
  ############ K-Nearest Neighbour ############
  #############################################
  
  gr.km <- kmeans(Y, F, nstart =10)
  recovRandIndex.km <<- adjustedRandIndex(c.true,as.factor(gr.km$cluster))
  label.train <- gr.km$cluster
  label.test <- knn(train = Y, test = Y.new, cl = label.train, k = F)
  predRandIndex.knear <<- adjustedRandIndex(c.true.new, label.test)
  
  
  ###### penCox ###################################################################
  ######## Penalized Cox PH with k-means clustering###########################################
  linear.cox <- c(0)
  for ( q in 1:F){
    ind <- which((gr.km$cluster) == q)
    time.tmp <- time[ind]
    censoring.tmp <- censoring[ind]
    Y.tmp <- Y[ind,]
    coxreg <- list(0)
    coxreg$x <- Y.tmp
    coxreg$time <- exp(time.tmp)
    coxreg$status <- censoring.tmp
    reg.pcox <- cv.glmnet(x = Y.tmp, y = Surv(coxreg$time, coxreg$status), family = "cox")
    linear.cox[ind] <- predict(object =reg.pcox, newx = Y.tmp, s= "lambda.min")
  }
  
  recovCIndex.km.pcox <<- as.numeric(survConcordance(smod ~ linear.cox)[1])
  ### Prediction with k-means + k-nearest neghbour
  linear.kkpcox.prediction <- c(0)
  for ( q in 1:F){
    ind <- which(label.train == q)
    ind.new <- which(label.test == q)
    reg.aft <- cv.glmnet(x = Y[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox")
    linear.kkpcox.prediction[ind.new] <- predict(object =reg.aft, newx = Y.new[ind.new,], s= "lambda.min")
  }
  predCIndex.kk.pcox <<- as.numeric(survConcordance(smod.new ~ linear.kkpcox.prediction)[1])
  
  
  
  ###### penAFT ###################################################################
  ######## Penalized AFT with k-means clustering ######################################################
  linear.aft <- c(0)
  for ( q in 1:F){
    ind <- which((gr.km$cluster) == q)
    L= length(ind)
    
    time.tmp <- time[ind]
    censoring.tmp <- censoring[ind]
    Y.tmp <- Y[ind,]
    
    reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
    coeff.pred <- coef(object =reg, newx = Y.tmp, s= "lambda.min")
    rel.coeff <- coeff.pred[2:(D+1)] 
    ind.rel <- which(rel.coeff !=0)
    linear.aft[ind] <- predict(object = reg, newx = Y.tmp, s = "lambda.min") 
  }
  recovCIndex.km.paft <<- as.numeric(survConcordance(smod ~ exp(-linear.aft))[1])
  
  
  
  #### prediction with penAFT ###################################################################
  linear.kkpaft.prediction <- c(0)
  for ( q in 1:F){
    ind <- which(label.train == q)
    ind.new <- which(label.test == q)
    reg.aft <- cv.glmnet(x = Y[ind,], y = time[ind], family = "gaussian")
    linear.kkpaft.prediction[ind.new] <- predict(object =reg.aft, newx = Y.new[ind.new,], s= "lambda.min")
  }
  predCIndex.kn.paft <<- as.numeric(survConcordance(smod.new ~ exp(-linear.kkpaft.prediction))[1])
  

  
  
  
  ######## Model fitting with Penalized Cox PH with TRue Clustering ###########################################
  ######## Model prediction with k-nn based on true clustering #####

  true.knn <- knn(train = Y, test = Y.new, cl = c.true, k = F)
  predRandIndex.true.knear <<- adjustedRandIndex(c.true.new, true.knn)
  
  
  linear.pred <- c(0)
  cox.pred <- c(0)
  for ( q in 1:F){
    ind <- which((c.true) == q)
    ind.new <- which(true.knn == q)
    time.tmp <- time[ind]
    censoring.tmp <- censoring[ind]
    Y.tmp <- Y[ind,1:D]
    Y.tmp.new <- Y.new[ind.new,1:D]
    coxreg <- list(0)
    coxreg$x <- Y.tmp
    coxreg$time <- exp(time.tmp)
    coxreg$status <- censoring.tmp
    reg.pcox <- cv.glmnet(x = Y.tmp, y = Surv(coxreg$time, coxreg$status), family = "cox")
    linear.pred[ind] <- predict(object =reg.pcox, newx = Y.tmp, s= "lambda.min")
    cox.pred[ind.new] <- predict(object =reg.pcox, newx = Y.tmp.new, s= "lambda.min")
    }

  recovCIndex.true.pcox <<- as.numeric(survConcordance(smod ~ linear.pred)[1])
  predCIndex.true.knn.pcox <<- as.numeric(survConcordance(smod.new ~ cox.pred)[1])
  
  
  ######## Penalized AFT with TRUE clustering ######################################################
  linear.aft <- c(0)
  pred.aft <- c(0)
  
  
  for ( q in 1:F){
    ind <- which((c.true) == q)
    ind.new <- which(true.knn == q)
    L= length(ind)
    
    time.tmp <- time[ind]
    censoring.tmp <- censoring[ind]
    Y.tmp <- Y[ind,]
    Y.tmp.new <- Y.new[ind.new,1:D]
    
    
    reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
    linear.pred <- predict(object =reg, newx = Y.tmp, s= "lambda.min")
    coeff.pred <- coef(object =reg, newx = Y.tmp, s= "lambda.min")
    rel.coeff <- coeff.pred[2:(D+1)] 
    ind.rel <- which(rel.coeff !=0)
    linear.aft[ind] <- predict(object = reg, newx = Y.tmp, s = "lambda.min") 
    pred.aft[ind.new] <- predict(object = reg, newx = Y.tmp.new, s = "lambda.min")
  }
  recovCIndex.true.paft <<- as.numeric(survConcordance(smod ~ exp(-linear.aft))[1])
  predCIndex.true.knn.paft <<- as.numeric(survConcordance(smod.new ~ exp(-pred.aft))[1])
  
  
  
  

  ############## Testing Mixture of Factor Analyzers ##########################3
  # mcfa.fit<- mcfa(Y, g= k, q=2, itmax=250, nkmeans=5, nrandom=5, tol=1.e-3)
  
  ########## Seeing if the PCA plot with the corresponding features with releevant features makes sense
  #randindexMCFA <<- adjustedRandIndex(mcfa.fit$clust, c.true)
  
  
  #############################################
  ########### sparse K-means #########################
  #############################################
  #############################################
  km.perm <- KMeansSparseCluster.permute(x = Y, K= k ,wbounds=c(1.5,2:6),nperms=5)
  km.out <- KMeansSparseCluster(x = Y, K=k,wbounds=km.perm$bestw)
  recovRandIndexSKM <<- adjustedRandIndex(km.out[[1]]$Cs, c.true)
  
  ###################################################
  ########### sparse hierarchical clustering #########################
  #############################################
  #############################################
  perm.out <- HierarchicalSparseCluster.permute(x = Y, wbounds=c(1.5,2:6), nperms=5)
  # Perform sparse hierarchical clustering
  sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists, wbound=perm.out$bestw, method="complete")
  recovRandIndexSHC <<- adjustedRandIndex(cutree(sparsehc$hc, k = k), c.true)
  
  
  
  ### Let's see if we can verify if the model likelihood for the real c is indeed the highest
  ### Use mixture of regression model using mixtools
  # library(mixtools)
  # out <- regmixEM.loc(y = time, x = Y[,1:D], lambda = NULL, beta = NULL, sigma = NULL,
  #              k = 2, addintercept = TRUE, kern.l = c("Gaussian"),
  #              epsilon = 1e-08, maxit = 10000, kernl.g = 0,
  #              kernl.h = 1, verb = TRUE)
  # 
  
  
  ############ Predicting New Class Labels using SVM #################################  
  # gr.km <- kmeans(Y, F, nstart =10)
  # label.train <- gr.km$cluster
  # svms <- sapply(2^(-10:14), function(cost) cross(ksvm(Y, factor(label.train), C=cost, kernel="vanilladot", kpar=list(), cross=5)))
  # mysvm <- ksvm(Y, factor(label.train), C=2^(-10:14)[which.min(svms)], kernel="vanilladot", kpar=list(), cross=10) # accuracy ~97%
  # pred.svm <-  predict(mysvm, Y.new)
  # predRandIndex.svm <<- adjustedRandIndex(c.true.new, pred.svm)
  
  
  
  

  
  
 
  
}
