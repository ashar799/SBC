##### This file concatenates the previous files of TRAIN and Test multi Comparison.R
##### The comparison also helps to initialize the model
##### Trying something more innovative for initializing the iSBC model 

### K-means + Penalized CoxPH
### K-means + Penalized AFT

### FlexMix +  CoxPH
### FlexMix +  AFT


### iCLUSTER +AFT


multiComparisonx = function(){
  
  Y <- cbind(Y1,Y2)
  Y.new <- cbind(Y1.test,Y2.test)
  D <- D1 + D2
  
  smod <-  Surv(exp(time), censoring)
  smod.new <- Surv(exp(time.new), censoring.new)
  
  ############ No CLUSTERING INFORMATION ############################################
  ### Fitting A Penalized Cox Proportional Hazard's Model
  reg.pcox <- cv.glmnet(x = Y, y = smod, family = "cox")
  lp <- predict(object =reg.pcox, newx = Y, s= "lambda.min")
  recovCIndex.na.pcox <<- as.numeric(survConcordance(smod ~lp)[1])
  linear.pred.cox <- predict(object =reg.pcox, newx = Y.new, s= "lambda.min")
  predCIndex.na.pcox <<- as.numeric(survConcordance(smod.new ~ linear.pred.cox)[1])
  
  
  #### Fitting A AFT Model #####
  reg.paft <- cv.glmnet(x = Y, y = time, family = "gaussian")
  linear.aft <- predict(object = reg.paft, newx = Y, s = "lambda.min") 
  recovCIndex.na.paft <<- as.numeric(survConcordance(smod ~ exp(-linear.aft))[1])
  linear.pred.paft <- predict(object = reg.paft, newx = Y.new, s= "lambda.min")
  predCIndex.na.paft <<- as.numeric(survConcordance(smod.new ~ exp(-linear.pred.paft))[1])
  
  #############################################
  ########### K-means #########################
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
  
  
  
  
  ##############################################################################
  ############### FlexMix #######################################################
  data <- data.frame(y =time, x =  Y)
  ## The cross validation folds for choosing lambda
  fo <- sample(rep(seq(10), length = nrow(data)))
  gr.flx <- flexmix(y ~ ., data = data, k = F, cluster = gr.km$cluster, model = FLXMRglmnet(foldid = fo, adaptive= FALSE, family = c("gaussian")), control = list(iter.max = 500))
  recovRandIndex.flx <<-  as.numeric(adjustedRandIndex(c.true,as.factor(clusters(gr.flx))))
  linear.recov.flx <- as.numeric(unlist(predict(gr.flx, newdata = data, aggregate = TRUE)))
  recovCIndex.flx <<- as.numeric(survConcordance(smod ~ exp(-linear.recov.flx))[1])
  
  
  ############### Penalized FlexMix #######################################################
  ################################################################################
  ## Fit a AFT model with FLXmix clustering
  
  # linear.flx <- c(0)
  # beta.flx <- matrix(0, nrow = D, ncol = F)
  # for ( q in 1:F){
  #   ind <- which(clusters(gr.flx) == q)
  #   L= length(ind)
  #   
  #   time.tmp <- time[ind]
  #   censoring.tmp <- censoring[ind]
  #   Y.tmp <- Y[ind,]
  #   
  #   reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
  #   
  #   coeff.pred <- coef(object =reg, newx = Y.tmp, s= "lambda.min")
  #   rel.coeff <- coeff.pred[2:(D+1)] 
  #   beta.flx[1:D,q] <- ((rel.coeff != 0)+0)
  #   
  #   linear.flx[ind] <- predict(object = reg, newx = Y.tmp, s = "lambda.min") 
  # }
  # recovCIndex.flx.aft  <<- as.numeric(survConcordance(smod ~ exp(-linear.flx))[1])
  # 
  # 
  
  ### Prediction ####
  smod.new <-  Surv(exp(time.new), censoring.new)
  
  ### Use flexmix clustering
  data.new <- data.frame(x = Y.new)
  linear.pred.flx <- as.numeric(unlist(predict(gr.flx, newdata = data.new, aggregate = TRUE)))
  predCIndex.flx <<- as.numeric(survConcordance(smod.new ~ exp(-linear.pred.flx))[1])
  
  
  
   ############## Using iCluster #######
  cv.fit <- tune.iClusterPlus(cpus=3,dt1 = Y1, dt2= Y2,
                              type=c("gaussian","gaussian"), K=F-1 ,n.lambda= 21,scale.lambda=c(1,1),
                              n.burnin=200,n.draw=200,maxiter=20,sdev=0.05,eps=1.0e-4)
  
  
  ### Now choosing that clustering that gives good separability and good c-index recovery and prediction
  ### Predicting the class labels is based on knn
  
  
  recov.CIndex.icl <- c(0)
  pred.CIndex.icl <- c(0)
  
  
 for ( i in 1:21){ 
  
  label.train <- cv.fit[[1]][[i]]$clusters
  label.test <- knn(train = Y, test = Y.new, cl = label.train, k = F)
  
  
  
  recov.icl <- c(0)
  pred.icl <- c(0)
  
  
  for ( q in 1:F){
    ind <- which(label.train == q)
    ind.new <- which(label.test == q)
    reg.aft <- cv.glmnet(x = Y[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox")
    recov.icl[ind] <- predict(object =reg.aft, newx = Y[ind,], s= "lambda.min")
    pred.icl[ind.new] <- predict(object =reg.aft, newx = Y.new[ind.new,], s= "lambda.min")
    
  }
  recov.CIndex.icl[i]  <- as.numeric(survConcordance(smod ~ recov.icl)[1])
  pred.CIndex.icl[i]  <- as.numeric(survConcordance(smod.new ~ pred.icl)[1])
 }

    
  recovCIndex.icl.pcox <<- recov.CIndex.icl
  predCIndex.icl.pcox <<- pred.CIndex.icl
  
}
