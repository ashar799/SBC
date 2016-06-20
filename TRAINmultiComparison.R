############ Ground Truth on TRAINING DATA ###################################
##############################################################
###########

### K-means + Penalized CoxPH
### K-means + Penalized AFT

### FlexMix +  CoxPH
### FlexMix +  AFT


### iCLUSTER
### k CCA 
### kmean sparse CCA


multigroundtruth = function(){
  
  Y <- cbind(Y1,Y2)
  D <- D1 + D2
  
  smod <-  Surv(exp(time), censoring)
  
  ############ No CLUSTERING INFORMATION ############################################
  ##### Both Data Sets put together
  ### Fitting A Penalized Cox Proportional Hazard's Model
  reg.pcox <- cv.glmnet(x = Y, y = smod, family = "cox")
  lp <- predict(object =reg.pcox, newx = Y, s= "lambda.min")
  cindex.pcox <- survConcordance(smod ~lp)[1]
  cindex.pen.cox <<- as.numeric(cindex.pcox)
  
  
  #### Fitting A AFT Model 
  reg <- cv.glmnet(x = Y, y = time, family = "gaussian")
  linear.aft <- predict(object = reg, newx = Y, s = "lambda.min") 
  cindex.paft <- survConcordance(smod ~ exp(-linear.aft))[1]
  cindex.pen.aft <<- as.numeric(cindex.paft)
  
  
  #############################################
  ########### K-means #########################
  #############################################
  #############################################
  
  
  gr.km <- kmeans(Y, F, nstart =10)
  gr.km.rand <- adjustedRandIndex(c.true,as.factor(gr.km$cluster))
  
  
  
  
  
  
  ################# Combined Data Set ########################
  
  cindex.km.pcox <-0
  cindex.km.paft <- 0
  ######## Penalized Cox PH ###########################################
  linear.pred <- c(0)
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
    linear.pred[ind] <- predict(object =reg.pcox, newx = Y.tmp, s= "lambda.min")
  }
  
  cindex.km.pcox <- survConcordance(smod ~ linear.pred)[1]
  
  
  
  
  ######## Penalized AFT ######################################################
  
  linear.aft <- c(0)
  for ( q in 1:F){
    ind <- which((gr.km$cluster) == q)
    L= length(ind)
    
    time.tmp <- time[ind]
    censoring.tmp <- censoring[ind]
    Y.tmp <- Y[ind,]
    
    reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
    linear.pred <- predict(object =reg, newx = Y.tmp, s= "lambda.min")
    coeff.pred <- coef(object =reg, newx = Y.tmp, s= "lambda.min")
    rel.coeff <- coeff.pred[2:(D+1)] 
    ind.rel <- which(rel.coeff !=0)
    linear.aft[ind] <- predict(object = reg, newx = Y.tmp, s = "lambda.min") 
  }
  cindex.km.paft <- survConcordance(smod ~ exp(-linear.aft))[1]
  
  
  
  ##### Save some Ground truth statistics
  
  gr.km.rand.final <<- gr.km.rand
  cindex.km.pcox.final <<- as.numeric(cindex.km.pcox)
  cindex.km.paft.final <<- as.numeric(cindex.km.paft)
  
  
  
  #################################################################################
  ##############################################################################
  ############### FlexMix #######################################################
  ################################################################################
  
  gr.flx <- flexmix(time ~ Y, k =F)
  gr.flx.rand <- adjustedRandIndex(c.true,clusters(gr.flx))
  
  
  
  ########## CoxPH #############################
  fit.cox.flx <- coxph(smod ~ Y[,1:D]*strata(as.factor(clusters(gr.flx))), data = as.data.frame(Y))
  ## C-Index
  cindex.flx.cox <- survConcordance(smod ~ predict(fit.cox.flx))[1]
  ## Brier Score
  fit.coxph <- survfit(fit.cox.flx, newdata = as.data.frame(Y[,1:D]))
  gr.flx.rand.final <<- gr.flx.rand
  cindex.flx.cox.final <<- as.numeric(cindex.flx.cox)
  
  
  ############## Using iCluster #######
  datas <- list(0)
  datas[[1]] <- Y1
  datas[[2]] <- Y2
  cv.fit <- tune.iCluster2(datas, k)
  fit <- iCluster2(datas, k= k, lambda= cv.fit$best.fit$lambda)
  randindexiCLUSTER <<- adjustedRandIndex(fit$clusters,c.true)
  
  ########### Using CCA ################
  fit.cc <- cc(Y1, Y2)
  y1 <- fit.cc$scores$xscores
  y2 <- fit.cc$scores$yscores
  f <- length(which(fit.cc$cor > 0.5))
  Y.CCA <- cbind(y1[,1:f],y2[,1:f])
  km <- kmeans(Y.CCA, centers =k, nstart =10)
  randindexCCA <<- adjustedRandIndex(c.true,as.factor(km$cluster))
  
  
  
  
}
