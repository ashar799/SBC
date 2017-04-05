### This function compares feature selection and C-Index recovery and prediction ###
### Let's see what comes up



ComparisionFLX = function(){
  
 
  gr.km <- kmeans(Y, F, nstart =10)
  gr.km.rand <- adjustedRandIndex(c.true,as.factor(gr.km$cluster))
  
  
  #################Flexmix Clustering ################################################################
  ##############################################################################
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
  
  linear.flx <- c(0)
  beta.flx <- matrix(0, nrow = D, ncol = F)
  for ( q in 1:F){
    ind <- which(clusters(gr.flx) == q)
    L= length(ind)
    
    time.tmp <- time[ind]
    censoring.tmp <- censoring[ind]
    Y.tmp <- Y[ind,]
    
    reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
    
    coeff.pred <- coef(object =reg, newx = Y.tmp, s= "lambda.min")
    rel.coeff <- coeff.pred[2:(D+1)] 
    beta.flx[1:D,q] <- ((rel.coeff != 0)+0)
    
    linear.flx[ind] <- predict(object = reg, newx = Y.tmp, s = "lambda.min") 
  }
  recovCIndex.flx.aft  <<- as.numeric(survConcordance(smod ~ exp(-linear.flx))[1])
  
  
  ### Use flexmix clustering
  data.new <- data.frame(x = Y.new)
  linear.pred.flx <- as.numeric(unlist(predict(gr.flx, newdata = data.new, aggregate = TRUE)))
  predCIndex.flx <<- as.numeric(survConcordance(smod.new ~ exp(-linear.pred.flx))[1])
  
  
  ### Use FLexmix to make predictions on future probabilities ####
  # cl.new <- posterior(gr.flx, newdata = data.new)
  # cl.pred.flx <- c(0)
  # 
  # for ( i in 1: dim(Y.new)[1]){
  #   cl.pred.flx[i] <- which.max(cl.new[i,])
  #   
  # }
  # predRandIndex.flx <<- adjustedRandIndex(c.true.new,as.factor(cl.pred.flx))
  # 
  # 
  # ### Use PAFT now to build the corresponding cluster specific PAFT models
  # pre.flx <- c(0)
  # 
  # for ( q in 1:F){
  #   ind <- which(clusters(gr.flx) == q)
  #   ind.new <- which(cl.pred.flx == q)
  #   
  #   time.tmp <- time[ind]
  #   
  #   Y.tmp <- Y[ind,]
  #   Y.tmp.new <- Y.new[ind.new,]
  #   
  #   reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
  #   pre.flx[ind.new] <- predict(object = reg, newx = Y.tmp.new, s = "lambda.min") 
  # }
  # predCIndex.flx.aft  <<- as.numeric(survConcordance(smod.new ~ exp(-pre.flx))[1])
  # 
  # 
  
  
  
  
  
  
  
}
