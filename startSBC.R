## This function was proposed on 24 March 2017 (3 years after the project) to initialize the model to a better starting point
## That initial clustering is chosen which recoverrs the origninal model well
### This function sets the very best values for the parameters right at the start of the Burnin period
### The idea is to Have the GMM parameters be influenced by the survival time at the very beginning


startSBC = function(){
  
  
  ################################# USE k-Means and Then FlxMIX  ###################################################
  ##### See which gives a better c-index recovery #######
  
  ## K-means clustering
  gr.km <- kmeans(Y, F, nstart =10)
  
 
  ## FlexMix Clustering
  data <- data.frame(time, Y)
  fo <- sample(rep(seq(10), length = nrow(data)))
  gr.flx <- flexmix(time ~ ., data = data, k = F, cluster = gr.km$cluster, model = FLXMRglmnet(foldid = fo, adaptive= FALSE, family = c("gaussian")), control = list(iter.max = 500))
  
  
  ## Fit a cluster specific AFT model with kmeans clustering
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
  recovCIndex.km.paft <- as.numeric(survConcordance(smod ~ exp(-linear.aft))[1])
  
  
  
  ## Fit a cluster specific AFT model with FLXmix clustering
  
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
    beta.flx[1:D,q] <- rel.coeff 
    
    linear.flx[ind] <- predict(object = reg, newx = Y.tmp, s = "lambda.min") 
  }
  recovCIndex.flx.paft  <- as.numeric(survConcordance(smod ~ exp(-linear.flx))[1])
  
  
  ### Change the cluster membership if the FLXMix clustering is better
  
  if(recovCIndex.flx.paft >  recovCIndex.km.paft){
    
    c <- clusters(gr.flx)
    
    prior.numclust <- table(factor(c, levels = 1:K))
    prior.activeclass<- which(prior.numclust!=0)
    
    ### The means are set using the k-means
    for ( i in 1:length(prior.activeclass)){
      
      ind <- which(c == prior.activeclass[i])
      
      Y.tmp <- Y[ind,]
      
      time.tmp <- time[ind]
      
      
      mu[prior.activeclass[i],1:D] <-  apply(Y.tmp,2,mean) 
      
      S[prior.activeclass[i],1:D,1:D] <-  solve(cov(Y.tmp) + diag(1,D))
      
      
      beta0[prior.activeclass[i]] <- mean(time.tmp)
      
      
      
      betahat[prior.activeclass[i], 1:D] <- beta.flx[1:D,prior.activeclass[i]]
      
      
      ##### For the hyper-parameters use the BLASSO function ####
      ##### This too disregards censoring information ###########
  
      
      reg.blas <- 0
      
      sum <- c(0)
      
      coeff <- 0
      
      Ytemp <-  matrix(NA, nrow = length(ind), ncol = D)
      
      Ytemp <- scale(Y[ind,1:D], center = TRUE, scale = TRUE)
      
      Ttemp <- as.vector(time[ind])
      
      ntemp <- length(ind)
      
      reg.blas <- blasso(Ytemp, Ttemp, T = 300,thin =  50, RJ = TRUE, mprior = 0.0 ,normalize = TRUE, verb = 0)
      
      sum <- summary(reg.blas, burnin= 100)
      
      ## Selecting those features which are relevant
      
      coeff <- unlist(lapply(strsplit(sum$coef[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
      
      ta <- unlist(lapply(strsplit(sum$tau2i[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
      
      ta.impute <- impute(ta)
      
      tau2[prior.activeclass[i],] <- ta.impute
      
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
    
  } 
    
    
    
  assign("c", c, envir = .GlobalEnv)
  assign("mu", mu, envir = .GlobalEnv)
  assign("S", S, envir = .GlobalEnv)
  assign("beta0", beta0, envir = .GlobalEnv)
  assign("betahat", betahat, envir = .GlobalEnv)
  assign("sigma2", sigma2, envir = .GlobalEnv)
  assign("lambda2", lambda2, envir = .GlobalEnv)
  assign("tau2", tau2, envir = .GlobalEnv)
  assign("lambda2", lambda2, envir = .GlobalEnv)
  
    
    
    
    
    
  }
  
 
  
  
  
  

