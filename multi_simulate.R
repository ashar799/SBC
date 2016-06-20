######################################### This function simulates Multi (2) views of DP mixture model #####################
############################################################################################################################

multi_simulate = function(){
  
  
  
  #################################################### FIRST DATA SET #################################
  #####################################################################################################
  ####################################################################################################
  
  ## Actual Number of Components and dimension  which are relevant
  rel.D1  <- as.integer(D1* (1-prob.noise.feature))
  ## Actual Number of Irrelevant Componenets
  irrel.D1 <-  D1 - rel.D1
  
  #################################################################################################
  ###############################################################################################
  ########## Molecular Data Generatioin ###########################################################
  ##################################################################################################
  #################################################################################################
  
  
  ## The data
  A <- MixSim(BarOmega = prob.overlap ,K = F, p = rel.D1, int =c(-1.0,1.0), lim = 1e09)
  data.mu = array(data = NA, dim =c(F,rel.D1))
  data.S = array(data = NA, dim =c(F,rel.D1,rel.D1))
  
  
  # Data coming from a hypothetical population NON diagoanl matrix
  for( i in 1:F){
    data.mu[i,1:rel.D1] <- A$Mu[i,1:rel.D1]
    data.S[i,1:rel.D1,1:rel.D1] <- A$S[1:rel.D1,1:rel.D1,i]
  }
  
  Y1.rel.list <- list(0)
  for ( i in 1:F){
    Y1.rel.list[[i]] <- mvrnorm(n = as.integer(N.train * p.dist[i]), mu = data.mu[i,1:rel.D1], Sigma = data.S[i,1:rel.D1,1:rel.D1])
  }
  
  Y1.rel.list.test <- list(0)
  for ( i in 1:F){
    Y1.rel.list.test[[i]] <- mvrnorm(n = as.integer(N.test * p.dist[i]), mu = data.mu[i,1:rel.D1], Sigma = data.S[i,1:rel.D1,1:rel.D1])
  }
  
  ## Scaling the Data as ONLY the scaled data will be used for generating the times
  Y1.rel.sc.list <- list(0)
  Y1.rel.sc.list.test <- list(0)
  for ( i in 1:F){
    obj <- scale(Y1.rel.list[[i]], center = TRUE, scale = TRUE)
    Y1.rel.sc.list[[i]] <- scale(Y1.rel.list[[i]], center = attr(obj,"scaled:center"), scale = attr(obj,"scaled:scale"))
    Y1.rel.sc.list.test[[i]] <- scale(Y1.rel.list.test[[i]], center = attr(obj,"scaled:center"), scale = attr(obj,"scaled:scale"))
  }
  
 ## Irrelevant features
  Y1.irrel.list <- list(0)
  for ( i in 1:F){
    mean <- runif(irrel.D1,-1.5,1.5)
    Y1.irrel.list[[i]] <- mvrnorm(n = as.integer(N.train * p.dist[i]), mu = mean, Sigma = diag(x =1, nrow = irrel.D1, ncol = irrel.D1))
  }
  
  Y1.irrel.list.test <- list(0)
  for ( i in 1:F){
    mean <- runif(irrel.D1,-1.5,1.5)
    Y1.irrel.list.test[[i]] <- mvrnorm(n = as.integer(N.test * p.dist[i]), mu = mean, Sigma = diag(x =1, nrow = irrel.D1, ncol = irrel.D1))
  }
  
  
  ### Combining the data with relevant and irrelevant columns
  data1.old <- list(0) 
  for (i in 1:F){
    data1.old[[i]] <-  cbind(Y1.rel.list[[i]], Y1.irrel.list[[i]]) 
  }
  
  data1.old.test <- list(0) 
  for (i in 1:F){
    data1.old.test[[i]] <-  cbind(Y1.rel.list.test[[i]], Y1.irrel.list.test[[i]]) 
  }
  
  
  ################### SECOND DATA SET #################################################################
  #####################################################################################################
  ## Actual Number of Components and dimension  which are relevant
  rel.D2  <- as.integer(D2* (1-prob.noise.feature))
  ## Actual Number of Irrelevant Componenets
  irrel.D2 <-  D2 - rel.D2
  
  #################################################################################################
  ###############################################################################################
  ########## Molecular Data Generatioin ###########################################################
  ##################################################################################################
  #################################################################################################
  
  
  ## The data
  A <- MixSim(BarOmega = prob.overlap ,K = F, p = rel.D2, int =c(-1.0,1.0), lim = 1e09)
  data.mu = array(data = NA, dim =c(F,rel.D2))
  data.S = array(data = NA, dim =c(F,rel.D2,rel.D2))
  
 
  # Data coming from a hypothetical population NON diagoanl matrix
  for( i in 1:F){
    data.mu[i,1:rel.D2] <- A$Mu[i,1:rel.D2]
    data.S[i,1:rel.D2,1:rel.D2] <- A$S[1:rel.D2,1:rel.D2,i]
  }
  
  Y2.rel.list <- list(0)
  for ( i in 1:F){
    Y2.rel.list[[i]] <- mvrnorm(n = as.integer(N.train * p.dist[i]), mu = data.mu[i,1:rel.D2], Sigma = data.S[i,1:rel.D2,1:rel.D2])
  }
  
  Y2.rel.list.test <- list(0)
  for ( i in 1:F){
    Y2.rel.list.test[[i]] <- mvrnorm(n = as.integer(N.test * p.dist[i]), mu = data.mu[i,1:rel.D2], Sigma = data.S[i,1:rel.D2,1:rel.D2])
  }
  
  
 ## Scaling the Data as ONLY the scaled data will be used for generating the times
 Y2.rel.sc.list <- list(0)
 Y2.rel.sc.list.test <- list(0)
 for ( i in 1:F){
   obj <- scale(Y2.rel.list[[i]], center = TRUE, scale = TRUE)
   Y2.rel.sc.list[[i]] <- scale(Y2.rel.list[[i]], center = attr(obj,"scaled:center"), scale = attr(obj,"scaled:scale"))
   Y2.rel.sc.list.test[[i]] <- scale(Y2.rel.list.test[[i]], center = attr(obj,"scaled:center"), scale = attr(obj,"scaled:scale"))
 }
  
  
  ## Irrelevant features
  Y2.irrel.list <- list(0)
  for ( i in 1:F){
    mean <- runif(irrel.D2,-1.5,1.5)
    Y2.irrel.list[[i]] <- mvrnorm(n = as.integer(N.train * p.dist[i]), mu = mean, Sigma = diag(x =1, nrow = irrel.D2, ncol = irrel.D2))
  }
  
  Y2.irrel.list.test <- list(0)
  for ( i in 1:F){
    mean <- runif(irrel.D2,-1.5,1.5)
    Y2.irrel.list.test[[i]] <- mvrnorm(n = as.integer(N.test * p.dist[i]), mu = mean, Sigma = diag(x =1, nrow = irrel.D2, ncol = irrel.D2))
  }
  
  
  ### Combining the data with relevant and irrelevant columns
  data2.old <- list(0) 
  for (i in 1:F){
    data2.old[[i]] <-  cbind(Y2.rel.list[[i]], Y2.irrel.list[[i]]) 
  }
  
  data2.old.test <- list(0) 
  for (i in 1:F){
    data2.old.test[[i]] <-  cbind(Y2.rel.list.test[[i]], Y2.irrel.list.test[[i]]) 
  }
  
 ######### Defining some lists
 data1.new <- list(0)
 data2.new <- list(0)
 data1.new.test <- list(0)
 data2.new.test <- list(0)
 
 
  
 ################## Dealing with indices #####################################################
 ###############################################################################################
 D <- D1+ D2
 rel.D <- rel.D1+ rel.D2
 irrel.D <- irrel.D1 +irrel.D2
 
 
 ind.data1 <- 1:rel.D1
 irrel.ind.data1 <- (rel.D1+ rel.D2 +1):(rel.D1+ rel.D2 +irrel.D1)
 ind.data2 <- (rel.D1+1):(rel.D1+ rel.D2)
 irrel.ind.data2 <- (max(c(ind.data1,ind.data2,irrel.ind.data1)) +1) :D
 
 ind.rel.D <- c(ind.data1, ind.data2)
 
 
  ######################################################################################### Making the irrelevant features independent from the dependent features #############
  ############ Training Data ###################################################################################################################################################
  ###############################################################################################################################################################################  
  for (f in 1:F){
    
    
  X <- cbind(data1.old[[f]][,1:rel.D1],data2.old[[f]][,1:rel.D2],data1.old[[f]][,(rel.D1+1):D1], data2.old[[f]][,(rel.D2+1):D2])
  
  rel.X <- as.matrix(X[,ind.rel.D])
  
  obj.qr <- qr(X)
  
  rk <- obj.qr$rank
  
  alpha <- qr.Q(obj.qr)[,1:rel.D]
  
  gamma <- qr.Q(obj.qr)[,(1+rel.D):rk]
  
  matT <- matrix(runif(n = rel.D*(rk -rel.D), min = -0.005, max= 0.005), nrow = rel.D, ncol = (rk -rel.D))
  
  matP <- t(matT) %*% matT
  
  max.eig <- eigen(matP)$values[1]
  
  max.corr <- sqrt(max.eig)/sqrt(1 + max.eig)
  
  linear.space <- gamma + alpha %*% matT
  
  irrel.X <- matrix(NA, nrow = nrow(X), ncol = irrel.D)
  
  for ( i in 1: irrel.D){
    
    matTemp <- matrix(runif(n = (rk -rel.D), min = -1.5, max= 1.5),  nrow = (rk-rel.D), ncol =1)
    irrel.X[,i] <- as.vector(linear.space %*% matTemp)
    
  }
  
  ## Checking if the covariance is indeed small
  
  cov.mat <- cov(rel.X,irrel.X)
  
  boxplot(cov.mat)
  
  ## Building the full data matrix
  
  X.full <- cbind(rel.X, irrel.X)
  
  
  levelplot(cov(X.full))
  
  data1.new[[f]] <- X.full[,c(ind.data1,irrel.ind.data1)]
  data2.new[[f]] <- X.full[,c(ind.data2,irrel.ind.data2)]
 
}
 
 ############################################### MAKING Y from the clusters data #####################3
 Y1 <- c(0)
 for (i in 1:F){
   Y1 <- rbind(Y1, data1.new[[i]])
 } 
 Y1 <- Y1[-1,]
 
 
 ############################################### MAKING Y from the clusters data #####################3
 Y2 <- c(0)
 for (i in 1:F){
   Y2 <- rbind(Y2 ,data2.new[[i]])
 } 
 Y2 <- Y2[-1,]
 
 
 
 
 
  ################################################################################################
  ############################### TESTING DATA ####################################################
  
for (f in 1:F){
 X <- cbind(data1.old.test[[f]][,1:rel.D1],data2.old.test[[f]][,1:rel.D2], data1.old.test[[f]][,(rel.D1+1):D1], data2.old.test[[f]][,(rel.D2+1):D2])
 
 rel.X <- as.matrix(X[,ind.rel.D])
 
 obj.qr <- qr(X)
 
 rk <- obj.qr$rank
 
 alpha <- qr.Q(obj.qr)[,1:rel.D]
 
 gamma <- qr.Q(obj.qr)[,(1+rel.D):rk]
 
 matT <- matrix(runif(n = rel.D*(rk -rel.D), min = -0.005, max= 0.005), nrow = rel.D, ncol = (rk -rel.D))
 
 matP <- t(matT) %*% matT
 
 max.eig <- eigen(matP)$values[1]
 
 max.corr <- sqrt(max.eig)/sqrt(1 + max.eig)
 
 linear.space <- gamma + alpha %*% matT
 
 irrel.X <- matrix(NA, nrow = nrow(X), ncol = irrel.D)
 
 for ( i in 1: irrel.D){
   
   matTemp <- matrix(runif(n = (rk -rel.D), min = -1.5, max= 1.5),  nrow = (rk-rel.D), ncol =1)
   irrel.X[,i] <- as.vector(linear.space %*% matTemp)
   
 }
 
 ## Checking if the covariance is indeed small
 
 cov.mat <- cov(rel.X,irrel.X)
 
 boxplot(cov.mat)
 
 ## Building the full data matrix
 
 X.full <- cbind(rel.X, irrel.X)
 
 
 levelplot(cov(X.full))
 
 data1.new.test[[f]] <- X.full[,c(ind.data1,irrel.ind.data1)]
 data2.new.test[[f]] <- X.full[,c(ind.data2,irrel.ind.data2)]
 
}

 
 ############################################### MAKING Y.test from the clusters data #####################3
 
Y1.test <- c(0)
for (i in 1:F){
  Y1.test <- rbind(Y1.test, data1.new.test[[i]])
} 
Y1.test <- Y1.test[-1,]


Y2.test <- c(0)
 for (i in 1:F){
   Y2.test <- rbind(Y2.test, data2.new.test[[i]])
 } 
 Y2.test <- Y2.test[-1,]
 
 
 
 
 
 ###########################################################################################
 ## True Labels for the points
 ### Training Data
 c.true <- c(0)
 for ( i in 1:F){
   c.true <- rbind(as.matrix(c.true) , as.matrix(c(rep(i, as.integer(N.train * p.dist[i])))))  
 }
 c.true <- as.factor(c.true[-1,])
 ### Testing Data Set
 c.true.new <- c(0)
 for ( i in 1:F){
   c.true.new <- rbind(as.matrix(c.true.new) , as.matrix(c(rep(i, as.integer(N.test * p.dist[i])))))  
 }
 c.true.new <- as.factor(c.true.new[-1,])
 
 
## All features PCA plots
  pc <- prcomp(Y1)
  pc.pred <- predict(pc,newdata = Y1)
  plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c.true, main = 'Main Training Data 1')
  
  pc <- prcomp(Y1.test)
  pc.pred <- predict(pc,newdata = Y1.test)
  plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c.true.new, main = 'Main Testing Data 1')
 
 ## All features PCA plots
 pc <- prcomp(Y2)
 pc.pred <- predict(pc,newdata = Y2)
 plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c.true, main = 'Main Training Data 2')
 
 pc <- prcomp(Y2.test)
 pc.pred <- predict(pc,newdata = Y2.test)
 plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c.true.new, main = 'Main Testing Data 2')
  
#####################################################################################################
#####################################################################################################  
#####################################################################################################
#####################################################################################################
######## Pure Survival Times ########################3
####Now WE DEAL WITH CLUSTERED DATA AND GENERATE NON RELEVANT FEATURES INDEPENENTLY ###
## Selcting the beta co-efficients
## The Co-efficients have to be obtained from uniform distribution between [-3,3]
  beta1.list <- list(0)
  for ( i in 1:F){
    beta1.list[[i]] <- runif(rel.D1,-3,3)
  }
  
  beta2.list <- list(0)
  for ( i in 1:F){
  beta2.list[[i]] <- runif(rel.D2,-3,3)
  }

  
  ## The pure time is generated for TRAINING DATA
  time1.pur.list <- list(0)
  for ( i in 1:F){
    time1.pur.list[[i]] <- t(beta1.list[[i]]) %*% t(Y1.rel.sc.list[[i]])
  }
  
  ## The pure time is generated for TESTING DATA
  time1.pur.list.test <- list(0)
  for ( i in 1:F){
    time1.pur.list.test[[i]] <- t(beta1.list[[i]]) %*% t(Y1.rel.sc.list.test[[i]])
  }
  
## The pure time is generated for TRAINING DATA
time2.pur.list <- list(0)
for ( i in 1:F){
  time2.pur.list[[i]] <- t(beta2.list[[i]]) %*% t(Y2.rel.sc.list[[i]])
}

## The pure time is generated for TESTING DATA
time2.pur.list.test <- list(0)
for ( i in 1:F){
  time2.pur.list.test[[i]] <- t(beta2.list[[i]]) %*% t(Y2.rel.sc.list.test[[i]])
}

  
  
  
  
 
  
assign("rel.D1", rel.D1, envir = .GlobalEnv)
assign("irrel.D1", irrel.D1, envir = .GlobalEnv)
assign("Y1", Y1, envir = .GlobalEnv)
assign("Y1.test", Y1.test, envir = .GlobalEnv)
assign("c.true", c.true, envir = .GlobalEnv)
assign("c.true.new", c.true.new, envir = .GlobalEnv)
assign("beta1.true", beta1.list, envir = .GlobalEnv)
assign("time1.pur.list", time1.pur.list, envir = .GlobalEnv)
assign("time1.pur.list.test", time1.pur.list.test, envir = .GlobalEnv)
assign("rel.D2", rel.D2, envir = .GlobalEnv)
assign("irrel.D2", irrel.D2, envir = .GlobalEnv)
assign("Y2", Y2, envir = .GlobalEnv)
assign("Y2.test", Y2.test, envir = .GlobalEnv)
assign("beta2.true", beta2.list, envir = .GlobalEnv)
assign("time2.pur.list", time2.pur.list, envir = .GlobalEnv)
assign("time2.pur.list.test", time2.pur.list.test, envir = .GlobalEnv)
  
  
 #########################################################################################
 ##########################################################################################
 ##### SURVIVAL TIMES  GENERATION #####################################################################
 ##########################################################################################
 ###########################################################################################
 
##### Noise Terms

### Training Data
time.noise.list1 <- list(0)
for ( i in 1:F){
  time.noise.list1[[i]] <- rnorm(length(which(c.true == i)), mean = 20*i, sd = i )
}

time.noise.list2 <- list(0)
for ( i in 1:F){
  time.noise.list2[[i]] <- rnorm(length(which(c.true == i)), mean = 20*i, sd = i )
}


###### Test Data
time.noise.list.test1 <- list(0)
for ( i in 1:F){
  time.noise.list.test1[[i]] <- rnorm(length(which(c.true.new == i)), mean = 20*i, sd = i )
}

time.noise.list.test2 <- list(0)
for ( i in 1:F){
  time.noise.list.test2[[i]] <- rnorm(length(which(c.true.new == i)), mean = 20*i, sd = i )
}

###### Calcuate ############## 
###### Times Per Data source
time.list1 <- list(0)
for ( i in 1:F){
  time.list1[[i]] <- time1.pur.list[[i]] + time.noise.list1[[i]]
}

time.list2 <- list(0)
for ( i in 1:F){
  time.list2[[i]] <- time2.pur.list[[i]] + time.noise.list2[[i]]
}


time.list1.test <- list(0)
for ( i in 1:F){
  time.list1.test[[i]] <- time1.pur.list.test[[i]] + time.noise.list.test1[[i]]
}

time.list2.test <- list(0)
for ( i in 1:F){
  time.list2.test[[i]] <- time2.pur.list.test[[i]] + time.noise.list.test2[[i]]
}


sdlist1 <- list(0)
 for ( j in 1:F){
   sdlist1[[j]] <- sd(as.vector(time.list1[[j]]))
 }
 
sdlist2 <- list(0)
for ( j in 1:F){
  sdlist2[[j]] <- sd(as.vector(time.list2[[j]]))
}


 sdlist1.test <- list(0)
 for ( j in 1:F){
   sdlist1.test[[j]] <- sd(as.vector(time.list1.test[[j]]))
 }
 
sdlist2.test <- list(0)
for ( j in 1:F){
  sdlist2.test[[j]] <- sd(as.vector(time.list2.test[[j]]))
}



## The actual training survival times
 time.final <- list(0)
 for ( i in 1:F){
   time.final[[i]] <- ((sdlist1[[i]])^2 * as.vector(time.list1[[i]]) + (sdlist2[[i]])^2 * as.vector(time.list2[[i]])) / ((sdlist1[[i]])^2  + (sdlist2[[i]])^2)
 }

## The actual testing survival times
time.final.test <- list(0)
for ( i in 1:F){
  time.final.test[[i]] <- ((sdlist1.test[[i]])^2 * as.vector(time.list1.test[[i]]) + (sdlist2.test[[i]])^2 * as.vector(time.list2.test[[i]])) / ((sdlist1.test[[i]])^2  + (sdlist2.test[[i]])^2)
}


#######################################MAKING TIME from cluster data ########################################################
## Real time without Censoring
  
time.real <- c(0)
time.real <- as.vector(unlist(time.final))
  
time.real.test <- c(0)
time.real.test <- as.vector(unlist(time.final.test))
  
time <- c(0)
censoring <- c(rep(1,N.train))
time <- time.real
  
time.test <- c(0)
censoring.test <- c(rep(1,N.test))
time.test <- time.real.test
  
  
  
  
  ## Boxplots for Vizualization of the time Data without censoring
  boxplot(time.final, main = "Training Times")
  boxplot(time.final.test, main = "Testing Times")
  
  
  ### A Quick ManWhittney U / Kruskal test  test to check if the time's of the two cluster are significantly different
  ## wilcox.test(as.vector(time.list[[1]]), as.vector(time.list[[2]]), alternative = "two.sided")
  
 ### Visualization of the different survival curves
  surv.ob <- Surv(time,censoring)
  survfit <- survfit(surv.ob ~ c.true)
  plot(survfit, main = "Training")
  
  surv.ob <- Surv(time.test,censoring.test)
  survfit <- survfit(surv.ob ~ c.true.new)
  plot(survfit, main = "Testing")
  

 assign("time.real", time.real, envir = .GlobalEnv)
 assign("time.real.new", time.real.test, envir = .GlobalEnv)
 assign("time", time, envir = .GlobalEnv)
 assign("time.new", time.test, envir = .GlobalEnv)
 assign("censoring", censoring, envir = .GlobalEnv)
 assign("censoring.new", censoring.test, envir = .GlobalEnv)
}