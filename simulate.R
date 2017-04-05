simulate = function(){
  

  ## Actual Number of Components and dimension  which are relevant
  rel.D  <- as.integer(D* (1-prob.noise.feature))
  ## Actual Number of Irrelevant Componenets
  irrel.D <-  D - rel.D
  
#################################################################################################
###############################################################################################
########## Molecular Data Generatioin ###########################################################
##################################################################################################
#################################################################################################
  
  
  ## The data
  A <- MixSim(BarOmega = prob.overlap ,K = F, p = rel.D, int =c(-1.0,1.0), lim = 1e09)
  data.mu = array(data = NA, dim =c(F,rel.D))
  data.S = array(data = NA, dim =c(F,rel.D,rel.D))
 
# ## Data coming from a hypothetical population diagonal  precision matrix
#  for( i in 1:F){
#   data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
#   data.S[i,1:rel.D,1:rel.D] <- diag(diag(A$S[1:rel.D,1:rel.D,i]))
# }
# Data coming from a hypothetical population NON diagoanl matrix
for( i in 1:F){
  data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
  data.S[i,1:rel.D,1:rel.D] <- A$S[1:rel.D,1:rel.D,i]
}



Y.rel.list <- list(0)
for ( i in 1:F){
  Y.rel.list[[i]] <- mvrnorm(n = as.integer(N.train * p.dist[i]), mu = data.mu[i,1:rel.D], Sigma = data.S[i,1:rel.D,1:rel.D])
}

Y.rel.list.test <- list(0)
for ( i in 1:F){
  Y.rel.list.test[[i]] <- mvrnorm(n = as.integer(N.test * p.dist[i]), mu = data.mu[i,1:rel.D], Sigma = data.S[i,1:rel.D,1:rel.D])
}




## Scaling the Data as ONLY the scaled data will be used for generating the times
Y.rel.sc.list <- list(0)
Y.rel.sc.list.test <- list(0)
for ( i in 1:F){
  obj <- scale(Y.rel.list[[i]], center = TRUE, scale = TRUE)
  Y.rel.sc.list[[i]] <- scale(Y.rel.list[[i]], center = attr(obj,"scaled:center"), scale = attr(obj,"scaled:scale"))
  Y.rel.sc.list.test[[i]] <- scale(Y.rel.list.test[[i]], center = attr(obj,"scaled:center"), scale = attr(obj,"scaled:scale"))
}




## Irrelevant features
Y.irrel.list <- list(0)
for ( i in 1:F){
  mean <- runif(irrel.D,-1.5,1.5)
  Y.irrel.list[[i]] <- mvrnorm(n = as.integer(N.train * p.dist[i]), mu = mean, Sigma = diag(x =1, nrow = irrel.D, ncol = irrel.D))
}

Y.irrel.list.test <- list(0)
for ( i in 1:F){
  mean <- runif(irrel.D,-1.5,1.5)
  Y.irrel.list.test[[i]] <- mvrnorm(n = as.integer(N.test * p.dist[i]), mu = mean, Sigma = diag(x =1, nrow = irrel.D, ncol = irrel.D))
}


### Combining the data with relevant and irrelevant columns
data.old <- list(0) 
for (i in 1:F){
  data.old[[i]] <-  cbind(Y.rel.list[[i]], Y.irrel.list[[i]]) 
}

data.old.test <- list(0) 
for (i in 1:F){
  data.old.test[[i]] <-  cbind(Y.rel.list.test[[i]], Y.irrel.list.test[[i]]) 
}

data.new <- list(0)
data.new.test <- list(0)


######################################################################################### Making the irrelevant features independent from the dependent features #############
############ Training Data #####################################################
##################################################################################
for (f in 1:F){

X <- data.old[[f]]

rel.X <- as.matrix(X[,1:rel.D])

obj.qr <- qr(X)

rk <- obj.qr$rank

alpha <- qr.Q(obj.qr)[,1:rel.D]

gamma <- qr.Q(obj.qr)[,(1+rel.D):rk]

matT <- matrix(runif(n = rel.D*(rk -rel.D), min = -0.0005, max= 0.0005), nrow = rel.D, ncol = (rk -rel.D))

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

data.new[[f]] <- X.full

}



################################################################################################
############################### TESTING DATA ####################################################
for (f in 1:F){
  
X <- data.old.test[[f]]


rel.X <- as.matrix(X[,1:rel.D])

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

#boxplot(cov.mat)

## Building the full data matrix

X.full <- cbind(rel.X, irrel.X)


levelplot(cov(X.full))

Y.test <- X.full

data.new.test[[f]] <- X.full
}


############################################### MAKING Y from the clusters data #####################
Y <- c(0)
for (i in 1:F){
  Y<- rbind(Y, data.new[[i]])
} 
Y<- Y[-1,]

Y <- as.matrix(Y)
############################################### MAKING Y from the clusters data #####################
Y.test <- c(0)
for (i in 1:F){
  Y.test <- rbind(Y.test, data.new.test[[i]])
} 
Y.test <- Y.test[-1,]
Y.test <- as.matrix(Y.test)


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
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c.true, main = 'Main Training Data')

pc <- prcomp(Y.test)
pc.pred <- predict(pc,newdata = Y.test)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c.true.new, main = 'Main Testing Data')
###################################################################################################


#########################################################################################
##########################################################################################
##### SURVIVAL TIMES #####################################################################
##########################################################################################
###########################################################################################


####Now WE DEAL WITH CLUSTERED DATA AND GENERATE NON RELEVANT FEATURES INDEPENENTLY ###


## Selcting the beta co-efficients

## The Co-efficients have to be obtained from uniform distribution between [-3,3]
beta.list <- list(0)
for ( i in 1:F){
  beta.list[[i]] <- runif(rel.D,-3,3)
}


## The pure time is generated for TRAINING DATA
time.pur.list <- list(0)
for ( i in 1:F){
  time.pur.list[[i]] <- t(beta.list[[i]]) %*% t(Y.rel.sc.list[[i]])
}

## The pure time is generated for TESTING DATA
time.pur.list.test <- list(0)
for ( i in 1:F){
  time.pur.list.test[[i]] <- t(beta.list[[i]]) %*% t(Y.rel.sc.list.test[[i]])
}



time.noise.list <- list(0)
for ( i in 1:F){
  time.noise.list[[i]] <- rnorm(length(which(c.true == i)), mean = 20*i, sd = 0.1 )
}
time.noise.list.test <- list(0)
for ( i in 1:F){
  time.noise.list.test[[i]] <- rnorm(length(which(c.true.new == i)), mean = 20*i, sd = 0.1 )
}




time.list <- list(0)
for ( i in 1:F){
  time.list[[i]] <- time.pur.list[[i]] + time.noise.list[[i]]
}

time.list.test <- list(0)
for ( i in 1:F){
  time.list.test[[i]] <- time.pur.list.test[[i]] + time.noise.list.test[[i]]
}

#######################################MAKING TIME from cluster data ########################################################
## Real time without Censoring
time.real <- c(0)
for (i in 1:F){
  time.real <- cbind(time.real, time.list[[i]])
} 
time.real <- time.real[,-1]
time.real <- as.vector(unlist(time.real))

time.real.test <- c(0)
for (i in 1:F){
  time.real.test <- cbind(time.real.test, time.list.test[[i]])
} 
time.real.test <- time.real.test[,-1]
time.real.test <- as.vector(unlist(time.real.test))




####################################### Adding CENSORING INFORMATION  ################################################################
## Adding the censoring information to the TIME

# censoring <- rbinom(n = NROW(Y), size =1, prob = 1- prob.censoring)
# right.censoring.time <- min(time.real)  
# 
# time <- time.real
# 
# index.time <- which(censoring==0)
# for ( q in 1:length(index.time)){
#   time[index.time[q]] <- right.censoring.time
#   
# }

time <- c(0)
censoring <- c(rep(1,N.train))
time <- time.real

time.test <- c(0)
censoring.test <- c(rep(1,N.test))
time.test <- time.real.test




## Boxplots for Vizualization of the time Data without censoring
boxplot(time.list, main = "Training Times")
boxplot(time.list.test, main = "Testing Times")


### A Quick ManWhittney U / Kruskal test  test to check if the time's of the two cluster are significantly different
## wilcox.test(as.vector(time.list[[1]]), as.vector(time.list[[2]]), alternative = "two.sided")

kruskal.test(time, c.true)
kruskal.test(time.test, c.true.new)


### Visualization of the different survival curves
surv.ob <- Surv(time,censoring)
survfit <- survfit(surv.ob ~ c.true)
plot(survfit, main = "Training")

surv.ob <- Surv(time.test,censoring.test)
survfit <- survfit(surv.ob ~ c.true.new)
plot(survfit, main = "Testing")

#### Check to see if the Actual Parameters are supplied if the c-index in 0.99
## Create the beta for the full data set
beta.list.updated <- list(0)
for ( i in 1:F ){
beta.list.updated[[i]] <- c(as.vector(beta.list[[i]]), rep(0,irrel.D))  
}

## Scaled Y matrix
Ytemp <-  matrix(NA, nrow = N.train, ncol = D)
for ( i in 1:F ){
clu <- which(c.true== i)
Ytemp[clu,1:D] <- scale(Y[clu,1:D], center = TRUE, scale = TRUE)
}

time.avg <- c(0)
for ( i in 1:F ){
  clu <- which(c.true== i)
  time.avg[i] <- mean(time[clu])
}
## List of predictors
predictor <- c(0)
for ( i in 1:N.train){
predictor[i] <- time.avg[c.true[i]] + as.numeric(beta.list.updated[[c.true[i]]] %*% Ytemp[i,1:D])
}
surv.aft <- Surv(exp(time),censoring)
cindex.true <-  survConcordance(surv.aft ~ exp(-predictor))[[1]]




assign("N",N.train, envir = .GlobalEnv)
assign("rel.D", rel.D, envir = .GlobalEnv)
assign("irrel.D", irrel.D, envir = .GlobalEnv)
assign("Y.dat", Y, envir = .GlobalEnv)
assign("Y.new.dat", Y.test, envir = .GlobalEnv)
assign("beta.list", beta.list, envir = .GlobalEnv)
assign("c.true", c.true, envir = .GlobalEnv)
assign("c.true.new", c.true.new, envir = .GlobalEnv)
assign("time.real", time.real, envir = .GlobalEnv)
assign("time.real.new", time.real.test, envir = .GlobalEnv)
assign("time", time, envir = .GlobalEnv)
assign("time.new", time.test, envir = .GlobalEnv)
assign("censoring", censoring, envir = .GlobalEnv)
assign("censoring.new", censoring.test, envir = .GlobalEnv)
assign("cindex.true", cindex.true, envir = .GlobalEnv)




}

