######## This function generates simulated data for the two data source case
generatedata = function(N,D,F,p.dist) {


  N = N
  
  ## Number of Clusters
  F = F
  
  ## Distribution of the points within three clusters
  
  p.dist = p.dist
  
  ## Total Number of features D
  
  D = D
  
  ## Overlap between Cluster of molecular Data of the relevant features
  
  prob.overlap = 0.05
  
  ## Percentage of Noise/Overlap in Time Data
  
  prob.noise.feature = 0.90
  
  ## Actual Number of Components and dimension  which are relevant
  rel.D = as.integer((D* (1-prob.noise.feature)))
  ## Actual Number of Irrelevant Componenets
  irrel.D = D - rel.D
  
  #######Generating the first data set
  
  
  ## Generating two Data with overlap only with the relevant features
  A <- MixSim(MaxOmega = prob.overlap ,K = F, p = rel.D, int =c(0,1),lim = 1e08)
  ## Generating 1st data set ##########################################################3

  data.mu = array(data = NA, dim =c(F,D))
  data.S = array(data = NA, dim =c(F,D,D))

for( i in 1:F){
  data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
  data.S[i,1:rel.D,1:rel.D] <- A$S[1:rel.D,1:rel.D,i]
}

## The relevant data is genereated first
Y.rel.list <- list(0)
for ( i in 1:F){
  Y.rel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = data.mu[i,1:rel.D], Sigma = data.S[i,1:rel.D,1:rel.D])
}

## Scaling the Data as ONLY the scaled data will be used for generating the times
Y.rel.sc.list <- list(0)
for ( i in 1:F){
  Y.rel.sc.list[[i]] <- scale(Y.rel.list[[i]], center = TRUE, scale = TRUE)
}

## Irrelevant features
Y.irrel.list <- list(0)
for ( i in 1:F){
  mean <- runif(irrel.D,-1.5,1.5)
  Y.irrel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = mean, Sigma = diag(x =1, nrow = irrel.D, ncol = irrel.D))
}

### Combining the data with relevant and irrelevant columns
data.old <- list(0) 
for (i in 1:F){
  data.old[[i]] <-  cbind(Y.rel.list[[i]], Y.irrel.list[[i]]) 
}

############################################### MAKING Y from the clusters data #####################3
Y.old <- c(0)
for (i in 1:F){
  Y.old <- rbind(Y.old, data.old[[i]])
} 
Y.old <- Y.old[-1,]

#########################################################################################
X <- Y.old

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

irrel.X <- matrix(NA, nrow = N, ncol = irrel.D)

for ( i in 1: irrel.D){
  
  matTemp <- matrix(runif(n = (rk -rel.D), min = -1.5, max= 1.5),  nrow = (rk-rel.D), ncol =1)
  irrel.X[,i] <- as.vector(linear.space %*% matTemp)
  
}

## Checking if the covariance is indeed small

cov.mat <- cov(rel.X,irrel.X)

boxplot(cov.mat)

## Building the full data matrix

X.full <- cbind(rel.X, irrel.X)


levelplot(cov(X.full[,1:D]))

Y1 <- X.full
## Selcting the beta co-efficients for the first data

## The Co-efficients have to be obtained from uniform distribution between [-3,3]
beta.list <- list(0)
half <- rel.D/2
ohalf <- rel.D - half
for ( i in 1:F){
  beta.list[[i]] <- as.vector(rbind(runif(half, min = -3, max = -0.1), runif(ohalf, min = 0.1, max = 3)))
}

## The pure time is generated
time.pur.list <- list(0)
for ( i in 1:F){
  time.pur.list[[i]] <- t(beta.list[[i]]) %*% t(Y.rel.sc.list[[i]]) 
}

beta.irrel.list <- list(0)
for ( i in 1:F){
  beta.irrel.list[[i]] <- rep(0,irrel.D)
}
                        
beta <- list(0)
for ( i in 1:F){
  beta[[i]] <- as.vector(c(beta.list[[i]], beta.irrel.list[[i]]))
}



list('Y' = Y.old, 'beta' = beta, 'timepur'= time.pur.list)
}
