### This is the Main Function and contains a simulation case
### It also checks the feature Importance


rm(list = ls())
#################################### SIMULATED DATA PROPERTIES ####################################################
## Number of points
N.test = 100
N.train = 100

## Number of Clusters
F = 2

## Distribution of the points within three clusters

p.dist = c(0.5,0.5)

## Total Number of features D

D = 30

## Total Percentage of irrelevant feature
prob.noise.feature = 0.2


## Overlap between Cluster of molecular Data of the relevant features
prob.overlap = 0.01




###### Get the Data #####################################

## Initialize the Training Data
source('simulate.R')
simulate()

### For Getting feature importance
response <- c(rep(1,rel.D),rep(0,irrel.D))


####### Assign training and testing data ###############
Y <- Y.dat


###############  Feature Importance with Penalized FlexMix #######################################################
################################################################################
data <- data.frame(time, Y)
## The cross validation folds for choosing lambda
fo <- sample(rep(seq(10), length = nrow(data)))
gr.km <- kmeans(Y, F, nstart =10)
gr.flx <- flexmix(time ~ ., data = data, k = F, cluster = gr.km$cluster, model = FLXMRglmnet(foldid = fo, adaptive= FALSE), control = list(iter.max = 500))
flx.features <- abs((unlist(parameters(gr.flx, component = 1, model = 1))[D+2])^-1 * unlist(parameters(gr.flx, component = 1, model = 1))[1:D] + (unlist(parameters(gr.flx, component = 2))[D+2])^-1 * unlist(parameters(gr.flx, component = 2))[1:D])
auc.flx <- as.numeric(auc(roc(response,flx.features)))




# ########### Feature Importance from sparse K-means #########################
# km.perm <- KMeansSparseCluster.permute(x = Y, K= F ,wbounds=c(1.5,2:6),nperms=20)
# km.out <- KMeansSparseCluster(x = Y, K=F,wbounds=km.perm$bestw)
# auc.SKM <-  as.numeric(auc(roc(response, km.out[[1]]$ws)))
# 
# 
# ########### Feature Importance from sparse hierarchical clustering #########################
# perm.out <- HierarchicalSparseCluster.permute(x = Y, wbounds=c(1.5,2:6), nperms=20)
# # Perform sparse hierarchical clustering
# sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists, wbound=perm.out$bestw, method="complete")
# auc.SHC <-  as.numeric(auc(roc(response, as.vector(sparsehc$ws))))








############################# PARAMETERS for GIBB's SAMPLING ####
iter = 50
iter.burnin = 50
iter.thin  = 2
k = F
Nps = iter.burnin/iter.thin

######################### Initialize the Parameters ##############################
source('initialize.R')
initialize()


##################### OPTIONAL COMPARISON WITH KNOWN METHODS ######################
######### BASIC METHODS + SOME ADVANCED METHODS ############################################
#source('TRAINComparisonx.R')
#TRAINComparison()


########### Train the Model #########################################
source('burninDPMM.R')
burninDPMM()

source('gibbsDPMM.R')
gibbsDPMM()

########## Analyze the fit ##########################################
### Good feature selection from heatmap plus cindex plus randindex
source('MCMCanalyze.R')
MCMCanalyze()



######### Feature Importance ###########################################
### SBC Model  ###############################
#### ############
Nps <- Nps -1
sum.diff <- matrix(0, nrow = D, ncol =Nps)
for ( i in 1:Nps){
  sum.diff[1:D,i] <- abs((mu.list[[i]][1,] - mu.list[[i]][2,])/diag(W.list[[i]]))  
}


auc.sbc <- c(0)
for( i in 1:Nps){
  auc.sbc[i]  <- as.numeric(auc(roc(response,sum.diff[,i])))
}

