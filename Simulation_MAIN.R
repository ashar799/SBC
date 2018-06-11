### This is the Main Function and contains a simulation case
### Also CHECK THE TIME REQUIRED FOR THE MODEL
setwd("~/Dropbox/Code/DPmixturemodel/SBC")

rm(list = ls())
#################################### SIMULATED DATA PROPERTIES ####################################################
## Number of points
N.test = 500
N.train = 500

## Number of Clusters
F = 2
k =F
N = N.train

## Distribution of the points within three clusters

p.dist = c(0.5,0.5)

## Total Number of features D

D = 20

## Total Percentage of irrelevant feature
prob.noise.feature = 0.50


## Overlap between Cluster of molecular Data of the relevant features
prob.overlap = 0.05

###### Get the Data #####################################

## Initialize the Training Data
source('simulate.R')
simulate()

####### Assign training and testing data ###############
Y <- Y.dat
Y.new <- Y.new.dat

smod <-  Surv(exp(time), censoring)
smod.new <- Surv(exp(time.new), censoring.new)


##################### STATE OF THE ART TECHNIQUES #################################
##################### BASIC METHODS + SOME ADVANCED METHODS ########################
source('Comparisonx.R')
Comparisonx()

source('ComparisionFLX.R')
ComparisionFLX()

source('ComparisionPReMiuM.R')
ComparisionPReMiuM()
setwd("~/Dropbox/Code/DPmixturemodel/SBC")




######################### Initialize the Parameters ##############################
source('initialize.R')
initialize()

###################### Start with a good configuration ###########################
source('startSBC.R')
startSBC()


############################# PARAMETERS for GIBB's SAMPLING ####
iter = 20
iter.burnin = 50
iter.thin  = 2
Nps = as.integer(iter/iter.thin)


########### Train the Model #########################################
source('burninDPMM.R')
burninDPMM()

source('gibbsDPMM.R')
gibbsDPMM()
  
########## Analyze the fit ##########################################
### Good feature selection from heatmap plus cindex plus randindex
source('MCMCanalyze.R')
MCMCanalyze()
recovRandIndex.sbc <<-  as.numeric(adjustedRandIndex(c.true, c.sbc))  

pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.sbc))) + ggtitle(" SBC Clustering \n Test Set") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 





######## Predict on New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
source('predictCLASS.R')
predictCLASS(Y.new)
## Check the predicted Rand Index 

predRandIndex.sbc <<-  as.numeric(adjustedRandIndex(c.true.new, c.sbc.new))  
c.sbc.new.knn <<- knn(train = Y, test = Y.new, cl = c.sbc, k = F)


pc <- prcomp(Y.new)
pc.pred <- predict(pc,newdata = Y.new)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.sbc.new))) + ggtitle(" SBC Clustering \n Test Set") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 

surv.fit <- survfit(smod.new ~ c.sbc.new)
logrank <- survdiff(smod.new ~ c.sbc.new)
p5 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2),labels = c('1', '2'))

#### Choose that configuration which has the highest difference in survival curves
predRandIndex.sbc <- c(0)
for (j in 1:Nps){
  predRandIndex.sbc[j] <-  adjustedRandIndex(c.true.new,c.matrix.new[,j])
}

### Use SBC + knn to make predictions on future probabilities ####
cl.old <-  c.sbc
cl.new <- c.sbc.new.knn

### Use PAFT now to build the corresponding cluster specific PAFT models
pre.sbc <- c(0)

for ( q in 1:F){
  ind <- which(cl.old == q)
  ind.new <- which(cl.new == q)

  time.tmp <- time[ind]

  Y.tmp <- Y[ind,]
  Y.tmp.new <- Y.new[ind.new,]

  reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
  pre.sbc[ind.new] <- predict(object = reg, newx = Y.tmp.new, s = "lambda.min")
}
predCIndex.sbc.knn.aft  <<- as.numeric(survConcordance(smod.new ~ exp(-pre.sbc))[1])







source('predictTIME.R')
predictchineseAFTtime(Y.new)




