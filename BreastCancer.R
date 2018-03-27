############ This file takes in training and testing data For the Breast Cancer #################
#################################################################################################
rm(list =ls())
################################################
######## Load the Data from the folder #########
load("/home/bit/ashar/ExpressionSets/ONE_VIEW/breastcancer/Final/DataVijverOWNsignatureBEST.RData")

Y <- relev$Y.train
time <- log(relev$time.train)
censoring <- relev$censoring.train
c.true <- relev$c.true


Y.new <- relev$Y.test
time.new <- log(relev$time.test)
censoring.new <- relev$censoring.test
c.true.new <- relev$c.true.new


####### Impute missing values
Y <- impute.knn(Y)$data
Y.new <- impute.knn(Y.new)$data

### Scaling the data could also be an option ###
# Y <- scale(Y, center = TRUE, scale = TRUE)
# obj <- scale(Y, center = TRUE, scale = TRUE)
# Y.new <- scale(Y.new, center = attr(obj,"scaled:center"), scale = (attr(obj,"scaled:scale")))


N <- nrow(Y)
D <- ncol(Y)
N.new <- nrow(Y.new)

smod <-  Surv(exp(time), censoring)
smod.new <- Surv(exp(time.new), censoring.new)

############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 100
iter.thin  = 5
k = 2
F =k
Nps = as.integer(iter/ iter.thin)



#### Checking distribution of training and Testing 
Y.tot <- rbind(Y,Y.new)
pc <- prcomp(Y.tot)
pc.pred <- predict(pc,newdata = Y.tot)
bt <- c(rep(0,N),rep(1,N.new))
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(bt))) + ggtitle(" Batch Effects \n Breast Cancer Data Set") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 
p1


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
#source('startSBC.R')
#startSBC()



########### Train the Model #########################################
source('burninDPMM.R')
burninDPMM()

source('gibbsDPMM.R')
gibbsDPMM()



########## Analyze the fit ##########################################
### Good feature selection from heatmap plus cindex plus randindex
source('MCMCanalyze.R')
MCMCanalyze()



###############################
### Some plots and analysis ####
#################################
#### Generating some plots and results ON Training data ###
logrank <- survdiff(smod ~ c.final)
1 - pchisq(unlist(logrank)$chisq,df =1)
surv.fit <- survfit(smod ~ c.final)
p5 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators \n Breast Cancer Data Set") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2),labels = c('1', '2'))
############ Generating some Plots ##########################
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.final))) + ggtitle(" SBC Clustering \n Breast Cancer Data Set") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 



######## Predict CLASS MEMBERSHIP on  New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
source('predictCLASS.R')
predictCLASS(Y.new)

### Use MPEAR APPROACH to get POINT ESTIMATES #############################################################################
### Build A consensus clustering based on the posterior matrix for both training and testing labels
psm2 <- comp.psm(t(c.matrix.new))
mpear2 <- maxpear(psm2)
adjustedRandIndex(c.true.new,mpear2$cl)
### Generally the MPEAR output needs post-processing
### If we build a cluster specific sbc approach
c.sbc.new <- mpear2$cl

#### Choose that configuration which has the highest difference in survival curves ####################
lr <- c(0)
for (j in 1:Nps){
  lr[j] <-  1 - pchisq(unlist(survdiff(smod.new ~ c.matrix.new[,j]))$chisq,df =1)
}
c.sbc.new <- c.matrix.new[,which.min(lr)]


logrank.new <- survdiff(smod.new ~ c.sbc.new)
surv.fit <- survfit(smod.new ~ c.sbc.new)
p2 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators \n Breast Cancer Test Data Set") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2),labels = c('1', '2'))
############ Generating some Plots ##########################
pc <- prcomp(Y.new)
pc.pred <- predict(pc,newdata = Y.new)
p3 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.sbc.new))) + ggtitle(" SBC Clustering \n Breast Cancer Test Data Set") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 


####### Use Ad-hoc method ##################
source('predictADHOCvijver.R')


######## Predict C_INDEX on  New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
source('predictTIME.R')
predictchineseAFTtime(Y.new)




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







