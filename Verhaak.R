############ This file takes in training and testing data for Glioblastoma #################
#################################################################################################
rm(list =ls())
################################################
######## Load the Data from the folder #########
load("/home/bit/ashar/ownCloud/DPMM_RESULTS/ONE_VIEW/Verhark/Final/DataVerhaak.RData")
Y.train.prelim <- relev$Y.train.prelim
Y.test.prelim  <- relev$Y.test.prelim
signature.vk <- relev$signature
signature.dpmm <- relev$signature.dpmm
pheno.train <- relev$pheno.train
pheno.test <- relev$pheno.test


Y <- t(Y.train.prelim[as.character(signature.dpmm),])
time <- pheno.train[,3]
censoring <- pheno.train[,2]
c.true <- pheno.train[,4]
levels(c.true) <- c(1,2,3,4)

Y.new <- t(Y.test.prelim[as.character(signature.dpmm),])
time.new <- pheno.test[,3]
censoring.new <- pheno.test[,2]
c.true.new <- pheno.test[,4]
levels(c.true.new) <- c(1,2,3,4)

Y <- scale(Y, center = TRUE, scale = TRUE)
obj <- scale(Y, center = TRUE, scale = TRUE)
Y.new <- scale(Y.new, center = attr(obj,"scaled:center"), scale = (attr(obj,"scaled:scale")))

N <- nrow(Y)
D <- ncol(Y)
N.new <- nrow(Y.new)

smod <-  Surv(exp(time), censoring)
smod.new <- Surv(exp(time.new), censoring.new)

k =4
F =k

Y.tot <- rbind(Y,Y.new)
pc <- prcomp(Y.tot)
pc.pred <- predict(pc,newdata = Y.tot)
bt <- c(rep(0,N),rep(1,N.new))
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(bt))) + ggtitle(" Batch Effects \n GBM I Cancer Data Set") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 
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


############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 100
iter.thin  = 5
Nps = as.integer(iter/ iter.thin)


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
1 - pchisq(unlist(logrank)$chisq,df =3)
surv.fit <- survfit(smod ~ c.final)
p5 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators \n Verhaak Cancer Data Set") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2),labels = c('1', '2'))
############ Generating some Plots ##########################
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.final))) + ggtitle(" SBC Clustering \n VerhaakCancer Data Set") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 


# ######## Predict on New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
# source('predictCLASS.R')
# predictCLASS(Y.new)
# ## Check the predicted Rand Index 
# 
# #### Choose that configuration which has the highest difference in survival curves ####################
# lr <- c(0)
# for (j in 1:Nps){
#   lr[j] <-  1 - pchisq(unlist(survdiff(smod.new ~ c.matrix.new[,j]))$chisq,df =2)
# }
# c.sbc.new <- c.matrix.new[,which.min(lr)]
# 
# ## If we fit cluster-specific models
# pre.sbc <- c(0)
# for ( q in 1:F){
#   ind <- which(c.sbc == q)
#   ind.new <- which(c.sbc.new == q)
#   
#   time.tmp <- time[ind]
#   
#   Y.tmp <- Y[ind,]
#   Y.tmp.new <- Y.new[ind.new,]
#   
#   reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
#   pre.sbc[ind.new] <- predict(object = reg, newx = Y.tmp.new, s = "lambda.min") 
# }
# predCIndex.sbc.aft  <<- as.numeric(survConcordance(smod.new ~ exp(-pre.sbc))[1])
# 

#### Use the adhoc-prediction algorithm ####################
source('predictADHOCverhaak.R')

surv.fit <- survfit(smod.new ~ c.sbc.new)
p2 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators \n Verhaak Cancer Test Data Set") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2,3,4),labels = c('1', '2','3','4'))
############ Generating some Plots ##########################
pc <- prcomp(Y.new)
pc.pred <- predict(pc,newdata = Y.new)
p3 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.sbc.new))) + ggtitle(" SBC Clustering \n Verhaak Cancer Test Data Set") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 


source('predictTIME.R')
predictchineseAFTtime(Y.new)

