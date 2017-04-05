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


N <- nrow(Y)
D <- ncol(Y)
N.new <- nrow(Y.new)
k =2
F =k


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
iter = 200
iter.burnin = 100
iter.thin  = 5
Nps = 40



######################### Initialize the Parameters ##############################
source('initialize.R')
initialize()

###################### Start with a good configuration ###########################
source('startSBC.R')
startSBC()



########### Train the Model #########################################
source('burninDPMM.R')
burninDPMM()

source('gibbsDPMM.R')
gibbsDPMM()



########## Analyze the fit ##########################################
### Good feature selection from heatmap plus cindex plus randindex
source('MCMCanalyze.R')
MCMCanalyze()

surv.ob <- Surv(exp(time),censoring)
logrank <- survdiff(surv.ob ~ c.final)

######## Predict on New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
source('predictCLASS.R')
predictCLASS(Y.new)
## Check the predicted Rand Index 




source('predictTIME.R')
predictchineseAFTtime(Y.new)
### Check of the Predicted C-index 
predCIndex.sbc <- as.numeric(survConcordance(Surv(exp(time.new),censoring.new) ~ exp(-post.time.avg))[1])

