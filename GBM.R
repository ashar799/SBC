###########################################################################################
### This applies TCGA data set for GBM on ourSBC model ###################################
rm(list = ls())

#### Load Data #########
load("/home/bit/ashar/ExpressionSets/TWO_VIEW/TCGA_GBM/FINAL/DataTCGA-GBM.RData")
Y1 <- relev$Y1
Y2 <- relev$Y2

Y1.test <- relev$Y1.test
Y2.test <- relev$Y2.test

pheno <- relev$pheno
pheno.test <- relev$pheno.test


## Number of points
N.train =  nrow(Y1)
N.test = nrow(Y1.test)

N <- N.train
## Number of Clusters
F = 2

######
D1 <- ncol(Y1)
D2 <- ncol(Y2)

####
time <- pheno$Survival
censoring <- pheno$Censoring

time.new <- pheno.test$Survival
censoring.new <- pheno.test$Censoring

###########
c.true <- pheno$Subtype
levels(c.true) <- c(1,2,3,4)

##########
c.true.new <- pheno.test$Subtype
levels(c.true.new) <- c(1,2,3,4)


############################# PARAMETERS for GIBB's SAMPLING ####
iter = 200
iter.burnin = 100
iter.thin  = 5
k = 4
K <-  as.integer(N)
Time <- cbind(time,censoring)


################# GroundTruth (by pasting togehter columns)
source('TRAINmultiComparison.R')
multigroundtruth()

########### Train the Model #########################################
begin.time <-  Sys.time()
source('burninmultiDPMM.R')
burninmultiDPMM()
end.time <- Sys.time()
########### Train the Model #########################################
source('gibbsmultiDPMM.R')
gibbsmultiDPMM()

##### Analyzing the Model #########################################
source('MCMCmultianalyze.R')
analyzemultiDPMM()


##### Predicting the Classes for new Data Points####################################
source('multipredictCLASS.R')
multipredictCLASS(Y1.test, Y2.test)
predicted.rindex <- adjustedRandIndex(c.true.new,apply(posteriorprob,1,which.max))


###### Predicting Survival Times Based on two Molecular Data Sources ####################################
source('multipredictTIME.R')
multipredictchineseAFTtime(Y1.test, Y2.test)
predicted.cindex <- survConcordance(Surv(exp(time.new),censoring.new) ~ exp(-post.time.avg))[1]

########### Check Prediction Ground Truth
source('TESTmultiComparison.R')
predictionGroundTruth()
