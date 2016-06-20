### This simulation Checks if the model works in 2 View Case
#### Copies from Simulation_Main from the one View Case

rm(list = ls())
#################################### SIMULATED DATA PROPERTIES ####################################################
## Number of points
N.test =  100
N.train = 100

N <- N.test
## Number of Clusters
F = 2

## Distribution of the points within three clusters

p.dist = c(0.5,0.5)

## Total Number of features D

D1 = 20
D2 = 20

## Total Percentage of irrelevant feature
prob.noise.feature = 0.20

## Overlap between Cluster of molecular Data of the relevant features
prob.overlap = 0.01

###### Get the Data #####################################

## Initialize the Training Data
source('multi_simulate.R')
multi_simulate()


############################# PARAMETERS for GIBB's SAMPLING ####
iter = 50
iter.burnin = 50
iter.thin  = 5
k = 2
K <-  as.integer(N)
Time <- cbind(time,censoring)


######################### Initialize the Parameters ################
source('multiinitialize.R')
multiinitialize()


################# GroundTruth (by pasting togehter columns)
source('TRAINmultiComparison.R')
multigroundtruth()

########### Train the Model #########################################
source('burninmultiDPMM.R')
burninmultiDPMM()

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


