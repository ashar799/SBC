### This is the Main Function and contains a simulation case
### Also CHECK THE TIME REQUIRED FOR THE MODEL


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

D = 20

## Total Percentage of irrelevant feature
prob.noise.feature = 0.2


## Overlap between Cluster of molecular Data of the relevant features
prob.overlap = 0.01

###### Get the Data #####################################

## Initialize the Training Data
source('simulate.R')
simulate()

####### Assign training and testing data ###############
Y <- Y.dat
Y.new <- Y.new.dat

############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 100
iter.thin  = 5
k = F
Nps = iter.burnin/iter.thin

######################### Initialize the Parameters ##############################
source('initialize.R')
initialize()


##################### OPTIONAL COMPARISON WITH KNOWN METHODS ######################
######### BASIC METHODS + SOME ADVANCED METHODS ############################################
source('TRAINComparisonx.R')
TRAINComparison()


########### Train the Model #########################################
source('burninDPMM.R')
burninDPMM()

source('gibbsDPMM.R')
gibbsDPMM()

########## Analyze the fit ##########################################
### Good feature selection from heatmap plus cindex plus randindex
source('MCMCanalyze.R')
MCMCanalyze()


######## Predict on New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
source('predictCLASS.R')
predictCLASS(Y.new)
## Check the predicted Rand Index 


source('predictTIME.R')
predictchineseAFTtime(Y.new)
### Check of the Predicted C-index 
predicted.cindex <- survConcordance(Surv(exp(time.new),censoring.new) ~ exp(-post.time.avg))[1]

##################### OPTIONAL COMPARISON WITH KNOWN METHODS ######################
######### BASIC METHODS + SOME ADVANCED METHODS ############################################
source('TESTCOMPARISON.R')
TESTCOMPARISON()



