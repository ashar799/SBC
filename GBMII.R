##############################################################################################
### This applies TCGA data set for GBM on our SBC model ######################################
### To avoid possible changes with the two data source case use CCA processed feature sets ###
rm(list = ls())

#### Load Data #########
load("/home/bit/ashar/ExpressionSets/TWO_VIEW/TCGA_GBM/CCA/DataTCGA-GBM-CCA.RData")

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
F = 4


Nps = as.integer(iter/ iter.thin)

######
D1 <- ncol(Y1)
D2 <- ncol(Y2)

####
time <- pheno$Survival
censoring <- pheno$Censoring
smod <- Surv(exp(time),censoring)



time.new <- pheno.test$Survival
censoring.new <- pheno.test$Censoring
smod.new <- Surv(exp(time.new), censoring.new)



###########
c.true <- pheno$Subtype
levels(c.true) <- c(1,2,3,4)

##########
c.true.new <- pheno.test$Subtype
levels(c.true.new) <- c(1,2,3,4)


############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 200
iter.thin  = 2
k = 4
K <-  as.integer(N)
Time <- cbind(time,censoring)

source('multiComparisonx.R')
multiComparisonx()

######################### Initialize the Parameters ################
source('multiinitialize.R')
multiinitialize()

########### Train the Model #########################################
source('burninmultiDPMM.R')
burninmultiDPMM()

########### Train the Model #########################################
source('gibbsmultiDPMM.R')
gibbsmultiDPMM()

##### Analyzing the Model #########################################
source('MCMCmultianalyze.R')
analyzemultiDPMM()


###############################
### Some plots and analysis ####
#################################
#### Generating some plots and results ON Training data ###
logrank <- survdiff(smod ~ c.final)
1 - pchisq(unlist(logrank)$chisq,df =3)
surv.fit <- survfit(smod ~ c.final)
p5 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators \n GBM II Cancer Data Set") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2,3,4),labels = c('1', '2','3','4'))

############ Generating some Plots ##########################
c.final <- c.kmeans

pc1 <- prcomp(Y1)
pc.pred1 <- predict(pc1,newdata = Y1)
p1 <- ggplot(as.data.frame(pc.pred1), aes(x=pc.pred1[,1], y= pc.pred1[,2], colour= as.factor(c.final))) + ggtitle(" SBC Clustering \n GBM II Cancer Data Set \n mRNA") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 

pc2 <- prcomp(Y2)
pc.pred2 <- predict(pc2,newdata = Y2)
p2 <- ggplot(as.data.frame(pc.pred2), aes(x=pc.pred2[,1], y= pc.pred2[,2], colour= as.factor(c.final))) + ggtitle(" SBC Clustering \n GBM II Cancer Data Set \n miRNA ") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 


source('multiplot.R')
multiplot(p1,p2,cols=2)

#### Use the adhoc-prediction algorithm ####################
source('predictADHOCGBMII.R')


logrank <- survdiff(smod.new ~ c.sbc.new)
surv.fit <- survfit(smod.new ~ c.sbc.new)
p2 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators \n GBM II Cancer Test Data Set") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2,3,4),labels = c('1', '2','3','4'))
############ Generating some Plots ##########################



###### Predicting Survival Times ####################################
source('multipredictTIME.R')
multipredictchineseAFTtime(Y1.test, Y2.test)

#### Choose that configuration which has the highest difference in survival curves ####################
lr <- c(0)
for (j in 1:Nps){
  lr[j] <-  1 - pchisq(unlist(survdiff(smod.new ~ c.matrix.new[,j]))$chisq,df =1)
}
c.sbc.new <- c.matrix.new[,which.min(lr)]





##### Predicting class labels #####################
source('multipredictCLASS.R')
multipredictCLASS(Y1.test,Y2.test)



