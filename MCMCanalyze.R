############# TRAININIG DATA ###########################################
#######################################################################
########################################################################
#### This function calculates some important metrices for the TRAINING DATA Data
#### Adjusted RandIndex
#### C-Index
#### Point Estimate of Clsuter Assignments
##### Brier Scores for feature Selection


MCMCanalyze = function(){
  
########## ANLAYSING THE MCMC samples AND CALCULATING METRICES #######################################################
source('linearprediction.R')
count <- Nps
final.rand <- c(0)

############ The Matrices that will store the results #################################################
##### Class Assignments ########################
c.matrix <- matrix(NA, nrow = N, ncol = count)
recovCIndex.sbc <- c(0)
recovCIndex.sbc.paft <- c(0)
################ Begin Analysig the MCMC samples #######################################################

require(Hmisc) 

for (h in 1:Nps){
### See C-Index (concordance index) returned by the model
surv.aft <- Surv(exp(time),censoring)
### Predict Time from the model
tem.tim <- predicttime(c.list[[h]], Y, That.list[[h]], Time, beta0.list[[h]], betahat.list[[h]], sigma2.list[[h]])$predicttime
library(Hmisc)
recovCIndex.sbc[h] <-  survConcordance(surv.aft ~ exp(-tem.tim))[[1]]

## try to fit a Cluster-specific AFT model on top of it ##

###### penAFT ###################################################################
######## Penalized AFT with k-means clustering ######################################################
sbc.aft <- c(0)
for ( q in 1:F){
  ind <- which((c.list[[h]]) == q)
  L= length(ind)
  
  time.tmp <- time[ind]
  censoring.tmp <- censoring[ind]
  Y.tmp <- Y[ind,]
  
  reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
  coeff.pred <- coef(object =reg, newx = Y.tmp, s= "lambda.min")
  sbc.aft[ind] <- predict(object = reg, newx = Y.tmp, s = "lambda.min") 
}
recovCIndex.sbc.paft[[h]] <- as.numeric(survConcordance(smod ~ exp(-sbc.aft))[1])

c.matrix[,h] <- c.list[[h]]

}




###############################################
###### Calculating POINT ESTIMATES ############
###############################################
psm <- comp.psm(t(c.matrix))
mpear <- maxpear(psm)

### If we build a cluster specific sbc approach
c.final <<- mpear$cl
c.sbc <<- mpear$cl

active <- as.numeric(rownames(table(c.final)))


############ Time Covariate Slopes FOR Relevant Clusters ############
list.betahat <- list(0)

count = Nps

for ( i in 1:count){
  list.betahat[[i]] <- (betahat.list[[i]][active,] != 0) +0
}

Q <- length(active)
matrix.betahat <- array(data = NA, dim =c(Q,count,D))

for ( z in 1:Q){
  for ( x  in 1:count){
    matrix.betahat[z,x,] <- list.betahat[[x]][z,]
    }
}

final.betahat <- apply(matrix.betahat,c(1,3),mean)
### Probability of betahat of genes FOR ONE SIMULATION
##colnames(final.betahat) =  c(rep("relevant",rel.D),rep("irrelevant",irrel.D))
heatmapdata <- as.data.frame(final.betahat)
heatmap.2(t(as.matrix(heatmapdata)),dendrogram="none", margins=c(6,10), main = "Posterior prob. \n for Selection \n ", cexCol = 0.85, cexRow = 0.7, Rowv = FALSE)

recovCIndex.sbc <<- recovCIndex.sbc
recovCIndex.sbc.paft <<- recovCIndex.sbc.paft
c.final <<- c.final
final.betahat <<- final.betahat 

}


