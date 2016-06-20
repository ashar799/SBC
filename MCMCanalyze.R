
############# TRAININIG DATA ###########################################
#######################################################################
########################################################################
#### This function calculates some important metrices for the TRAINING DATA Data
#### Adjusted RandIndex
#### C-Index
#### Brier Scores for the Survival Curves
#### Point Estimate of Clsuter Assignments
##### Brier Scores for feature Selection


MCMCanalyze = function(){
  
########## ANLAYSING THE MCMC samples AND CALCULATING METRICES #######################################################
source('linearprediction.R')


Nps = as.integer(iter/ iter.thin)
count <- Nps
final.rand <- c(0)

############ The Matrices that will store the results #################################################
final.rand <- c(0)
cindex.final <- c(0)
brier.final <- matrix(NA, nrow = Nps, ncol = F)

################ Begin Analysig the MCMC samples #######################################################

require(Hmisc) 

for (h in 1:Nps){
### Adjusted Rand Indices
# final.rand[h] <- adjustedRandIndex(c.list[[h]],as.factor(c.true))

### See C-Index (concordance index)
surv.aft <- Surv(exp(time),censoring)
library(Hmisc)
### Predict Time from the model
tem.tim <- as.vector(unlist(predicttime(c.list[[h]], Y, That.list[[h]], Time, beta0.list[[h]], betahat.list[[h]], sigma2.list[[h]])))
library(Hmisc)
cindex.final[h] <-  survConcordance(surv.aft ~ exp(-tem.tim))[[1]]

### Brier Scores
# for ( v in 1:F){
#   ind <- which((c.list[[h]] == v))
#   time.tmp <- time[ind]
#   censoring.tmp <- censoring[ind]
#   Y.tmp <- Y[ind,]
#   rownames(Y.tmp) <- as.character(c(1:nrow(Y.tmp)))
#   smod <-  Surv(exp(time.tmp), censoring.tmp)
#   L = length(ind)
#   linear.pred <- tem.tim[ind]
#   sigma.pred <- log(sqrt(sigma2.list[[h]][v]))
#   ### The survival function in AFT model
#   S1 <- function (times = NULL, lp = NULL, parms = sigma.pred) 
#   {
#     t.trans <- logb(times)
#     names(t.trans) <- format(times)
#     1 - pnorm((t.trans - lp)/exp(parms))
#   }
#   mat.tmp <- matrix(NA, nrow = L, ncol = L)
#   for (j in 1:L){
#     mat.tmp[,j] <- S1(exp(time.tmp[j]),lp =tem.tim[ind])
#   }
#   brier.final[h,v] <- sbrier(smod,mat.tmp, exp(time.tmp))[1] 
#   
# }

}
###############################################
###### Calculating POINT ESTIMATES ############
###############################################





##### Class Assignments ########################
c.matrix <- matrix(NA, nrow = N, ncol = count)
for ( i in 1:count){
  c.matrix[,i] <- c.list[[i]]
}


c.final <- c.matrix[,count/2]

active <- as.numeric(rownames(table(c.final)))

active <- active[(active %%1 ==0)]
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



final.rand <<- final.rand
cindex.final <<- cindex.final
brier.final <<- brier.final
c.final <<- c.final
final.betahat <<- final.betahat 

}


