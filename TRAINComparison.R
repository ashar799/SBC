############ Ground Truth on TRAINING DATA ###################################
#############################################################################
#############################################################################


### K-means + CoxPH
### K-means + AFT

### K-means + Penalized CoxPH
### K-means + Penalized AFT

### FlexMix +  CoxPH
### FlexMix +  AFT

### Mixture of Factor Analyzers
#### Sparse KMeans
#### Sparse Hierarchical clustering 


TRAINComparison = function(){
  
smod <-  Surv(exp(time), censoring)
  
# ######## Fitting a  Cox Proportional Hazards Model without Clustering
# fit.cox <- coxph(smod ~ Y[,1:D], data = as.data.frame(Y))
# cindex.cox <- survConcordance(smod ~ predict(fit.cox))[1]
# 
# ####### Fitting  AFT Model without Clustering #######################
# fit.aft <- survreg(smod ~ Y[,1:D] , dist="lognormal")
# cindex.aft <- rcorr.cens(predict(fit.aft), smod)[1]

### Fitting A Penalized Cox Proportional Hazard's Model
reg.pcox <- cv.glmnet(x = Y, y = smod, family = "cox")
lp <- predict(object =reg.pcox, newx = Y, s= "lambda.min")
cindex.pcox <- survConcordance(smod ~lp)[1]
cindex.pen.cox <<- as.numeric(cindex.pcox)


#### Fitting A AFT Model 
reg <- cv.glmnet(x = Y, y = time, family = "gaussian")
linear.pred <- predict(object =reg, newx = Y, s= "lambda.min")
coeff.pred <- coef(object =reg, newx = Y.tmp, s= "lambda.min")
linear.aft <- predict(object = reg, newx = Y, s = "lambda.min") 
cindex.paft <- survConcordance(smod ~ exp(-linear.aft))[1]
cindex.pen.aft <<- as.numeric(cindex.paft)

#############################################
########### K-means #########################
#############################################
#############################################
gr.km <- kmeans(Y, F, nstart =10)
gr.km.rand <- adjustedRandIndex(c.true,as.factor(gr.km$cluster))






# ########## CoxPH #############################
# fit.cox.km <- coxph(smod ~ Y[,1:D]*strata(as.factor(gr.km$cluster)) , data = as.data.frame(Y))
# 
# ## C-Index
# cindex.km.cox <- survConcordance(smod ~ predict(fit.cox.km))[1]
# ## Brier Score
# fit.coxph <- survfit(fit.cox.km, newdata = as.data.frame(Y[,1:D]))
# ##brier.km.cox <- sbrier(Surv(fit.coxph$time,fit.coxph$n.event), fit.coxph$surv)[[1]]
#                      
# 
# 
# ######## AFT ###################################
# fit.aft.km <- survreg(smod ~ Y[,1:D]*strata(as.factor(gr.km$cluster)) , dist="lognormal")
# cindex.km.aft <- rcorr.cens(predict(fit.aft.km), smod)[1]
# # predict.km.aft <- predict(fit.aft.km,newdata = as.data.frame(Y),type ="quantile",p=seq(.01,.99,by=.01))
# # predict.km.aft <- predict(fit.aft.km,newdata = as.data.frame(Y),type ="response")
# 


### Brier Score
# brier.km.aft <- c(0)
# for ( q in 1:F){
#   ind <- which((gr.km$cluster) == q)
#   time.tmp <- time[ind]
#   censoring.tmp <- censoring[ind]
#   Y.tmp <- Y[ind,]
#   rownames(Y.tmp) <- as.character(c(1:nrow(Y.tmp)))
#   smod.tmp <-  Surv(exp(time.tmp), censoring.tmp) 
#   f1 <- psm(smod.tmp ~ Y.tmp[,1:D]  , dist="lognormal")
#   S1 <- Survival(f1)
#   L = length(ind)
#   mat.tmp <- matrix(NA, nrow = L, ncol = L)
#   for (j in 1:L){
#     mat.tmp[,j] <- S1(exp(time.tmp[j]),f1$linear.predictors)
#   }
#   brier.km.aft[q] <- sbrier(smod.tmp,mat.tmp,exp(time.tmp))[1]
#   
# }
# 



# brier.km.pcox <- c(0)
# cindex.km.pcox <-0
# brier.km.paft <- c(0)
# cindex.km.paft <- 0
# cindex.true.pcox <- c(0)



# ######## Penalized Cox PH with TRue Clustering ###########################################
# linear.pred <- c(0)
# for ( q in 1:F){
#   ind <- which((c.true) == q)
#   time.tmp <- time[ind]
#   censoring.tmp <- censoring[ind]
#   Y.tmp <- Y[ind,1:rel.D]
#   coxreg <- list(0)
#   coxreg$x <- Y.tmp
#   coxreg$time <- exp(time.tmp)
#   coxreg$status <- censoring.tmp
#   # path <- coxpath(data = coxreg)
#   # f.reg <- predict(object = path, data = coxreg, s =5,type =  "coxph", mode = "lambda")
#   # fit.coxregph <-  survfit(f.reg, newdata = as.data.frame(Y.tmp[,1:D]))
#   # brier.km.pcox[q] <- sbrier(Surv(fit.coxregph$time,fit.coxregph$n.event), fit.coxregph$surv)[[1]]
#   # cindex.km.pcox[q]  <-  survConcordance(Surv(coxreg$time,coxreg$status) ~ predict(f.reg))[1]
#   ### see if we can use glmnet
#   reg.pcox <- cv.glmnet(x = Y.tmp, y = Surv(coxreg$time, coxreg$status), family = "cox")
#   linear.pred[ind] <- predict(object =reg.pcox, newx = Y.tmp, s= "lambda.min")
# }
# 
# cindex.true.pcox <<- as.numeric(survConcordance(smod ~ linear.pred)[1])
# 


######## Penalized Cox PH ###########################################
linear.cox <- c(0)
for ( q in 1:F){
ind <- which((gr.km$cluster) == q)
time.tmp <- time[ind]
censoring.tmp <- censoring[ind]
Y.tmp <- Y[ind,]
coxreg <- list(0)
coxreg$x <- Y.tmp
coxreg$time <- exp(time.tmp)
coxreg$status <- censoring.tmp
# path <- coxpath(data = coxreg)
# f.reg <- predict(object = path, data = coxreg, s =5,type =  "coxph", mode = "lambda")
# fit.coxregph <-  survfit(f.reg, newdata = as.data.frame(Y.tmp[,1:D]))
# brier.km.pcox[q] <- sbrier(Surv(fit.coxregph$time,fit.coxregph$n.event), fit.coxregph$surv)[[1]]
# cindex.km.pcox[q]  <-  survConcordance(Surv(coxreg$time,coxreg$status) ~ predict(f.reg))[1]
### see if we can use glmnet
reg.pcox <- cv.glmnet(x = Y.tmp, y = Surv(coxreg$time, coxreg$status), family = "cox")
linear.cox[ind] <- predict(object =reg.pcox, newx = Y.tmp, s= "lambda.min")
}

cindex.km.pcox <- survConcordance(smod ~ linear.cox)[1]


# ####### Penalized AFT with true clustering ############
# linear.aft <- c(0)
# for ( q in 1:F){
#   ind <- which((c.true) == q)
#   L= length(ind)
#   
#   time.tmp <- time[ind]
#   censoring.tmp <- censoring[ind]
#   Y.tmp <- Y[ind,1:rel.D]
#   
#   reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
#   linear.pred <- predict(object =reg, newx = Y.tmp, s= "lambda.min")
#   coeff.pred <- coef(object =reg, newx = Y.tmp, s= "lambda.min")
#   rel.coeff <- coeff.pred[2:(D+1)] 
#   ind.rel <- which(rel.coeff !=0)
#   linear.aft[ind] <- predict(object = reg, newx = Y.tmp, s = "lambda.min") 
# }
# cindex.true.paft <<- as.numeric(survConcordance(smod ~ exp(-linear.aft))[1])
# 



######## Penalized AFT ######################################################

linear.aft <- c(0)
for ( q in 1:F){
ind <- which((gr.km$cluster) == q)
L= length(ind)

time.tmp <- time[ind]
censoring.tmp <- censoring[ind]
Y.tmp <- Y[ind,]
  
reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
linear.pred <- predict(object =reg, newx = Y.tmp, s= "lambda.min")
coeff.pred <- coef(object =reg, newx = Y.tmp, s= "lambda.min")
rel.coeff <- coeff.pred[2:(D+1)] 
ind.rel <- which(rel.coeff !=0)
linear.aft[ind] <- predict(object = reg, newx = Y.tmp, s = "lambda.min") 
}
cindex.km.paft <- survConcordance(smod ~ exp(-linear.aft))[1]

# sum.err <- c(0)
# for ( j in 1:L){
#   sum.err[j] <- (time.tmp[j] - predicted.tmp[j])^2 
# }
# variance = sqrt(median(sum.err))
# 
# ## This is just to get a survival function1 - pnorm((t.trans - lp)/exp(parms))
# aft_survival = function (times = NULL, lp = NULL, parms = variance) 
# {
#   t.trans <- logb(times)
#   names(t.trans) <- format(times)
#   1 - pnorm((t.trans - lp)/(parms))
# }
# 
# mat.reg.tmp <- matrix(NA, nrow = L, ncol = L)
# for (j in 1:L){
#   mat.reg.tmp[,j] <- aft_survival(exp(time.tmp[j]),predicted.tmp)
# }
# brier.km.paft[q] <- sbrier(Surv(time.tmp,censoring.tmp),mat.reg.tmp,exp(time.tmp))[1]
# }



##### Save some Ground truth statistics

gr.km.rand.final <<- gr.km.rand
# cindex.km.cox.final <<- as.numeric(cindex.km.cox)
# cindex.km.aft.final <<-  as.numeric(cindex.km.aft)

#brier.km.cox.final <<- as.numeric(brier.km.cox)
#brier.km.aft.final <<- mean(brier.km.aft)

cindex.km.pcox.final <<- as.numeric(cindex.km.pcox)
cindex.km.paft.final <<- as.numeric(cindex.km.paft)

#brier.km.pcox.final <<- mean(unlist(brier.km.pcox))
#brier.km.paft.final <<- mean(unlist(brier.km.paft))



#################################################################################
##############################################################################
############### Penalized FlexMix #######################################################
################################################################################

data <- data.frame(time, Y)

## The cross validation folds for choosing lambda
fo <- sample(rep(seq(10), length = nrow(data)))
gr.flx <- flexmix(time ~ ., data = data, k = F, cluster = gr.km$cluster, model = FLXMRglmnet(foldid = fo, adaptive= FALSE), control = list(iter.max = 500))
linear.flxmix.recovery <- unlist(predict(gr.flx, newdata = data, aggregate = TRUE))
cindex.flx <- as.numeric(survConcordance(smod ~ exp(-linear.flxmix.recovery))[1])
rand.flx <-  adjustedRandIndex(c.true,as.factor(clusters(gr.flx)))

gr.flx.rand.final <<- rand.flx
cindex.flx.final <<- as.numeric(cindex.flx)



#brier.flx.cox.final <<-  brier.flx.cox

#cindex.cox <<- as.numeric(cindex.cox)

## cindex.aft <<- as.numeric(cindex.aft)


############## Testing Mixture of Factor Analyzers ##########################3

mcfa.fit<- mcfa(Y, g= k, q=2, itmax=250, nkmeans=5, nrandom=5, tol=1.e-3)

########## Seeing if the PCA plot with the corresponding features with releevant features makes sense
randindexMCFA <<- adjustedRandIndex(mcfa.fit$clust, c.true)



#############################################
########### sparse K-means #########################
#############################################
#############################################
km.perm <- KMeansSparseCluster.permute(x = Y, K= k ,wbounds=c(1.5,2:6),nperms=5)
km.out <- KMeansSparseCluster(x = Y, K=k,wbounds=km.perm$bestw)
randindexSKM <<- adjustedRandIndex(km.out[[1]]$Cs, c.true)

###################################################
########### sparse hierarchical clustering #########################
#############################################
#############################################
perm.out <- HierarchicalSparseCluster.permute(x = Y, wbounds=c(1.5,2:6), nperms=5)
# Perform sparse hierarchical clustering
sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists, wbound=perm.out$bestw, method="complete")
randindexSHC <<- adjustedRandIndex(cutree(sparsehc$hc, k = k), c.true)





}

