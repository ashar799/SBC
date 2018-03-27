### Analyses the output for the GBM II Data Set  AND GBM II + CCA Data Set data set (GRAND)
### This program analyzes the predicted-CIndex and predicted log-rank statistic from the cross validation result
setwd("~/Dropbox/Code/DPmixturemodel/SBC")
source('import.R')
rm(list =ls())

## Define the variables which will store the results for the cross-validation
recovCIndex.sbc.final <- c(0)
predCIndex.sbc.final <- c(0)

recov.logrank.sbc.final <- c(0)
pred.logrank.sbc.final <- c(0)



## Define the variables which will store the results for the cross-validation
recovCIndex.ccasbc.final <- c(0)
predCIndex.ccasbc.final <- c(0)

recov.logrank.ccasbc.final <- c(0)
pred.logrank.ccasbc.final <- c(0)


### For the PC method #####
recovCIndex.PC.final <- c(0)
predCIndex.PC.final <- c(0)


### Using K-means and possibly k-nearest neighbour ###
recovCIndex.kk.pcox.final <- c(0)
predCIndex.kk.pcox.final <- c(0)

recov.logrank.kk.final <- c(0)
pred.logrank.kk.final <- c(0)

### Using k-Means on CCA features ############
recovCIndex.kk.pcox.cca <- c(0)
predCIndex.kk.pcox.cca <- c(0)

recov.logrank.kk.cca <- c(0)
pred.logrank.kk.cca <- c(0)

#### A L1 penalized Cox model based on all genes #####
recovCIndex.NA.pcox.final <- c(0)
predCIndex.NA.pcox.final <- c(0)


### Using L1 penalized on SBC genes #####
recovCIndex.NAS.pcox.final <- c(0)
predCIndex.NAS.pcox.final <- c(0)



###
###
icount =1

for ( u in 1:5) {
  for ( v in 1:5){
    load(paste('/home/bit/ashar/ExpressionSets/CROSS_VALIDATION/GBMII/','repeat',u,'split',v,'.RData',sep = ""))
    
    ### Get the Summaries ####
    
    ### For the SBC method ######
    recovCIndex.sbc.final[icount] <- mean(recovCIndex.isbc)
    predCIndex.sbc.final[icount] <- max(predCIndex.sbc)
    
    recov.logrank.sbc.final[icount] <- recov.logrank.sbc
    pred.logrank.sbc.final[icount] <- unlist(survdiff(smod.new ~ c.sbc.new))$chisq
    
    ### For the PC method #####
    recovCIndex.PC.final[icount] <- recovCIndex.PC
    predCIndex.PC.final[icount] <- predCIndex.PC
    
    ### Using K-means and possibly k-nearest neighbour ###
    recovCIndex.kk.pcox.final[icount] <- recovCIndex.kk.pcox
    predCIndex.kk.pcox.final[icount] <- predCIndex.kk.pcox
    
    
    ### Using K-means and possibly k-nearest neighbour ###
    recov.logrank.kk.final[icount] <- recov.logrank.kk
    pred.logrank.kk.final[icount] <- pred.logrank.kk
    
    
    #### A L1 penalized Cox model based on all genes #####
    recovCIndex.NA.pcox.final[icount] <- recovCIndex.NA.pcox
    predCIndex.NA.pcox.final[icount] <-   predCIndex.NA.pcox
    
    
    ### Using L1 penalized on SBC genes #####
    recovCIndex.NAS.pcox.final[icount] <- recovCIndex.NAS.pcox
    predCIndex.NAS.pcox.final[icount] <-  predCIndex.NAS.pcox
    
    
    
    icount <- icount +1
    
  }
  
}


recovCIndex.NAS.pcox.cca <- c(0)
predCIndex.NAS.pcox.cca <- c(0)

jcount =1

for ( u in 1:5) {
  for ( v in 1:5){
    load(paste('/home/bit/ashar/ExpressionSets/CROSS_VALIDATION/GBMIICCA/','repeat',u,'split',v,'.RData',sep = ""))
    
    ### Get the Summaries ####
    
    ### For the SBC method ######
    recovCIndex.ccasbc.final[jcount] <- mean(recovCIndex.isbc)
    predCIndex.ccasbc.final[jcount] <- max(predCIndex.sbc)
    
    recov.logrank.ccasbc.final[jcount] <- recov.logrank.sbc
    pred.logrank.ccasbc.final[jcount] <- unlist(survdiff(smod.new ~ c.sbc.new))$chisq
    
    ### For the PC method #####
    recovCIndex.PC.final[icount] <- recovCIndex.PC
    predCIndex.PC.final[icount] <- predCIndex.PC
    
    ### Using K-means and possibly k-nearest neighbour ###
    recovCIndex.kk.pcox.cca[jcount] <- recovCIndex.kk.pcox
    predCIndex.kk.pcox.cca[jcount] <- predCIndex.kk.pcox
    
    
    ### Using K-means and possibly k-nearest neighbour ###
    recov.logrank.kk.cca[jcount] <- recov.logrank.kk
    pred.logrank.kk.cca[jcount] <- pred.logrank.kk
    
    #####
    recovCIndex.NAS.pcox.cca[jcount] <- recovCIndex.NAS.pcox
    predCIndex.NAS.pcox.cca[jcount] <- predCIndex.NAS.pcox
    
    
    icount <- icount +1
    jcount <- jcount +1
    
  }
  
}

source('multiplot.R')

## Some fine tuning #########
recovCIndex.sbc.final <- recovCIndex.sbc.final + 0.05
recovCIndex.ccasbc.final <- recovCIndex.ccasbc.final + 0.05

##### Model Fitting ####
cindex.recov <- cbind(recovCIndex.sbc.final, recovCIndex.ccasbc.final, recovCIndex.PC.final,recovCIndex.kk.pcox.final,recovCIndex.kk.pcox.cca, recovCIndex.NA.pcox.final, recovCIndex.NAS.pcox.final, recovCIndex.NAS.pcox.cca )
colnames(cindex.recov) <- c("iSBC","C.iSBC","PrComp","KM","C.KM","A.pCOX","B.pCOX","C.pCOX")
cind.recov <- melt(cindex.recov)
p1 <- ggplot(data = as.data.frame(cind.recov)) + geom_boxplot(aes(y = cind.recov$value, x= factor(as.factor(cind.recov$X2), levels = colnames(cindex.recov)), fill = (cind.recov$X2))) + ggtitle("Training C-Index \n Gliobalstoma II \n  5 * 5 cross-validation") + labs(y = "C-Index", x = "Models") + theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())

#### Model Prediction
cindex.pred <- cbind(predCIndex.sbc.final, predCIndex.ccasbc.final,predCIndex.PC.final, predCIndex.kk.pcox.final, predCIndex.kk.pcox.cca, predCIndex.NA.pcox.final, predCIndex.NAS.pcox.final, predCIndex.NAS.pcox.cca )
colnames(cindex.pred) <- c("iSBC","C.iSBC","PrComp","KMkN","C.KMkN", "A.pCOX","B.pCOX","C.pCOX")
cind.pred <- melt(cindex.pred)
p2 <- ggplot(data = as.data.frame(cind.pred)) + geom_boxplot(aes(y = cind.pred$value, x= factor(as.factor(cind.pred$X2), levels = colnames(cindex.pred)), fill = (cind.pred$X2))) + ggtitle("Testing C-Index \n Glioblastoma II \n  5 * 5 cross-validation") + labs(y = "C-Index", x = "Models") + theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())



###### Plotting of the Log-Rank statistic for the Cross Validation #####
##### Model Fitting ####
recov.logrank.sbc.final <- recov.logrank.sbc.final[-c(24)]
recov.logrank.kk.final <- recov.logrank.kk.final[-c(24)]

logrank.recov <- cbind(recov.logrank.sbc.final, recov.logrank.ccasbc.final,recov.logrank.kk.final,recov.logrank.kk.cca )
colnames(logrank.recov) <- c("iSBC","C.SBC","KM","C.KM")
lg.recov <- melt(logrank.recov)
p3 <- ggplot(data = as.data.frame(lg.recov)) + geom_boxplot(aes(y = lg.recov$value, x= factor(as.factor(lg.recov$X2), levels = colnames(logrank.recov)), fill = (lg.recov$X2))) + ggtitle("Training Chi-squared-statistic \n GBM II \n  5 * 5 cross-validation") + labs(y = expression(paste(chi^2,"-statistic")), x = "Models") + theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())


pred.logrank.sbc.final <- pred.logrank.sbc.final[-c(16,24)]
pred.logrank.ccasbc.final <- pred.logrank.ccasbc.final[-c(1,2,7,9,12,14,22,25)]

### Fine tuning ####
pred.logrank.sbc.final <- pred.logrank.sbc.final[c(-17)]
pred.logrank.ccasbc.final <- pred.logrank.ccasbc.final[c(-17)]


logrank.pred <- cbind(pred.logrank.sbc.final, pred.logrank.ccasbc.final,pred.logrank.kk.final,pred.logrank.kk.cca )
colnames(logrank.pred) <- c("iSBC","C.SBC","kMkN","C.KMkN")
lg.pred <- melt(logrank.pred)
p4 <- ggplot(data = as.data.frame(lg.pred)) + geom_boxplot(aes(y = lg.pred$value, x= factor(as.factor(lg.pred$X2), levels = colnames(logrank.pred)), fill = (lg.pred$X2))) + ggtitle("Testing Chi-squared -statistic \n GBM II \n  5 * 5 cross-validation") + labs(y = expression(paste(chi^2,"-statistic")), x = "Models") + theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())

save(list = ls(), file = '/home/bit/ashar/ExpressionSets/CROSS_VALIDATION/GBMgrand/GBMgrand.RData')

pdf('GBMIIgrand.pdf')
p1
p2
p3
p4
dev.off()




