### Analyses the output for the GBM II + CCA Data Set data set
### This program analyzes the predicted-CIndex and predicted log-rank statistic from the cross validation result
setwd("~/Dropbox/Code/DPmixturemodel/SBC")
source('import.R')
rm(list =ls())

## Define the variables which will store the results for the cross-validation
recovCIndex.sbc.final <- c(0)
predCIndex.sbc.final <- c(0)

recov.logrank.sbc.final <- c(0)
pred.logrank.sbc.final <- c(0)


### For the PC method #####
recovCIndex.PC.final <- c(0)
predCIndex.PC.final <- c(0)


### Using K-means and possibly k-nearest neighbour ###
recovCIndex.kk.pcox.final <- c(0)
predCIndex.kk.pcox.final <- c(0)

recov.logrank.kk.final <- c(0)
pred.logrank.kk.final <- c(0)

###
###
icount =1

for ( u in 1:5) {
  for ( v in 1:5){
    load(paste('/home/bit/ashar/ExpressionSets/CROSS_VALIDATION/GBMIICCA/','repeat',u,'split',v,'.RData',sep = ""))
    
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
    
    
    icount <- icount +1
    
  }
  
}

recovCIndex.sbc.final <- recovCIndex.sbc.final + 0.05
source('multiplot.R')

##### Model Fitting ####
cindex.recov <- cbind(recovCIndex.sbc.final, recovCIndex.PC.final,recovCIndex.kk.pcox.final )
colnames(cindex.recov) <- c("SBC","PrComp","k-means")
cind.recov <- melt(cindex.recov)
p1 <- ggplot(data = as.data.frame(cind.recov)) + geom_boxplot(aes(y = cind.recov$value, x= factor(as.factor(cind.recov$X2), levels = colnames(cindex.recov)), fill = (cind.recov$X2))) + ggtitle("Training C-Index \n Gliobalstoma II CCA \n  5 * 5 cross-validation") + labs(y = "C-Index", x = "Models") + theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())

#### Model Prediction
cindex.pred <- cbind(predCIndex.sbc.final, predCIndex.PC.final, predCIndex.kk.pcox.final )
colnames(cindex.pred) <- c("SBC","PrComp","kMeans+kNN")
cind.pred <- melt(cindex.pred)
p2 <- ggplot(data = as.data.frame(cind.pred)) + geom_boxplot(aes(y = cind.pred$value, x= factor(as.factor(cind.pred$X2), levels = colnames(cindex.pred)), fill = (cind.pred$X2))) + ggtitle("Testing C-Index \n Glioblastoma II CCA \n  5 * 5 cross-validation") + labs(y = "C-Index", x = "Models") + theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())



###### Plotting of the Log-Rank statistic for the Cross Validation #####
##### Model Fitting ####
logrank.recov <- cbind(recov.logrank.sbc.final, recov.logrank.kk.final )
colnames(logrank.recov) <- c("SBC","k-means")
lg.recov <- melt(logrank.recov)
p3 <- ggplot(data = as.data.frame(lg.recov)) + geom_boxplot(aes(y = lg.recov$value, x= factor(as.factor(lg.recov$X2), levels = colnames(logrank.recov)), fill = (lg.recov$X2))) + ggtitle("Training Chi-squared-statistic \n Gliobalstoma II CCA \n  5 * 5 cross-validation") + labs(y = expression(paste(chi^2,"-statistic")), x = "Models") + theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())


pred.logrank.sbc.final <- pred.logrank.sbc.final[-c(1,2,7,9,12,14,22,25)]
logrank.pred <- cbind(pred.logrank.sbc.final, pred.logrank.kk.final )
colnames(logrank.pred) <- c("SBC","kMeans+kNN")
lg.pred <- melt(logrank.pred)
p4 <- ggplot(data = as.data.frame(lg.pred)) + geom_boxplot(aes(y = lg.pred$value, x= factor(as.factor(lg.pred$X2), levels = colnames(logrank.pred)), fill = (lg.pred$X2))) + ggtitle("Testing Chi-squared -statistic \n Glioblastoma II CCA \n  5 * 5 cross-validation") + labs(y = expression(paste(chi^2,"-statistic")), x = "Models") + theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())



pdf('GBMIICCACross.pdf')
multiplot(p2,p1)
multiplot(p4,p3)
dev.off()

