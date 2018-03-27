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

### The Vijver Classification Train ###
recov.logrank.vijver.final <- c(0)
recovCIndex.vv.pcox.final <- c(0)


### The Vijver Classification Test ###
pred.logrank.vijver.final <- c(0)
predCIndex.vv.pcox.final <- c(0)

## The VIJVER Classification train followed by k-nearest neighbour
pred.logrank.vv.kk.final <- c(0)
predCIndex.vv.kk.pcox.final <- c(0)


### Using K-means and possibly k-nearest neighbour ###
recov.logrank.kk.final <- c(0)
pred.logrank.kk.final <- c(0)

recovCIndex.kk.pcox.final <- c(0)
predCIndex.kk.pcox.final <- c(0)



#### A L1 penalized Cox model based on all genes #####
recovCIndex.NA.pcox.final <- c(0)
predCIndex.NA.pcox.final <- c(0)


### Using L1 penalized on SBC genes #####
recovCIndex.NAS.pcox.final <- c(0)
predCIndex.NAS.pcox.final <- c(0)



icount =1


for ( u in 1:5) {
     for ( v in 1:5){
        if (!((u == 3) & (v ==5))){
   load(paste('/home/bit/ashar/ExpressionSets/CROSS_VALIDATION/Breast/','repeat',u,'split',v,'.RData',sep = ""))
    
    ### Get the Summaries ####
    
    ### For the SBC method ######
    recovCIndex.sbc.final[icount] <- mean(recovCIndex.sbc)
    predCIndex.sbc.final[icount] <- max(predCIndex.sbc)
    
    recov.logrank.sbc.final[icount] <- recov.logrank.sbc
    pred.logrank.sbc.final[icount] <- unlist(survdiff(smod.new ~ c.sbc.new))$chisq
    
    
    ### For the PC method #####
    recovCIndex.PC.final[icount] <- recovCIndex.PC
    predCIndex.PC.final[icount] <- predCIndex.PC
    
    ### The Vijver Classification Train ##
    recov.logrank.vijver.final[icount] <- recov.logrank.vijver
    recovCIndex.vv.pcox.final[icount] <-  recovCIndex.vv.pcox
    
    
    ### The Vijver Classification Test ###
    pred.logrank.vijver.final[icount] <- pred.logrank.vijver
    predCIndex.vv.pcox.final[icount] <- predCIndex.vv.pcox
    
    ## The VIJVER Classification train followed by k-nearest neighbour
    pred.logrank.vv.kk.final[icount] <- pred.logrank.vv.kk
    predCIndex.vv.kk.pcox.final[icount] <- predCIndex.vv.kk.pcox
    
    
    ### Using K-means and possibly k-nearest neighbour ###
    recov.logrank.kk.final[icount] <- recov.logrank.kk
    pred.logrank.kk.final[icount] <- pred.logrank.kk
    
    recovCIndex.kk.pcox.final[icount] <- recovCIndex.kk.pcox
    predCIndex.kk.pcox.final[icount] <- predCIndex.kk.pcox
    
    
    
    recovCIndex.NA.pcox.final[icount] <- recovCIndex.NA.pcox
    predCIndex.NA.pcox.final[icount] <- predCIndex.NA.pcox
    
    
    
    recovCIndex.NAS.pcox.final[icount] <- recovCIndex.NAS.pcox
    predCIndex.NAS.pcox.final[icount] <- predCIndex.NAS.pcox
    
    
    icount <- icount +1
    }
   
}
}

source('multiplot.R')

##### Model Fitting ####
cindex.recov <- cbind(recovCIndex.sbc.final,recovCIndex.PC.final,recovCIndex.vv.pcox.final,recovCIndex.kk.pcox.final, recovCIndex.NA.pcox.final, recovCIndex.NAS.pcox.final)
colnames(cindex.recov) <- c("SBC","PrComp","VV","HC","ALL.pCOX","SBC.pCOX")
cind.recov <- melt(cindex.recov)
p1 <- ggplot(data = as.data.frame(cind.recov)) + geom_boxplot(aes(y = cind.recov$value, x= factor(as.factor(cind.recov$X2), levels = colnames(cindex.recov)), fill = (cind.recov$X2))) + ggtitle("Training C-Index \n Breast Cancer \n  5 * 5 cross-validation") + labs(y = "C-Index", x = "Models") + theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())

#### Model Prediction
cindex.pred <- cbind(predCIndex.sbc.final, predCIndex.PC.final, predCIndex.vv.pcox.final, predCIndex.vv.kk.pcox.final, predCIndex.kk.pcox.final, predCIndex.NA.pcox.final, predCIndex.NAS.pcox.final)
colnames(cindex.pred) <- c("SBC","PrComp","VV","VV+kNN","HC+kNN", "ALL.pCOX","SBC.pCOX")
cind.pred <- melt(cindex.pred)
p2 <- ggplot(data = as.data.frame(cind.pred)) + geom_boxplot(aes(y = cind.pred$value, x= factor(as.factor(cind.pred$X2), levels = colnames(cindex.pred)), fill = (cind.pred$X2))) + ggtitle("Testing C-Index \n Breast Cancer \n  5 * 5 cross-validation") + labs(y = "C-Index", x = "Models") + theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())

#### 

###### Plotting of the Log-Rank statistic for the Cross Validation #####
##### Model Fitting ####
logrank.recov <- cbind(recov.logrank.sbc.final, recov.logrank.vijver.final, recov.logrank.kk.final )
colnames(logrank.recov) <- c("SBC","VV","HC")
lg.recov <- melt(logrank.recov)
p3 <- ggplot(data = as.data.frame(lg.recov)) + geom_boxplot(aes(y = lg.recov$value, x= factor(as.factor(lg.recov$X2), levels = colnames(logrank.recov)), fill = (lg.recov$X2))) + ggtitle("Training Chi-squared-statistic \n Breast Cancer \n  5 * 5 cross-validation") + labs(y = expression(paste(chi^2,"-statistic")), x = "Models") + theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())


pred.logrank.sbc.final <- pred.logrank.sbc.final[-c(11,22)]
logrank.pred <- cbind(pred.logrank.sbc.final, pred.logrank.vijver.final, pred.logrank.vv.kk.final, pred.logrank.kk.final )
colnames(logrank.pred) <- c("SBC","VV","VV+kNN","HC")
lg.pred <- melt(logrank.pred)
p4 <- ggplot(data = as.data.frame(lg.pred)) + geom_boxplot(aes(y = lg.pred$value, x= factor(as.factor(lg.pred$X2), levels = colnames(logrank.pred)), fill = (lg.pred$X2))) + ggtitle("Testing Chi-squared -statistic \n Breast Cancer \n  5 * 5 cross-validation") + labs(y = expression(paste(chi^2,"-statistic")), x = "Models") + theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())

pdf('BreastCancerCross.pdf')
p1
p2
p3
p4
dev.off()


