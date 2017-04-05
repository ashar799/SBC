#### This file analyses the feature Importance for each feature ######################
rm(list =ls())
setwd('/home/bit/ashar/ownCloud/DPMM_RESULTS/DPMMSIMULATIONS/FeatureImportance')

### Rand-Indices
auc.sbc10 <- matrix(0, nrow = 24, ncol = 4)
auc.flx10 <- c(0)
auc.sbc30 <- matrix(0, nrow = 24, ncol = 3)
auc.flx30 <- c(0)
auc.sbc50 <- matrix(0, nrow = 24, ncol = 2)
auc.flx50 <- c(0)

setwd('/home/bit/ashar/ownCloud/DPMM_RESULTS/DPMMSIMULATIONS/FeatureImportance/D=10,LowNoise')
for (b in 1:4){
  filename <- paste("Sim",b,".RData",sep = "")
  load(filename)
  auc.sbc10[1:24,b] <- auc.sbc 
  auc.flx10[b] <- auc.flx
  }

setwd('/home/bit/ashar/ownCloud/DPMM_RESULTS/DPMMSIMULATIONS/FeatureImportance/D=30,LowNoise')
for (b in 1:3){
  filename <- paste("Sim",b,".RData",sep = "")
  load(filename)
  auc.sbc30[1:24,b] <- auc.sbc 
  auc.flx30[b] <- auc.flx
}

setwd('/home/bit/ashar/ownCloud/DPMM_RESULTS/DPMMSIMULATIONS/FeatureImportance/D=50,LowNoise')
for (b in 1:2){
  filename <- paste("Sim",b,".RData",sep = "")
  load(filename)
  auc.sbc50[1:24,b] <- auc.sbc
  auc.flx50[b] <- auc.flx
}



auc.sbc10.final <- apply(auc.sbc10,2,mean)
auc.flx10 

auc.sbc30.final <- apply(auc.sbc30,2,mean)
auc.flx30

auc.sbc50.final <- apply(auc.sbc50,2,mean)
auc.flx50

aflx <- cbind(auc.flx10, auc.flx30, auc.flx50)
asbc <- cbind(auc.sbc10.final, auc.sbc30.final, auc.sbc50.final)
colnames(aflx) <- c("D=10","D=30","D=50")
rownames(aflx) <- rep("FLxmix",5)
colnames(asbc) <- c("D=10","D=30","D=50")
rownames(asbc) <- rep("SBC",5)

combo <- rbind(melt(aflx),melt(asbc))
p1 <- ggplot(data = as.data.frame(combo)) + geom_boxplot(aes(y = combo$value, x= factor(as.factor(combo$X2)), fill = combo$X1)) + ggtitle("Low Noise \n 5 Simulations AUC \n Feature Importance") + labs(y = "AUC", x = "Dimensions") + theme(legend.title = element_blank())

##### Model Fitting #####################################################
randindex.comp <- cbind(rand.sbc, rand.km, rand.SKM, rand.SHC, rand.flx )
colnames(randindex.comp) <- c("SBC","kmeans","SparsKM","SparsHC","FLXmix")
rind <- melt(randindex.comp)
p1 <- ggplot(data = as.data.frame(rind)) + geom_boxplot(aes(y = rind$value, x= factor(as.factor(rind$X2), levels = colnames(randindex.comp)), fill = (rind$X2))) + ggtitle("Simulation Adj. Rand-Index") + labs(y = "Adjusted Rand-Index", x = "Models") + theme(legend.title = element_blank())

##### Model Fitting ####
cindex.comp <- cbind(cindex.sbc, cindex.flx, cindex.kpaft, cindex.kpcox, cindex.pcox, cindex.paft )
colnames(cindex.comp) <- c("SBC","FLXmix","K-PAFT","K-PCOX","N-PCOX","N-PAFT")
cind <- melt(cindex.comp)
p2 <- ggplot(data = as.data.frame(cind)) + geom_boxplot(aes(y = cind$value, x= factor(as.factor(cind$X2), levels = colnames(cindex.comp)), fill = (cind$X2))) + ggtitle("Simulation C-Index") + labs(y = "C-Index", x = "Models") + theme(legend.title = element_blank())

pdf("Model.Recovery.LESSNOISE.pdf")
p1
p2
dev.off()




