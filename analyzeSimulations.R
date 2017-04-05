## Analyze the 5 simulations
rm(list =ls())
setwd('/home/bit/ashar/ownCloud/DPMM_RESULTS/DPMMSIMULATIONS/Dimension/D30Over10oise50')

### Rand-Indices
rand.km <- rep(0,5)
rand.flx <- rep(0,5)
rand.sbc <- rep(0,5)

########
rand.SHC <- rep(0,5)
rand.SKM <- rep(0,5)


cindex.sbc <- rep(0,5)
cindex.flx <- rep(0,5)
cindex.paft <- rep(0,5)
cindex.pcox <- rep(0,5)
cindex.kpaft <- rep(0,5)
cindex.kpcox <- rep(0,5)


for (i in 1:5){
  filename <- paste("Sim",i,".RData",sep = "")
  load(filename)
  rand.km[i] <- gr.km.rand.final
  rand.flx[i] <- gr.flx.rand.final
  rand.SHC[i] <- randindexSHC
  rand.SKM[i] <- randindexSHC
  
  ra <- c(0)
  for (u in 1:Nps){
  ra[u] <- adjustedRandIndex(c.list[[u]],c.true)
  }
  rand.sbc <- max(ra)
  
  cindex.flx[i] <- cindex.flx.final 
  cindex.paft[i] <- cindex.pen.aft
  cindex.pcox[i] <- cindex.pen.cox
  
  cindex.kpaft[i] <- cindex.km.paft.final
  cindex.kpcox[i] <- cindex.km.pcox.final
  
  cindex.sbc[i] <- max(cindex.final)
  
  
}

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


