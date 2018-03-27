### This function augments the situations in MAIN_SIMULATION
### After training

for ( i in 1:Nps){
cluster1 <- rMVN(5000, mean=mu.list[[i]][2,], Q=S.list[[i]][2,1:D,1:D])$x
cluster2 <- Y.new[1:250,1:D]
cluster3 <- Y.new[251:500,1:D]
cluster.all <- rbind(cluster1,cluster2, cluster3)
label.all <- c(rep(1,5000),rep(2,250), rep(3,250))
pc <- prcomp(cluster.all)
pc.pred <- predict(pc,newdata = cluster.all)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(label.all))) + ggtitle("Class Separability") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 
p1

}

hmcols<-colorRampPalette(c("white","black"))(128)
heatmap.2(t(as.matrix(weights[1:2,])),dendrogram="none", margins=c(6,10), col =hmcols, main = "Posterior prob. \n for assignment \n ", cexCol = 0.85, cexRow = 0.7, Rowv = FALSE, trace = 'none')

#### Let's see if we build classifiers based on the present class labels of training data can we obtain good predictions

### k-NN




### Let's use a SVM
library(e1071)
library(kernlab)
datasvm <- data.frame(Y)
label.test <- as.factor(c.sbc)
fsig =  sigest(c.sbc ~ ., data = datasvm,  frac = 0.50)
range.sig =  seq(fsig[1],fsig[3],(fsig[3]- fsig[1])/10)
obj.base <- tune(svm, train.x = Y, train.y = factor(c.sbc), ranges = list(gamma = (range.sig)^-1, cost = 2^(-10:14)), tunecontrol = tune.control(sampling = "cross"))

sigma.final <- (obj.base$best.parameters$gamma) ^-1
cost.final <- obj.base$best.parameters$cost
## Training
ksv = ksvm(c.sbc ~ ., data = datasvm, kernel = "rbfdot",kpar = list(sigma =sigma.final),C= cost.final, prob.model =TRUE)
## Predicting
ksv.predict = predict(ksv,Y.new, type = "response")
#
