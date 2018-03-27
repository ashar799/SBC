### This function takes the posterior parameters AND predicts CLUSTER MEMEBERSHIP for the new points
### The posterior parameters being the RECOVERED CLUSTER MEMBERSHIP
#### THEN IT FITS A LOGISTIC/MULTINOMIAL regression FOR LINEAR SEPARABILITY
#### WE CAN ALSO EXPLORE NON-LINEAR Classifiers as k-SVM (RBF kernel)
### Check to see if this gives differences to state of the art

### POSSIBLE APPROACHES TO GET POINT ESTIMATES #####



### POSSIBILITY 1
### Use MPEAR APPROACH to get POINT ESTIMATES #############################################################################
psm2 <- comp.psm(t(c.matrix.new))
mpear2 <- maxpear(psm2)
### Generally the MPEAR output needs post-processing
c.sbc.new.mpear <- mpear2$cl




### POSSIBILITY 2
### Use Logistic regression to get the labels
reg <- cv.glmnet(x = Y, y = as.factor(c.sbc), family = "binomial")
reg.new <- predict(object = reg, newx = Y.new, s = "lambda.min", type="class") 
c.sbc.new.log <- as.numeric(reg.new)


### POSSIBILITY 3
### Use k-Nearest neighbour
label.train <- as.factor(c.sbc)
### One has to to tune the k-NN classifier for k ###
fitControl <- trainControl(method = "repeatedcv", number = 5,repeats = 5)
### Tune the parameter k 
knnFit <- train(x = Y, y = label.train, method =  "knn", trControl = fitControl, tuneLength = 5)
c.sbc.new.knn <- predict(knnFit,newdata = Y.new )
#### Choose that configuration which has the highest difference in survival curves ####################

### POSSIBILITY 4
### Use k-svm 

library(e1071)
library(kernlab)
df <- data.frame(cbind(c.sbc,Y))
sig <- sigest(c.sbc ~.,data = df)
range.sig <-  seq(from = sig[1],to = sig[3], length.out = 25)
obj.base <- tune(svm, train.x = Y, train.y = factor(c.sbc), ranges = list(gamma = 2^(-10:14), cost = 2^(-10:14)), tunecontrol = tune.control(sampling = "cross"))

sigma.final <- obj.base$best.parameters$gamma
cost.final <- obj.base$best.parameters$cost

# Training the Classifier
ksvm.verk = ksvm(Y, factor(c.sbc), kernel = "rbfdot",kpar = list(sigma =sigma.final),C= cost.final, prob.model =TRUE)
## Predicting
labels.ksvm <- predict(ksvm.verk ,Y.new, type = "response")
pred.logrank.rbf <- survdiff(smod.new ~ labels.ksvm)
c.sbc.new.rbf <- labels.ksvm


### POSSIBILITY 5
### Use one of the actual cluster assignments

c.matrix.new.sbc <- matrix(0, nrow = N.new, ncol = Nps + 4)
c.matrix.new.sbc[,1:Nps] <- c.matrix.new[,1:Nps]
c.matrix.new.sbc[,(Nps+1)] <- c.sbc.new.mpear
c.matrix.new.sbc[,(Nps+2)] <- c.sbc.new.log
c.matrix.new.sbc[,(Nps+3)] <- c.sbc.new.knn 
c.matrix.new.sbc[,(Nps+4)] <- c.sbc.new.rbf


Nps.mod <- Nps +4
lr <- c(0)
for (j in 1:Nps.mod){
  lr[j] <-  1 - pchisq(unlist(survdiff(smod.new ~ c.matrix.new.sbc[,j]))$chisq,df = nlevels(as.factor(c.matrix.new.sbc[,i])) -1 )
}

c.sbc.new <- c.matrix.new.sbc[,which.min(lr)]















