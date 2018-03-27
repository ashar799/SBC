### This function takes the posterior parameters AND predicts CLUSTER MEMEBERSHIP for the new points
### The posterior parameters being the RECOVERED CLUSTER MEMBERSHIP
#### THEN IT FITS A LOGISTIC/MULTINOMIAL regression FOR LINEAR SEPARABILITY
#### WE CAN ALSO EXPLORE NON-LINEAR Classifiers as k-SVM or k-DA
### Check to see if this gives differences to state of the art

Y <- cbind(Y1,Y2)
Y.new <- cbind(Y1.test,Y2.test)
label.train <- c.sbc

#### Fitting a penalized logistic regression Model to make predictions ####
reg <- cv.glmnet(x = Y, y = as.factor(label.train), family = "multinomial")
reg.new <- predict(object = reg, newx = Y.new, s = "lambda.min", type="class") 
label.test <- as.numeric(reg.new)
survdiff(smod.new ~ label.test)

### Fit a penalized Cox regression model
linear.sbclog.prediction <- c(0)
for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.aft <- cv.glmnet(x = Y[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox", alpha = 0.0)
  linear.sbclog.prediction[ind.new] <- predict(object =reg.aft, newx = Y.new[ind.new,], s= "lambda.min")
}
predCIndex.sbclog.pcox <- as.numeric(survConcordance(smod.new ~ linear.sbclog.prediction)[1])






#### Using k-NN to make predictions ####
#### Use L2 penalized Cox regression ####
label.test <- knn(train = Y, test = Y.new, cl = label.train, k = F)
survdiff(smod.new ~ label.test)

linear.sbcknn.prediction <- c(0)
for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.aft <- cv.glmnet(x = Y[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox", alpha = 0.0)
  linear.sbcknn.prediction[ind.new] <- predict(object =reg.aft, newx = Y.new[ind.new,], s= "lambda.min")
}
predCIndex.sbcknn.pcox <- as.numeric(survConcordance(smod.new ~ linear.sbcknn.prediction)[1])




### Fitting a kSVM #####
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
survdiff(smod.new ~ labels.ksvm)
label.test <- labels.ksvm


### Now fitting penalized Cox regression in R ###
linear.sbcrbf.prediction <- c(0)



for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.aft <- cv.glmnet(x = Y[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox", alpha = 0.5)
  linear.sbcrbf.prediction[ind.new] <- predict(object =reg.aft, newx = Y.new[ind.new,], s= "lambda.min")
}
predCIndex.sbcrbf.pcox <- as.numeric(survConcordance(smod.new ~ linear.sbcrbf.prediction)[1])

##### Now fitting penalized AFT ###
linear.sbcrbf.prediction <- c(0)   
for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.aft <- cv.glmnet(x = Y[ind,], y = time[ind], family = "gaussian", alpha = 0.5)
  linear.sbcrbf.prediction[ind.new] <- predict(object =reg.aft, newx = Y.new[ind.new,], s= "lambda.min")
}

predCIndex.sbcrbf.aft <- as.numeric(survConcordance(smod.new ~ linear.sbcrbf.prediction)[1])




##### Choose that cluster prediction method that gives significant differences in Survival curves
c.sbc.new <- label.test






