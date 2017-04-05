### This function compares feature selection and C-Index recovery and prediction ###
### Let's see what comes up



ComparisionPReMiuM = function(){
  
  smod <-  Surv(exp(time), censoring)
  ### Prediction ####
  smod.new <-  Surv(exp(time.new), censoring.new)
  
  ####### PREMIUM ############################
  direc <- as.character('/home/bit/ashar/ownCloud/Research/PreMiuM')
  setwd(direc)
  data.premium <- list(NULL)
  data.premium$inputData <- as.data.frame(cbind(Y,time))
  data.premium$covNames <- paste0(rep("Variable",D),as.character(c(1:D)))
  names(data.premium$inputData) <- c(data.premium$covNames,"outcome")
  data.premium$xModel <- "Normal"
  data.premium$yModel <- "Normal"
  data.premium$ncovariates <- D
  data.premium$outputData <- as.data.frame(Y.new)
  names(data.premium$outputData) <- paste0(rep("Variable",D),as.character(c(1:D)))


  ####### Run the Profile Regression

  runInfoObj <-profRegr(yModel=data.premium$yModel, xModel= data.premium$xModel, nSweeps= 200, nBurn= 100, data= data.premium$inputData, output= "outcome", nFilter= 5 ,covNames= data.premium$covNames,nClusInit= k,reportBurnIn=FALSE, predict = data.premium$outputData)
  dissimObj<- calcDissimilarityMatrix(runInfoObj)
  clusObj<-calcOptimalClustering(dissimObj)

  
  ###### Use this clustering to fit AFT survival models
  linear.premium <- c(0)
  linear.premium.pred <- c(0)
  for ( q in 1:F){
    ind <- which(clusObj$clustering == q)
    ind.new <- which(clusObj$clusteringPred == q)
   
    Y.tmp <- Y[ind,]
    time.tmp <- time[ind]

   
    Y.tmp.new <- Y.new[ind.new,]
    reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
    linear.premium[ind] <- predict(object = reg, newx = Y.tmp, s = "lambda.min") 
    linear.premium.pred[ind.new] <- predict(object = reg, newx = Y.tmp.new, s = "lambda.min") 
    
    }
  recovCIndex.premium.aft  <<- as.numeric(survConcordance(smod ~ exp(-linear.premium))[1])
  predCIndex.premium.aft  <<- as.numeric(survConcordance(smod.new ~ exp(-linear.premium.pred))[1])
  
  
  
  
  
  ######## Rand Indices ##############################
  recovRandIndex.premium <<- adjustedRandIndex(c.true,clusObj$clustering)
  predRandIndex.premium <<-  adjustedRandIndex(c.true.new,clusObj$clusteringPred)

  ########## Making predictions #####################
  riskProfileObj <- calcAvgRiskAndProfile(clusObj)

  output_predictions <- calcPredictions(riskProfileObj)

  predCIndex.premium <<-  as.numeric(survConcordance(smod.new ~ exp(-output_predictions$predictedY))[1])

  
  }
