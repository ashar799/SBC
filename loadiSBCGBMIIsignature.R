loadiSBCGBMIIsignature = function(Y1.pre.train, Y2.pre.train,time, censoring){

#################################################################################################################################################################
############################ PRE - FILTERED FEATURE SELECTION ####################################################################################################
  ######## Prefiltering of the Genes ############################### ###########################
  ######### Using Univariate Cox's Model for Ranking of the Survival Power of the Genes ########
  surv.obj <- Surv(exp(time),censoring)
  coeff.sig <- c(0)
  
  pvalue1.sig <- c(0)
  pvalue2.sig <- c(0)
  
  calcCox = function(x){
    q1 <- unlist(summary(coxph(surv.obj ~ ., data = as.data.frame(x))))
    return(q1$logtest.pvalue)
  }
  
  
  pvalue1.sig <- apply(Y1.pre.train,2,calcCox)
  pvalue2.sig <- apply(Y2.pre.train,2,calcCox)
  
  ###### Adjusting p-values for Multiple Test Correction
  pvalue1.sig.adj <- p.adjust(pvalue1.sig, method = "fdr")
  pvalue2.sig.adj <- p.adjust(pvalue2.sig, method = "fdr")
  
  ##### GENE EXPRESSION #################
  #### As the number of features are quite variable choose first a very loose cut-off 
  signature1.loose <- colnames(Y1.pre.train)[(pvalue1.sig < 0.2)] 
  ### Combined the P-values
  pvalue1.combined <- (pvalue1.sig.adj) 
  names(pvalue1.combined) <- colnames(Y1.pre.train)
  ## Sort it
  pvalue1.combined.sort <- sort(pvalue1.combined)
  ## Only select those genes which are loosely in the signature
  pvalue1.combined.adj <- pvalue1.combined.sort[names(pvalue1.combined.sort) %in% signature1.loose]
  ### Take the top 20 genes ####
  probes1.signature <- names(pvalue1.combined.adj[1:20])
  
  ##### mi-RNA EXPRESSION #################
  #### As the number of features are quite variable choose first a very loose cut-off 
  signature2.loose <- colnames(Y2.pre.train)[(pvalue2.sig < 0.2)] 
  ### Combined the P-values
  pvalue2.combined <- (pvalue2.sig.adj) 
  names(pvalue2.combined) <- colnames(Y2.pre.train)
  ## Sort it
  pvalue2.combined.sort <- sort(pvalue2.combined)
  ## Only select those genes which are loosely in the signature
  pvalue2.combined.adj <- pvalue2.combined.sort[names(pvalue2.combined.sort) %in% signature2.loose]
  ### Take the top 20 genes ####
  probes2.signature <- names(pvalue2.combined.adj[1:20])
  ### Sending the signature
  relev <- list('signature1.sbc'= probes1.signature,'signature2.sbc'= probes2.signature)

}

