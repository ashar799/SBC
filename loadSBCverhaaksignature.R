### This file loads SBC gene signature ####
loadSBCverhaaksignature = function(Y.pre.train, time, censoring){
  
  
  
  ######## Prefiltering of the Genes ############################### ###########################
  ######### Using Univariate Cox's Model for Ranking of the Survival Power of the Genes ########
  surv.obj <- Surv(time,censoring)
  coeff.sig <- c(0)
  
  pvalue.sig <- c(0)
  
  
  calcCox = function(x){
    q1 <- unlist(summary(coxph(surv.obj ~ ., data = as.data.frame(x))))
    return(q1$logtest.pvalue)
  }
  
  
  pvalue.sig <- apply(Y.pre.train,2,calcCox)
  
  
  ###### Adjusting p-values for Multiple Test Correction
  pvalue.sig.adj <- p.adjust(pvalue.sig, method = "fdr")
  
  #### As the number of features are quite variable choose first a very loose cut-off 
  
  signature.loose <- colnames(Y.pre.train)[(pvalue.sig.adj < 0.5)] 
  
  ### Combined the P-values
  pvalue.combined <- (pvalue.sig.adj) 
  names(pvalue.combined) <- colnames(Y.pre.train)
  ## Sort it
  pvalue.combined.sort <- sort(pvalue.combined)
  ## Only select those genes which are loosely in the signature
  pvalue.combined.adj <- pvalue.combined.sort[names(pvalue.combined.sort) %in% signature.loose]
  ### Take the top 50 genes ####
  probes.signature <- names(pvalue.combined.adj[1:50])
  
  relev <- list('signature.sbc'= probes.signature)
  
}