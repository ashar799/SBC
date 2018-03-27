### This file loads SBC gene signature ####
loadSBCvijversignature = function(Y.pre.train, time, censoring){
  
  
  # train.index <- pheno.vijver$Label_Traing_and_Validation == 'Training'
  # test.index <- pheno.vijver$Label_Traing_and_Validation == 'Validation'
   ######## Prefiltering of the Genes ############################### ###########################
  ######### Using Univariate Cox's Model for Ranking of the Survival Power of the Genes ########
  ######## Using T-test (or K way Anova) for Ranking of Clustering Power of the Genes ###########
  ######## The Class Labels ARE WHETHER METASTATIS OCCURED OR NOT ###############################
  
  surv.obj <- Surv(time,censoring)
  coeff.sig <- c(0)
  
  pvalue.sig <- c(0)
  pvalue.anova <- c(0)
  
  calcCox = function(x){
    q1 <- unlist(summary(coxph(surv.obj ~ ., data = as.data.frame(x))))
    return(q1$logtest.pvalue)
    }
  calcANOVA = function(x){
    q2 <- unlist(summary(aov(as.numeric(censoring) ~ ., data=as.data.frame(x)) ))
   return(as.numeric(q2[9])) 
  }
  
  pvalue.sig <- apply(Y.pre.train,2,calcCox)
  pvalue.anova <- apply(Y.pre.train,2,calcANOVA)
  
  ###### Adjusting p-values for Multiple Test Correction
  pvalue.sig.adj <- p.adjust(pvalue.sig, method = "fdr")
  pvalue.anova.adj <-  p.adjust(pvalue.anova, method = "fdr")
  
  
  #### As the number of features are quite variable choose first a very loose cut-off 
  
  signature.loose <- colnames(Y.pre.train)[(pvalue.anova.adj < 0.05) & (pvalue.sig.adj < 0.05)] 
  
  ### Combined the P-values
  pvalue.combined <- (pvalue.sig.adj + pvalue.anova.adj) 
  names(pvalue.combined) <- colnames(Y.pre.train)
  ## Sort it
  pvalue.combined.sort <- sort(pvalue.combined)
  ## Only select those genes which are loosely in the signature
  pvalue.combined.adj <- pvalue.combined.sort[names(pvalue.combined.sort) %in% signature.loose]
  ### Take the top 50 genes ####
  probes.signature <- names(pvalue.combined.adj[1:50])
  
  
  
  relev <- list('signature.sbc'= probes.signature)
  
}