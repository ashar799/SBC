### This file loads SBC gene signature ####
loadSBC-BC-signature = function(Y.ashar.pretrain, time, censoring){
  
  
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
  
  for ( i in 1:ncol(Y.ashar.pretrain)){
    q1 <- unlist(summary(coxph(surv.obj ~ Y.ashar.pretrain[,i], data = as.data.frame(Y.ashar.pretrain))))
    pvalue.sig[i] <- q1$logtest.pvalue 
    q2 <- unlist(summary(aov(as.numeric(censoring) ~ as.numeric(Y.ashar.pretrain[,i]), data=as.data.frame(Y.ashar.pretrain)) ))
    pvalue.anova[i] <- as.numeric(q2[9]) 
    
  }
  
  ###### Adjusting p-values for Multiple Test Correction
  pvalue.sig.adj <- p.adjust(pvalue.sig, method = "fdr")
  pvalue.anova.adj <-  p.adjust(pvalue.anova, method = "fdr")
  
  
  
  sum(((pvalue.anova.adj < 0.04) & (pvalue.sig.adj < 0.05)) + 0)
  
  signature <-  ((pvalue.anova.adj < 0.04) & (pvalue.sig.adj < 0.05))
  
  probes.signature <- colnames(Y.ashar.pretrain)[signature]
  
  
  
  
  relev <- list('signature.sbc'= signature)
  
}