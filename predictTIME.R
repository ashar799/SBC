### This function takes the posterior parameters AND predicts the time for the new points
#### The fundamental assumption is that EACH NEW TEST POINT IS CONDITIONALLY INDEPENDENT on the OTHER POINTS
#### We predict value of one point at a time
### The final output is Time for the new samples, ONE AT A TIME
#### when making a final comparison it's a good idea to weigh the points according to the model likelihood

predictchineseAFTtime = function(Y.input){
  
  source('priorPARAMETERS.R')
  
  Y.new  <- Y.input
  N.new <- nrow(Y.new)
  c.new.list <- list(0)
  ## The number of posterior samples
  #That.new <- time.new 
  post.time  = matrix(NA,nrow = nrow(Y.new), ncol = Nps)
  print("GOING THROUGH MCMC Samples")
  pb <- txtProgressBar(min = 1, max = Nps , style = 3)
  
  cind <- c(0)
  
  # modelweights <- c(0)
  
  predCIndex.sbc <- c(0)
  
  
  for (count in 1:Nps){
    
    ## Assign the parameters to the posterior sample
    ctemp <- c.list[[count]]
    mu <- mu.list[[count]]
    S <- S.list[[count]]
    beta0 <- beta0.list[[count]]
    betahat  <- betahat.list[[count]] 
    sigma2  <- sigma2.list[[count]] 
    g <- table(factor(ctemp, levels = 1:K))
    activeclass <- which(g!=0)
    ## The table function helps converting the data point specific indicator variables to class specific indicator variables
    kminus <- length(activeclass)
  # active <- activeclass
    #Two Auxilary Variables
    #The name of the auxilary variables are taken to be one and two more than the maximum value in the already active cluster set
    activeclass <- append(activeclass, max(activeclass)+1)
    activeclass <- append(activeclass, max(activeclass)+1)
    active <- activeclass 
    


    ### Assigning values to parameters 
  
    priortwo <- NA
    priortwo <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
    mu[active[kminus+1],1:D]  <- priortwo$mu  
    S[active[kminus+1],1:D,1:D]  <- priortwo$Sigma[1:D,1:D]  
    beta0[active[kminus+1]] <- priortwo$beta0 
    sigma2[active[kminus+1]] <- priortwo$sigma2
    betahat[active[kminus+1],1:D] <- priortwo$betahat 
   
    
    
   
    priorthree <- NA
    priorthree <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
    mu[active[kminus+2],1:D]  <- priorthree$mu  
    S[active[kminus+2],1:D,1:D]  <- priorthree$Sigma[1:D,1:D]  
    beta0[active[kminus+2]] <- priorthree$beta0 
    sigma2[active[kminus+2]] <- priorthree$sigma2
    betahat[active[kminus+2],1:D] <- priorthree$betahat 
   
    
    #######################################################
    ########### Scaling the data ########################
    Y.new.scaled.list <- list(0)
    for (j in 1:kminus) {
     clust <- which(ctemp == active[j])
     
     if(length(clust) > 1){
     obj.t <- scale(Y[clust,1:D], center = TRUE, scale = TRUE)
     Y.new.scaled.list[[j]] <- scale(Y.new, center = attr(obj.t,"scaled:center"), scale = (attr(obj.t,"scaled:scale")))
     } else {
       obj.t <- scale(Y[,1:D], center = TRUE, scale = TRUE)
      Y.new.scaled.list[[j]] <- scale(Y.new, center = attr(obj.t,"scaled:center"), scale = (attr(obj.t,"scaled:scale")))
       
     }
  
    
    }
    
    for (j in (kminus+1):(kminus+2)) {
     obj.t <- scale(Y[,1:D], center = TRUE, scale = TRUE)
  
     Y.new.scaled.list[[j]] <- scale(Y.new, center = attr(obj.t,"scaled:center"), scale = (attr(obj.t,"scaled:scale")))
     }

    
  ###### Some quantities used to store probabilities  
  posteriortime <- matrix(0, nrow = length(active), ncol = N.new)
  posteriortimeweight <- matrix(0, nrow = length(active), ncol = N.new)
  weights <- matrix(0, nrow = length(active), ncol = N.new)
 
  

    ## This can't be parallelized !!!!!
    for(l in 1:N.new)  {
      
      ## Calculating the Expectations and also the normalization constant for the Expectation
      for (j in 1:kminus) {
        
        posteriortime[j,l] <-  beta0[active[j]] + betahat[active[j],1:D] %*% as.vector(t(Y.new.scaled.list[[j]][l,])) 
        
        posteriortimeweight[j,l] <- log(g[active[j]]/ (N-1+alpha))  +  dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[j],1:D], Q = S[active[j],1:D,1:D], log =TRUE) 
      }
      
      res <- try(dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]), silent=TRUE)
      if (class(res) == "try-error"){
        posteriortime[kminus+1,l] <- 0
        posteriortimeweight[kminus+1,l] <- -Inf
      } else{
        posteriortime[kminus+1,l] <-  beta0[active[kminus+1]] + betahat[active[kminus+1],1:D] %*% as.vector(t(Y.new.scaled.list[[j]][l,]))
        posteriortimeweight[kminus+1,l] <-  log(alpha/ (N-1+alpha)) +  dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D], log = TRUE) 
      }
      
      res2 <- try(dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+2],1:D], Q= S[active[kminus+2],1:D,1:D]), silent=TRUE)
      if (class(res) == "try-error"){
        posteriortime[kminus+2,l] <- 0
        posteriortimeweight[kminus+2,l] <- -Inf
        
      } else{
        posteriortime[kminus+2,l] <- beta0[active[kminus+2]] + betahat[active[kminus+2],1:D] %*% as.vector(t(Y.new.scaled.list[[j]][l,])) 
        posteriortimeweight[kminus+2,l] <-  log(alpha/ (N-1+alpha))  +  dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+2],1:D], Q= S[active[kminus+2],1:D,1:D], log = TRUE)      
      }
      
    

        weights[,l] <- exp(posteriortimeweight[,l])/sum(exp(posteriortimeweight[,l]))
    }

    
        
   for ( l in 1:N.new){
     post.time[l,count] <- as.numeric(t(posteriortime[,l]) %*% weights[,l])
     
   } 
    
 # modelweights[count] <- sum(exp((1/N.new) *apply(posteriortimeweight,1,sum)))
    
  predCIndex.sbc[count] <- as.numeric(survConcordance(smod.new ~ exp(-post.time[,count]))[1]) 
# 
#    print(cind[count])

## Calculting the Model Weight

    Sys.sleep(0.1)
    setTxtProgressBar(pb, count)
    
    
  }
  
#   #### To calculate average values over MCMC samples
#   modelweight.norm <- modelweights/(sum(modelweights))
# 
#  post.time.corrected <- post.time 
#  for ( i in 1:Nps){
#  post.time.corrected[,i] <- post.time[,i] *modelweight.norm[i] 
#  }
# post.time.avg <<- apply(post.time.corrected[,1:Nps],1,sum)
  predCIndex.sbc <<- predCIndex.sbc
  post.time <<- post.time
}
