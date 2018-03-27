### This function takes the posterior parameters AND predicts CLUSTER MEMEBERSHIP for the new points
#### The fundamental assumption is that EACH NEW TEST POINT IS CONDITIONALLY INDEPENDENT on the OTHER POINTS
#### We predict value of one point GIVEN ONLY ITS MOLECULAR DATA
### The final output is Time for the new samples, ONE AT A TIME


predictCLASS = function(Y.input){
  
  
  
  source('priorPARAMETERS.R')
  Y.new  <- Y.input
  N.new <- nrow(Y.new)
  c.new.list <- list(0)
  ## The number of posterior samples

  print("GOING THROUGH MCMC Samples")
  pb <- txtProgressBar(min = 1, max = Nps , style = 3)
  
  

  ctemp.new <- c(0)
  modelweights <- c(0)
  
  
  for (count in 1:Nps){
    
    ## Assign the parameters to the posterior sample
    ctemp <- c.list[[count]]
    mu <- mu.list[[count]]
    S <- S.list[[count]]
    
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
   
    
    
    
    
    priorthree <- NA
    priorthree <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
    mu[active[kminus+2],1:D]  <- priorthree$mu  
    S[active[kminus+2],1:D,1:D]  <- priorthree$Sigma[1:D,1:D]  
   
    
    
    ###### Some quantities used to store probabilities  
    posteriorweight <- matrix(0, nrow = length(active), ncol = N.new)
    weights <- matrix(0, nrow = length(active), ncol = N.new)
    weights.final <- c(0)
    ctemp.new <- c(0)
    ## This can't be parallelized !!!!!
    for(l in 1:N.new)  {
      
      ## Calculating the Expectations and also the normalization constant for the Expectation
      for (j in 1:kminus) {
        posteriorweight[j,l] <- log(g[active[j]]/ (N-1+alpha))  +  dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[j],1:D], Q = S[active[j],1:D,1:D], log =TRUE) 
      }
      
      res <- try(dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]), silent=TRUE)
      if (class(res) == "try-error"){
        posteriorweight[kminus+1,l] <- -Inf
      } else{
        posteriorweight[kminus+1,l] <-  log(alpha/ (N-1+alpha)) +  dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D], log = TRUE) 
      }
      
      res2 <- try(dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+2],1:D], Q= S[active[kminus+2],1:D,1:D]), silent=TRUE)
      if (class(res) == "try-error"){
        posteriorweight[kminus+2,l] <- -Inf
        
      } else{
        posteriorweight[kminus+2,l] <-  log(alpha/ (N-1+alpha))  +  dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[kminus+2],1:D], Q= S[active[kminus+2],1:D,1:D], log = TRUE)      
      }
      
       weights[,l] <- exp(posteriorweight[,l])/sum(exp(posteriorweight[,l]))
    
      
      if (sum(exp(posteriorweight[,l])) < 1e-200){
        ctemp.new[l] <- sample(active, 1, prob = rep(1,length(active)), replace = TRUE)
      } else {  
        ctemp.new[l] <- sample(active, 1, prob= weights[,l], replace = TRUE)
        
      }
      
     weights.final[l] <- dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[ctemp.new[l],1:D], Q = S[ctemp.new[l],1:D,1:D], log =TRUE) 
       
     }
    
    
    modelweights[count] <- sum(weights.final)
      
    c.new.list[[count]] <- ctemp.new
    Sys.sleep(0.1)
    setTxtProgressBar(pb, count)
    
    }
  
  
  
  
  ## Converting the list to a matrix
  
  c.matrix.new <- matrix(NA, nrow = N.new, ncol = Nps)
  for( h in 1:Nps){
  c.matrix.new[,h] <- c.new.list[[h]]
  }
  
  
  c.matrix.new <<- c.matrix.new
  
  
  
  
  
  
  # ## As the clusters are different we switch the labels
  # for ( i in 1:N.new){
  #   if(mpear2$cl[i] ==1){
  #     c.sbc.new[i] = 1
  #   } else {
  #     c.sbc.new[i] = 2
  #   }}
  # 
  # 
  test.modelweights <<- modelweights
  
 }
