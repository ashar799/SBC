posteriorGMMparametrs = function(c,Y,mu,S, alpha,K, epsilon, W, beta, ro, N, D) {
  
  numclust <- table(factor(c, levels = 1:K))
  activeclust <- which(numclust!=0)
  
  for (j in 1:length(activeclust)) {
    clust <- which(c==activeclust[j])
    
    if (length(clust)==1){
      tempmean <- as.vector(Y[clust,1:D])
      temp.number <- 1
      tempmat2 <- ((ro* temp.number)/(ro + temp.number)) * (as.numeric(tempmean - epsilon ) %o% as.numeric(tempmean - epsilon))
      Zi <- beta*W +  tempmat2
      ## Posterior estimate for precision matrix
      res <- try(rWISHART(1, (beta + temp.number), solve(Zi)), silent=TRUE)
      if (class(res) == "try-error"){
        S[activeclust[j],1:D,1:D]= S[activeclust[j],1:D,1:D]
      } else{
        S[activeclust[j],1:D,1:D] <- rWISHART(1, (beta + temp.number), solve(Zi))
      }
      
      S[activeclust[j],1:D,1:D] <- rWISHART(1, (beta + temp.number), solve((Zi)))
      ## Posterior estimate for the mean
      mu[activeclust[j],1:D] <- rMVN(n=1, mean = (ro*epsilon + (temp.number)* tempmean )/(ro + temp.number ), Q = (ro+ temp.number)*S[activeclust[j],1:D,1:D])$x
    }
    
    
    else{
      tempmatrix <- as.matrix(Y[clust,1:D])
      tempmean <- apply(tempmatrix,2,mean)
      ## Some temporary matrices for calcuating the posterior
      ## The number of elements in the class
      temp.number <- as.numeric(numclust[activeclust[j]])
      tempmat1 <-  matrix(apply(as.matrix(apply(tempmatrix,1,function(x) (x -tempmean)%o% (x-tempmean))),1,sum), nrow=D,ncol=D)
      tempmat2 <- ((ro* temp.number)/(ro + temp.number)) * ((tempmean - epsilon ) %o% (tempmean - epsilon))
      Zi <- beta*W + tempmat1 + tempmat2
      ## Posterior estimate for precision matrix
      res <- try(rWISHART(1, (beta + temp.number), solve(Zi)), silent=TRUE)
      if (class(res) == "try-error"){
        S[activeclust[j],1:D,1:D]= S[activeclust[j],1:D,1:D]
      } else{
        S[activeclust[j],1:D,1:D] <- rWISHART(1, (beta + temp.number), solve(Zi))
      }
      # S[activeclust[j],1:D,1:D] <- rWISHART(1, (beta + temp.number), solve(Zi))
      ## Posterior estimate for the mean
      mu[activeclust[j],1:D] <- rMVN(n=1, mean = (ro*epsilon + (temp.number)* tempmean )/(ro + temp.number ), Q = (ro+ temp.number)*S[activeclust[j],1:D,1:D])$x
    }
  } 
  
list('mean' = mu,'precision' = S )
}
