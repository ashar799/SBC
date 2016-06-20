posterioralpha = function(c, N, alpha, shape.alpha, rate.alpha) {
  
  ## Calculate the number of Active clusters
  ## Calculate the number of Active clusters
  
  g <- table(factor(c, levels = 1:K))
  active <- which(g!=0)
  Kplus <- length(active)
  
  
#   f = function(x, N = N, Kp = Kplus){
#     (Kp - 1.5)* log (x) + lgamma(x) - lgamma(N+x) - (0.5/(x))
#     
#   }
#   
#   fprima = function(x, N = N, Kp = Kplus){
#     (Kp - 1.5)* (1/x) + digamma(x) - digamma(N+x) + (0.5/ (x)^2)
#     
#   }
#   
#   ars(1, f, fprima, x = alpha , m =1, N = N, Kp = Kplus ) 
#   

### Introduce Eta
 eta <- rbeta(1,alpha+1,N)
  p1 <- (shape.alpha + Kplus -1 )/(shape.alpha + Kplus -1+ N*(rate.alpha - log(eta)))
 p2 <- 1- p1
indt <- sample(c(1,2), 1, prob = c(p1,p2), replace = TRUE)
if (indt == 2){
  alpha <- rgamma(1, shape = shape.alpha + Kplus -1, rate = rate.alpha - log(eta))
} else{
  alpha <-  rgamma(1, shape = shape.alpha + Kplus, rate = rate.alpha - log(eta))
}

alpha
}

