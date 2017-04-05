
### Based on the DPpackage for Bayesian Density Regression
### This function uses the DP package to estimate some parameters
CompareDP = function(){
  
  
  ### Example from the package ###
  
  dtrue <- function(grid, x) {
    exp(-2 * x) * dnorm(grid, mean = x, sd = sqrt(0.01)) +(1 - exp(-2 * x)) * dnorm(grid, mean = x^4, sd = sqrt(0.04))
     }
  
  mtrue <- function(x) exp(-2 * x) * x + (1 - exp(-2 * x)) * x^4
  
  prior <- list(a0 = 10, b0 = 1, nu1 = 4, nu2 = 4, s2 = 0.5 * wcov,
                m2 = wbar, psiinv2 = 2 * solve(wcov), tau1 = 6.01, tau2 = 3.01)
  mcmc <- list(nburn = 5000, nsave = 5000, nskip = 3, ndisplay = 1000)
  
  nrec <- 500
  x <- runif(nrec)
  y1 <- x + rnorm(nrec, 0, sqrt(0.01))
  y2 <- x^4 + rnorm(nrec, 0, sqrt(0.04))
  u <- runif(nrec)
  prob <- exp(-2 * x)
  y <- ifelse(u < prob, y1, y2)
  w <- cbind(y, x)
  wbar <- apply(w, 2, mean)
  wcov <- var(w)
  
  fitWDDP <- DPcdensity(y = y, x = x, xpred = seq(0, 1, 0.02),
              ngrid = 100, compute.band = TRUE, type.band = "HPD",
               prior = prior, mcmc = mcmc, state = NULL, status = TRUE)
  
  
  ### My own example did not work !!! :(())
  w <- cbind(time, Y)
  wbar <- apply(w, 2, mean)
  wcov <- var(w)
  
  prior <- list(a0 = 10, b0 = 1, nu1 = 60, nu2 = 60, s2 = 0.5 * wcov,
                m2 = wbar, psiinv2 = 2*solve(wcov), tau1 = 6.01, tau2 = 3.01)
  
  mcmc <- list(nburn = 5000, nsave = 5000, nskip = 3, ndisplay = 1000)
  
  fitWDDP <- DPcdensity(y = time, x = Y, xpred = Y.new, 
                        prior = prior, mcmc = mcmc, state = NULL, status = TRUE,
                        compute.band = TRUE, type.band = "PD")
  
  
  
}
