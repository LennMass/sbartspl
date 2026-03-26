aceBB<-function(x,nboot=250){
  # see supplementary material of Wang et al. (2015)
  
  #For βαY, we assume by default a uniform prior, though in specific applications one can also 
  # specify informative priors.
  diffpotmat<-matrix(rep(x,nboot),nrow=nboot,byrow=T)
  
  # Prior sampling
  # For the prior of θ, we assume π(θ) ∝ ΠKk=1θ−1k (Rubin, 1981). 
  # The posterior follows a Dirichlet distribution D(n1, . . . , nK )
  # p(θ|V ) ∝ ΠKk=1θnk −1k 
  # 
  dirichlet_sample <- matrix( rexp(length(x) * nboot, 1) , nrow=nboot, byrow = TRUE)
  dirichlet_sample <- dirichlet_sample / rowSums(dirichlet_sample) # normalize to sum to 1
  
  # Calculate posterior for each bootstrap sample
  # p(βαY, θ|D) = p(βαY|D)p(θ|D)
  posteriorace<-rowSums(dirichlet_sample*diffpotmat)
  
  ## sample from the posterior ##
  MCMCace<-sample(posteriorace,size=1)
  
  return(MCMCace)
}