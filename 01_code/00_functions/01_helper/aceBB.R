# Bayesian Bootstrap for the Average Causal Effect
# ------------------------------------------------
# Draws one posterior sample of the ACE using Rubin's (1981) Bayesian
# bootstrap, following the supplementary material of Wang et al. (2015).
# Weights come from a Dirichlet posterior (Rubin's improper
# prior, multiplicity-one data) applied to the unit-level contrasts in x.
#
# Args:
#   x     : numeric vector of individual treatment effects (e.g. Y(1) - Y(0))
#   nboot : number of bootstrap replicates (default 250)
#
# Returns:
#   A single numeric draw from the ACE posterior.
#
# References:
#   Rubin (1981), Annals of Statistics 9(1), 130-134.
#   Wang et al. (2015), supplementary material.
#   Nethery et al. (2019), Annals of Applied Statistics

aceBB<-function(x,nboot=250){

	# Replicate unit-level effects across bootstrap rows
  diffpotmat<-matrix(rep(x,nboot),nrow=nboot,byrow=T)
  
  # Dirichlet(1, ..., 1) weights via the Exp(1) / sum trick
  dirichlet_sample <- matrix( rexp(length(x) * nboot, 1) , nrow=nboot, byrow = TRUE)
  dirichlet_sample <- dirichlet_sample / rowSums(dirichlet_sample) # normalize to sum to 1
  
  # Posterior of the ACE: weighted average per replicate
  posteriorace<-rowSums(dirichlet_sample*diffpotmat)
  
  # Return one draw
  MCMCace<-sample(posteriorace,size=1)
  
  return(MCMCace)
}