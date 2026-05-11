# BART-Alone Baseline for the Average Causal Effect
# -------------------------------------------------
# Plain BART (no spline extrapolation) used as a baseline against the
# BART+SPL / SBART+SPL methods. Fits BART on the training data, draws
# from the posterior predictive on the counterfactual test design, and
# returns posterior summaries of the ITE plus the ACE posterior via the
# Bayesian bootstrap.
#
# Args:
#   xtr : training design matrix; first column is the binary exposure
#   ytr : observed outcome on the training units
#   xte : test design matrix with the exposure flipped (counterfactuals)
#
# Depends on dbarts and the project-internal helper aceBB().
#
# Returns a named list:
#   iceavg : posterior-mean ITE per unit
#   icelw  : 2.5%  posterior quantile of the per-unit posterior predictive
#   icehi  : 97.5% posterior quantile of the per-unit posterior predictive
#   ace_pd : ACE posterior draws via the Bayesian bootstrap

bartalone<-function(xtr,ytr,xte){
  
  bartps<-dbarts::bart(x.train=xtr,y.train=ytr,x.test=xte)
  ppd_test<-t(apply(bartps$yhat.test,1,function(x) rnorm(n=length(x),mean=x,sd=bartps$sigma)))
  ppd_test_mean<-apply(ppd_test,2,mean)
  
  ## individual causal effects ##
  iceavg<-rep(NA,length(ytr))
  iceavg[which(xtr[,1]==1)]<-ytr[which(xtr[,1]==1)]-ppd_test_mean[which(xtr[,1]==1)]
  iceavg[which(xtr[,1]==0)]<-ppd_test_mean[which(xtr[,1]==0)]-ytr[which(xtr[,1]==0)]
  icelw<-apply(ppd_test,2,quantile,probs=0.025)
  icehi<-apply(ppd_test,2,quantile,probs=0.975)
  
  ## average causal effects ##
  ppd_ice<-matrix(NA,nrow=nrow(ppd_test),ncol=length(ytr))
  for (j in 1:length(ytr)){
    if (xtr[j,1]==1) ppd_ice[,j]<-ytr[j]-ppd_test[,j]
    else ppd_ice[,j]<-ppd_test[,j]-ytr[j]
  }
  ## get ACE posterior using the Bayesian bootstrap ##
  ace_pd<-apply(ppd_ice,1,aceBB)
  
  return(list(iceavg=iceavg,icelw=icelw,icehi=icehi,ace_pd=ace_pd))
}