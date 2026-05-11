# Gutman and Rubin's (2015) Method
# -------------------------------------------------
# Multiple-imputation-style ATE estimator that strata-checks the
# propensity-score distribution, fits separate Bayesian linear regressions
# for the treated and control outcomes (on a natural spline of logit-PS
# plus PS-orthogonalised covariates), and imputes the missing potential
# outcomes from posterior draws. Returns posterior individual and average
# causal effects along with a 0/1 `use` flag.
#
# The flag also lets the function double as a *usability check* on a
# simulated draw: if any PS stratum has fewer than 3 treated or 3 control
# units, the function aborts early and returns just `list(0)` so the
# caller can reject the draw. This is how dgp_bartspl() / dgp_bcf() use it.
#
# Args:
#   Y   : numeric outcome vector, length n
#   trt : 0/1 treatment indicator, length n
#   ps  : estimated propensity scores, length n
#   X   : covariate matrix, n rows
#   M   : number of posterior draws kept (thinned from 10000 MCMC samples)
#   qps : quantile cutpoints used to bin ps into PS strata (e.g.
#         quantile(ps, c(0, .3, .4, .5, .6, .7, 1)))
#
# Depends on MCMCpack::MCMCregress and splines::ns.
#
# Returns:
#   If any PS stratum has < 3 treated or < 3 control units:
#     list(0)              -- positional, signals "reject this draw"
#   Otherwise, a named list:
#     use    : 1
#     iceavg : posterior-mean ITE per unit (length n)
#     icelw  : 2.5%  posterior quantile of the per-unit ITE
#     icehi  : 97.5% posterior quantile of the per-unit ITE
#     ace    : ACE posterior draws (length M, mean across units per draw)


gr<-function(Y,trt,ps,X,M,qps){
	## step 1: create subclasses based on the PS ##
	ps.<-as.numeric(cut(ps,qps,include.lowest = T,right=F))
	
	use<-1
	for (i in 1:length(unique(ps.))){
		if (length(which(ps.==i & trt==0))<3 | length(which(ps.==i & trt==1))<3) use<-0
	}
	
	if (use==0){
		return(list(use))
	}
	else{
		## steps 2 & 3: estimate splines separately for treated and controls & sample M times from the posterior ##
		pstrans<-log(ps/(1-ps))
		psspline<-ns(pstrans,knots = qps[-c(1,length(qps))])
		xort<-NULL
		for (i in 1:ncol(X)){
			xort<-rbind(xort,lm(pstrans~X[,i])$resid)
		}
		xort<-t(xort)
		
		Ytrt<-Y[which(trt==1)]
		pssplinetrt<-psspline[which(trt==1)]
		xorttrt<-xort[which(trt==1)]
		Yctl<-Y[which(trt==0)]
		pssplinectl<-psspline[which(trt==0)]
		xortctl<-xort[which(trt==0)]
		
		keep<-sample(1:10000,M)
		
		outtrt<-as.matrix(MCMCregress(Ytrt~pssplinetrt+xorttrt,burnin = 1000,mcmc=10000,verbose=F))[keep,]
		outctl<-as.matrix(MCMCregress(Yctl~pssplinectl+xortctl,burnin = 1000,mcmc=10000,verbose=F))[keep,]
		
		## step 4: impute the missing potential outcomes for each sample from the posterior##
		imptrt<-list()
		impctl<-list()
		for (i in 1:M){
			imptrt<-c(imptrt,list(rowSums(matrix(outctl[i,(-ncol(outctl))],nrow=length(Ytrt),ncol=ncol(outctl)-1,byrow=T)*cbind(1,pssplinetrt,xorttrt))))
			impctl<-c(impctl,list(rowSums(matrix(outtrt[i,(-ncol(outtrt))],nrow=length(Yctl),ncol=ncol(outtrt)-1,byrow=T)*cbind(1,pssplinectl,xortctl))))
		}
		imptrt<-as.matrix(as.data.frame(imptrt))
		impctl<-as.matrix(as.data.frame(impctl))
		
		## step 6: estimate the treatment effect and its sample variance for each dataset ##
		indtrteff<-t(rbind(matrix(Ytrt,nrow=length(Ytrt),ncol=M,byrow=F)-imptrt,
											 impctl-matrix(Yctl,nrow=length(Yctl),ncol=M,byrow=F)))
		
		iceavg<-apply(indtrteff,2,mean)
		icelw<-apply(indtrteff,2,quantile,probs=0.025)
		icehi<-apply(indtrteff,2,quantile,probs=0.975)
		ace<-rowMeans(indtrteff)
		
		
		
		return(list(use=use,iceavg=iceavg,icelw=icelw,icehi=icehi,ace=ace))
	}
}