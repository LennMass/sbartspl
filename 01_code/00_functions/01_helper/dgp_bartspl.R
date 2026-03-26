
dgp_bartspl <- function(N = 500, #as.numeric(args[4]) #total N#
												ncovs = 25, #-as.numeric(args[5]) #number of additional "potential confounders"#
												a_0 = .1,  # range of the overlap set
												b_0 = 7 # points in the overlap set
												){
	
	n1<-N/2
	n0<-N/2
	
	
	max_iter <- 100
	iter <- 0
	
	repeat {
		iter <- iter + 1
		
		
		x<-c(rep(1,n1),rep(0,n0))
		
		
		## construct 10 true confounders ##
		covs<-matrix(c(rbinom(n1,size=1,prob=.45),rbinom(n0,size=1,prob=.4),rnorm(n=n1,mean=2,sd=2),rnorm(n=n0,mean=1.3)),nrow=N,byrow = F)
		for (i in 1:4){
			covs<-cbind(covs,c(rbinom(n1,size=1,prob=.45),rbinom(n0,size=1,prob=.4)),c(rnorm(n=n1,mean=2,sd=2),rnorm(n=n0,mean=1.3)))
		}
		
		if (ncovs>0){
			## make additional ncov covariates ##
			covs<-cbind(covs,matrix(rnorm(ncovs*N),nrow=N,ncol=ncovs))
		}
		
		## estimate misspecified PS ##
		emod<-stats::glm(x~covs,family=binomial())
		ps_mis<-emod$fitted
		
		## find RO and RN ##
		RO_mis<-pw_overlap(ps=ps_mis,E=x,a=a_0,b=b_0)
		
		RO_mean <- mean(RO_mis)
		
		## true potential outcomes and causal effects ##
		Y0<-rowSums(.5*covs[,c(1,3,5,7,9)])+(15/(1+exp(-(covs[,2]*8-1))))+rowSums(covs[,c(4,6,8,10)])-5
		Y1<-rowSums(covs[,c(1,3,5,7,9)]-.5*covs[,c(2,4,6,8,10)])
		ce_true<-Y1-Y0
		
		## observed outcome data ##
		Yobs<-rep(NA,N)
		Yobs[which(x==0)]<-Y0[which(x==0)]
		Yobs[which(x==1)]<-Y1[which(x==1)]
		
		## create a dataset with PS and all confounders ##
		datall_mis<-data.frame(Yobs,x,ps_mis,covs)
		#datpso_mis<-data.frame(Yobs,x,ps_mis)
		X_untrim <- data.frame(x, ps_mis, covs)
		
		## create a dataset with PS and all confounders, but trimmed ##
		datall_mis_trim <- datall_mis[which(RO_mis==1), ]
		X_trim <- data.frame(x, ps_mis, covs)[which(RO_mis==1), ]
		
		
		# check GR
		ce_gr_mis<-gr(Y=datall_mis$Yobs,trt=datall_mis$x,ps=datall_mis$ps_mis,X=covs,M=500,qps=quantile(datall_mis$ps_mis,probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))
		keepsim_mis<-ce_gr_mis[[1]]
		
		
		
		# safeguards 
		if( (sum(1-RO_mis)>0) && keepsim_mis==1) {
			
			break
			
		}
		
		if (iter >= max_iter) {
			stop("Timeout: condition not met after ", max_iter, " iterations")
		}
		
	}
	
	
	return(list(untrimmed_dat = datall_mis,
							trimmed_dat = datall_mis_trim,
							covs=covs,
							X_trim=X_trim, 
							X_untrim=X_untrim,
							ate_true=mean(ce_true),
							ite_true=ce_true,
							RO=RO_mis,
							RO_share=RO_mean, 
							p_noise = ncovs 
	))
	
}

