# SBART+SPL function for application results

sbartspl<-function(datall,RO,nburn=10000,nsamp=5000, 
									 var_infl=0, # no variance inflation within non-overlap region
									 n_forest=4, # number of separate/parallel chains
									 probs.rcs.ps=c(.1,.33,.66,.9), # knots for restriced cubic splines for ps
									 probs.rcs.Y=c(.2,.4,.6,.8) # knots for restriced cubic splines for Y
){
	
	## this function implements SBART+SPL method (for continuous outcomes) based on BART+SPL (Nethery et al. 2019)
	## first column in datall should be the observed outcome variable
	## second column should be the binary exposure indicator
	## third variable should be the PS or confounder on which to base the non-overlap
	## any other variables to be included in the model are in the fourth column and beyond
	
	names(datall)[1:3]<-c('Yobs','x','ps')
	if (ncol(datall)>3) names(datall)[4:ncol(datall)]<-paste('u',1:(ncol(datall)-3),sep='')
	
	#### adapt nburn/nsamp according to number of parallel chains
	nsamp <- nsamp/n_forest
	nburn <- nburn/n_forest
	
	#### Prepare data for SoftBART sampler
	# preprocess data for softbart_sampler
	prep_X <- preprocess_df(datall[, !(colnames(datall) == "Yobs")])
	quant_X <- quantile_normalize_bart(prep_X$X)
	group_X <- prep_X$group
	scale_Y <- normalize_bart(datall$Yobs)
	# create overlap dataset
	datov_scale <- data.frame(Yobs=scale_Y, quant_X)[which(RO==1),]
	## create a test overlap dataset used to predict counterfactuals with BART ##
	datovtest_scale<-datov_scale
	datovtest_scale$x<-1-datovtest_scale$x
	
	## create an overlap dataset ##
	datov<-datall[which(RO==1),]
	## create a test overlap dataset used to predict counterfactuals with BART ##
	datovtest<-datov
	datovtest$x<-1-datovtest$x
	## create an overlap dataset for plugging into the spline model ##
	sp_ind<-which(datov$ps>quantile(datov$ps,.05) & datov$ps<quantile(datov$ps,.95))
	datov_sp<-datov[sp_ind,]
	## create a non-overlap dataset ##
	datno<-datall[which(RO==0),]
	## dataset containing units in the RN with E=0 (untreated) ##
	datno0<-datall[which(RO==0 & datall$x==0),]
	## dataset containing units in the RN with E=1 (treated) ##
	datno1<-datall[which(RO==0 & datall$x==1),]
	
	## for everyone in the RN, find the distance from their PS to the nearest PS in the RO ##
	ROdist<-NULL
	for (i in 1:nrow(datno)){
		psi<-datno$ps[i]
		ROdist<-c(ROdist,min(abs(datov$ps-psi)))
	}
	
	#### Initialize MakeForest function for SBART model
	# Default specs for sbart model
	sb_forest <- function(tau_rate_sb=10,
						  alpha_sb=1, # for smooth bart: ncol(datov_scale[, !(colnames(datov_scale) == "Yobs")])
						  update_tau_sb = TRUE, # for dart: FALSE
						  update_tau_mean_sb = TRUE, # for dart: FALSE
						  update_s_sb = TRUE, # for smooth bart: FALSE
						  update_alpha_sb = TRUE # for smooth bart: FALSE
	){
		SoftBart::MakeForest(hypers=SoftBart::Hypers(Y=datov_scale$Yobs,
																								 X=datov_scale[, !(colnames(datov_scale) == "Yobs")], # design matrix X normalized between 0 and 1
																								 group=group_X,
																								 tau_rate=tau_rate_sb,
																								 alpha=alpha_sb
		), 
		opts=SoftBart::Opts(num_burn=0, # because we do burn-in and saving later on our own with nburn, nsamp
												num_save=1, # because we do burn-in and saving later on our own with nburn, nsamp
												update_tau=update_tau_sb,
												update_tau_mean=update_tau_mean_sb,
												update_s = update_s_sb, 
												update_alpha = update_alpha_sb 
		)) 
	}
	## initialize MakeForest sampler to implement SBART in the range of overlap
	## initialize separate chains (4 in total)
	sb_chain1 <- sb_forest()
	sb_chain2 <- sb_forest()
	sb_chain3 <- sb_forest()
	sb_chain4 <- sb_forest()
	
	
	## hyperparameters for the spline ##
	# Restricted Cubic Spline Design Matrix
	# rcspline.eval returns a matrix with predictor variable (if inclx=TRUE) followed by nk−2 nonlinear terms. 
	rcs_ps_p <- 1 + (length(probs.rcs.ps) - 2) 
	rcs_Y_p <- 1 + (length(probs.rcs.Y) - 2)
	# dims of X_spl
	p_spl1<-1 + rcs_ps_p + rcs_Y_p +(ncol(datall)-3)
	p_spl0<-1 + rcs_ps_p + rcs_Y_p +(ncol(datall)-3)
	mu0<-matrix(0,nrow=p_spl1,ncol=1)
	Sigma0<-10000*diag(p_spl1)
	Sigma0_inv<-solve(Sigma0)
	a0<-1
	b0<-1
	
	## initialize parameters etc ##
	delta<-matrix(0,nrow=nrow(datov),ncol=1)
	Y1s<-datov_sp$Yobs
	Y0s<-datov_sp$Yobs
	beta1<-matrix(0,nrow=p_spl1,ncol=1)
	sigma_spl1<-1
	beta0<-matrix(0,nrow=p_spl0,ncol=1)
	sigma_spl0<-1
	
	## matrices to store posterior samples ##
	delta_star_save<-matrix(NA,nrow=nsamp,ncol=nrow(datall))
	
	# storing matrices for the sampler
	s_temp_test <- matrix(NaN, ncol=n_forest, nrow=nrow(datovtest)) 
	s_temp_sigma <- matrix(NaN, ncol=n_forest, nrow=1)
	s_temp_pred <- matrix(NaN, ncol=n_forest, nrow=nrow(datovtest))
	
	## run the sampler ##
	for (i in 1:(nburn+nsamp)){
		
		##############################################################################
		## 1. run the SoftBART sampler once and take a sample from the SoftBART ppd ##
		##############################################################################
		
		
		# for softbart_sampler
		# chain 1
		s_temp_test[,1]<- sb_chain1$do_gibbs(X=data.matrix(datov_scale[, !(colnames(datov_scale) == "Yobs")]), 
																				 Y=data.matrix(datov_scale$Yobs),
																				 #X_test=data.matrix(datov[, !(colnames(datov) == "Yobs")]),
																				 X_test=data.matrix(datovtest_scale[,!(colnames(datovtest_scale) == "Yobs")]), 
																				 i=1) %>%
			t()
		s_temp_sigma[,1] <- sb_chain1$get_sigma() # get_sigma() gives error standard deviation of the forest
		# chain 2
		s_temp_test[,2]<- sb_chain2$do_gibbs(X=data.matrix(datov_scale[, !(colnames(datov_scale) == "Yobs")]), 
																				 Y=data.matrix(datov_scale$Yobs),
																				 #X_test=data.matrix(datov[, !(colnames(datov) == "Yobs")]),
																				 X_test=data.matrix(datovtest_scale[,!(colnames(datovtest_scale) == "Yobs")]), 
																				 i=1) %>%
			t()
		s_temp_sigma[,2] <- sb_chain2$get_sigma() # get_sigma() gives error standard deviation of the forest
		# chain 3
		s_temp_test[,3]<- sb_chain3$do_gibbs(X=data.matrix(datov_scale[, !(colnames(datov_scale) == "Yobs")]), 
																				 Y=data.matrix(datov_scale$Yobs),
																				 #X_test=data.matrix(datov[, !(colnames(datov) == "Yobs")]),
																				 X_test=data.matrix(datovtest_scale[,!(colnames(datovtest_scale) == "Yobs")]), 
																				 i=1) %>%
			t()
		s_temp_sigma[,3] <- sb_chain3$get_sigma() # get_sigma() gives error standard deviation of the forest
		# chain 4
		s_temp_test[,4]<- sb_chain4$do_gibbs(X=data.matrix(datov_scale[, !(colnames(datov_scale) == "Yobs")]), 
																				 Y=data.matrix(datov_scale$Yobs),
																				 #X_test=data.matrix(datov[, !(colnames(datov) == "Yobs")]),
																				 X_test=data.matrix(datovtest_scale[,!(colnames(datovtest_scale) == "Yobs")]), 
																				 i=1) %>%
			t()
		s_temp_sigma[,4] <- sb_chain4$get_sigma() # get_sigma() gives error standard deviation of the forest
		
		# construct unnormalized ppd
		ppd_test<-rnorm(n=length(s_temp_test),
										mean=s_temp_test,
										sd=s_temp_sigma) 
		# re-transform 
		ppd_test <- unnormalize_bart(ppd_test, min(datall$Yobs), max(datall$Yobs))
		
		## form the predicted individual causal effects in the RO from the BART predictions ##
		delta[which(datov$x==1),]<-datov$Yobs[which(datov$x==1)]-ppd_test[which(datov$x==1)]
		delta[which(datov$x==0),]<-ppd_test[which(datov$x==0)]-datov$Yobs[which(datov$x==0)]
		delta_sp<-delta[sp_ind,]
		
		## update Y1s ##
		foosp<-s_temp_test[sp_ind]
		Y1s[which(datov_sp$x==0)]<-foosp[which(datov_sp$x==0)]
		Y0s[which(datov_sp$x==1)]<-foosp[which(datov_sp$x==1)]
		
		## fit spline to the estimated causal effects in the RO from BART ##
		delta_star<-rep(NA,sum(RO==0))
		bs_ps<-rcspline.eval(datov_sp$ps,knots=quantile(datov_sp$ps,probs=probs.rcs.ps),inclx = T)
		
		## spline for treated units (E=1) in the RN ##
		if (nrow(datno1)>0){
			
			####################################
			## 2. run the spline sampler once ##
			####################################
			
			bs_y1s<-rcspline.eval(Y1s,knots=quantile(Y1s,probs=probs.rcs.Y),inclx = T)
			if (ncol(datov_sp)<=3){
				X_spl<-as.matrix(cbind(bs_ps,bs_y1s))
				X_spl<-quantile_normalize_bart(X_spl)
				check_singularity <- sapply(as.data.frame(X_spl), function(x) length(unique(x)) == 1)
				if(sum(check_singularity)>0){
					X_spl[, check_singularity] <- rnorm(dim(X_spl)[1], 0, 0.00000001)
				}
				X_spl<-as.matrix(cbind(1,X_spl))
			}else{
				X_spl<-as.matrix(cbind(bs_ps,bs_y1s,datov_sp[,4:ncol(datov_sp)]))
				X_spl<-quantile_normalize_bart(X_spl)
				check_singularity <- sapply(as.data.frame(X_spl), function(x) length(unique(x)) == 1)
				if(sum(check_singularity)>0){
					X_spl[, check_singularity] <- rnorm(dim(X_spl)[1], 0, 0.00000001)
				}
				X_spl<-as.matrix(cbind(1,X_spl))
			}
			
			## sample the betas ##
			Vbeta<-solve(Sigma0_inv+((1/sigma_spl1)*t(X_spl)%*%X_spl))
			Ebeta<-Vbeta%*%(Sigma0_inv%*%mu0+((1/sigma_spl1)*t(X_spl)%*%delta_sp))
			beta1<-matrix(mvrnorm(n=1,mu=Ebeta,Sigma=Vbeta),nrow=p_spl1,ncol=1)
			
			## sample the sigma_spls ##
			a<-a0+(nrow(datov_sp)/2)
			b<-b0+((1/2)*sum((delta_sp-(X_spl%*%beta1))^2))
			sigma_spl1<-rinvgamma(n=1,shape=a,scale=b)
			
			######################################################################
			## 3. draw from the posterior predictive for the non-overlap region ##
			######################################################################
			bs_ps_star<-rcspline.eval(datno1$ps,knots=quantile(datov_sp$ps,probs=probs.rcs.ps),inclx = T)
			bs_y1s_star<-rcspline.eval(datno1$Yobs,knots=quantile(Y1s,probs=probs.rcs.Y),inclx=T)
			
			
			if (ncol(datno)<=3){
				X_spl_star<-as.matrix(cbind(bs_ps_star,bs_y1s_star))
				X_spl_star<-quantile_normalize_bart(X_spl_star)
				check_singularity <- sapply(as.data.frame(X_spl_star), function(x) length(unique(x)) == 1)
				if(sum(check_singularity)>0){
					X_spl_star[, check_singularity] <- rnorm(dim(X_spl_star)[1], 0, 0.00000001)
				}
				X_spl_star<-as.matrix(cbind(1,X_spl_star))
			}else{
				X_spl_star<-as.matrix(cbind(bs_ps_star,bs_y1s_star,datno1[,4:ncol(datno1)]))
				X_spl_star<-quantile_normalize_bart(X_spl_star)
				check_singularity <- sapply(as.data.frame(X_spl_star), function(x) length(unique(x)) == 1)
				if(sum(check_singularity)>0){
					X_spl_star[, check_singularity] <- rnorm(dim(X_spl_star)[1], 0, 0.00000001)
				}
				X_spl_star<-as.matrix(cbind(1,X_spl_star))
			}
			Eppd<-X_spl_star%*%beta1
			delta_star[which(datno$x==1)]<-rnorm(n=nrow(datno1),mean=Eppd,sd=sqrt(sigma_spl1+ROdist[which(datno$x==1)]*var_infl*(max(delta)-min(delta))))
		}
		
		## spline for untreated units (E=0) in the RN ##
		if (nrow(datno0)>0){
			
			####################################
			## 2. run the spline sampler once ##
			####################################
			
			bs_y0s<-rcspline.eval(Y0s,knots=quantile(Y0s,probs=probs.rcs.Y),inclx=T)
			if (ncol(datov_sp)<=3){
				X_spl<-as.matrix(cbind(bs_ps,bs_y0s))
				X_spl<-quantile_normalize_bart(X_spl)
				check_singularity <- sapply(as.data.frame(X_spl), function(x) length(unique(x)) == 1)
				if(sum(check_singularity)>0){
					X_spl[, check_singularity] <- rnorm(dim(X_spl)[1], 0, 0.00000001)
				}
				X_spl<-as.matrix(cbind(1,X_spl))
				}
			else{
				X_spl<-as.matrix(cbind(bs_ps,bs_y0s,datov_sp[,4:ncol(datov_sp)]))
				X_spl<-quantile_normalize_bart(X_spl)
				check_singularity <- sapply(as.data.frame(X_spl), function(x) length(unique(x)) == 1)
				if(sum(check_singularity)>0){
					X_spl[, check_singularity] <- rnorm(dim(X_spl)[1], 0, 0.00000001)
				}
				X_spl<-as.matrix(cbind(1,X_spl))
				}
			
			## sample the betas ##
			Vbeta<-solve(Sigma0_inv+((1/sigma_spl0)*t(X_spl)%*%X_spl))
			Ebeta<-Vbeta%*%(Sigma0_inv%*%mu0+((1/sigma_spl0)*t(X_spl)%*%delta_sp))
			beta0<-matrix(mvrnorm(n=1,mu=Ebeta,Sigma=Vbeta),nrow=p_spl0,ncol=1)
			
			## sample the sigma_spls ##
			a<-a0+(nrow(datov_sp)/2)
			b<-b0+((1/2)*sum((delta_sp-(X_spl%*%beta0))^2))
			sigma_spl0<-rinvgamma(n=1,shape=a,scale=b)
			
			######################################################################
			## 3. draw from the posterior predictive for the non-overlap region ##
			######################################################################
			bs_ps_star<-rcspline.eval(datno0$ps,knots=quantile(datov_sp$ps,probs=probs.rcs.ps),inclx = T)
			bs_y0s_star<-rcspline.eval(datno0$Yobs,knots=quantile(Y0s,probs=probs.rcs.Y),inclx=T)
			if (ncol(datno)<=3){
				X_spl_star<-as.matrix(cbind(bs_ps_star,bs_y0s_star))
				X_spl_star<-quantile_normalize_bart(X_spl_star)
				check_singularity <- sapply(as.data.frame(X_spl_star), function(x) length(unique(x)) == 1)
				if(sum(check_singularity)>0){
					X_spl_star[, check_singularity] <- rnorm(dim(X_spl_star)[1], 0, 0.00000001)
				}
				X_spl_star<-as.matrix(cbind(1,X_spl_star))
			}
			else{
				X_spl_star<-as.matrix(cbind(bs_ps_star,bs_y0s_star,datno0[,4:ncol(datno0)]))
				X_spl_star<-quantile_normalize_bart(X_spl_star)
				check_singularity <- sapply(as.data.frame(X_spl_star), function(x) length(unique(x)) == 1)
				if(sum(check_singularity)>0){
					X_spl_star[, check_singularity] <- rnorm(dim(X_spl_star)[1], 0, 0.00000001)
				}
				X_spl_star<-as.matrix(cbind(1,X_spl_star))
			}
			Eppd<-X_spl_star%*%beta0
			delta_star[which(datno$x==0)]<-rnorm(n=nrow(datno0),mean=Eppd,sd=sqrt(sigma_spl0+ROdist[which(datno$x==0)]*10*(max(delta)-min(delta))))
		}
		
		
		###############################################
		## 4. if we're past burn-in, save the output ##
		###############################################
		if (i>nburn){
			delta_star_save[(i-nburn),which(RO==1)]<-c(delta)
			delta_star_save[(i-nburn),which(RO==0)]<-delta_star
		}
		
	}
	
	## use the Bayesian bootstrap to get ACE posterior ##
	ace_pd<-apply(delta_star_save,1,aceBB)
	
	return(list(ace_pd,apply(delta_star_save,2,mean),apply(delta_star_save,2,quantile,probs=.025),apply(delta_star_save,2,quantile,probs=.975)))
	
}
