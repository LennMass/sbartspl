# BART+SPL / SBART+SPL for Continuous Outcomes
# ------------------------------------------------------
# Extends Nethery et al. (2019)'s BART+SPL to allow SoftBART, DART, and
# Smooth BART in the region of overlap (RO), with a restricted cubic
# spline (RCS) extrapolation in the region of non-overlap (RN). Returns
# the posterior of the average causal effect (ACE) and per-unit ITE
# summaries.
#
# Input data layout (datall):
#   col 1         : observed outcome (Yobs)
#   col 2         : binary exposure indicator (x)
#   col 3         : propensity score (ps) - basis of the RO/RN split
#   col 4 onwards : additional covariates included in both fits
# RO is the 0/1 overlap indicator from pw_overlap(), aligned with datall.
#
# Args:
#   datall           : data.frame in the layout above
#   RO               : 0/1 overlap indicator, length nrow(datall)
#   nburn, nsamp     : MCMC burn-in and stored draws (default 10000, 5000)
#   bart_model       : "dbarts" | "sbart" | "dart" | "smooth_bart"
#   softbart_sampler : TRUE for any SoftBart-family model
#                      (sbart / dart / smooth_bart)
#   n_forest         : number of SoftBart chains run in parallel (default 4)
#   var_infl         : variance inflation for the RN spline extrapolation
#                      when softbart_sampler = TRUE (default 0)
#   tau_rate_dart    : only used for "bart_model=dart". large bandwidth control
#											 param for splitting if using DART instead of SBART
#   probs.rcs.ps     : knot quantiles for the PS spline basis 
#   probs.rcs.Y      : knot quantiles for the Y  spline basis 
#
# Depends on the project-internal helpers preprocess_df(),
# quantile_normalize_bart(), normalize_bart(), unnormalize_bart(),
# and aceBB(), and on the packages dbarts, SoftBart, Hmisc, MASS, invgamma.
#
# Returns a list with:
#   ace_pd          : posterior draws of the ACE (length nsamp)
#   delta_star_save : nsamp x n matrix of ITE draws
#   ice_mean/median/sd/lw/up : column-wise summaries of delta_star_save
#                              (lw, up = 2.5% and 97.5% quantiles)
#   share_NA        : share of NA cells in delta_star_save

bartspl<-function(datall,RO,
                  nburn=10000, 
                  nsamp=5000, 
                  bart_model="dbarts", 
                  softbart_sampler=FALSE, 
                  n_forest=4, 
                  var_infl = 0, 
                  tau_rate_dart=1e+18,  
									probs.rcs.ps=c(.1,.25,.5,.75,.9), 
									probs.rcs.Y=c(.2,.4,.6,.8) 
                  ){
  
  
  #### prepare datasets ####
  
  names(datall)[1:3]<-c('Yobs','x','ps')
  if (ncol(datall)>3) names(datall)[4:ncol(datall)]<-paste('u',1:(ncol(datall)-3),sep='')
  
  if(softbart_sampler==T){
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
  }
  
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
  
  # MakeForest function for SBART/DART/SmoothBART models
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
  											 opts=SoftBart::Opts(num_burn=0, # because we do burn-in and saving later on our own by nburn, nsamp
  											 										num_save=1, # because we do burn-in and saving later on our own by nburn, nsamp
  											 										update_tau=update_tau_sb,
  											 										update_tau_mean=update_tau_mean_sb,
  											 										update_s = update_s_sb, 
  											 										update_alpha = update_alpha_sb 
  	)) 
  }
  
  #### initialize samplers ####
  if(bart_model=="dbarts"){
    ## initialize the dbarts sampler to implement BART in the range of overlap ##
    dbfit<-dbarts::dbarts(formula=Yobs~.,data=datov,test=datovtest[,2:ncol(datovtest)], verbose = F)
  }else if(bart_model=="sbart"){
    ## initialize MakeForest sampler to implement SBART in the range of overlap
  	
    ##initialize chains (4 in total)
    sb_chain1 <- sb_forest()
    sb_chain2 <- sb_forest()
    sb_chain3 <- sb_forest()
    sb_chain4 <- sb_forest()
    
    
  }else if(bart_model=="dart"){
    ## initialize MakeForest sampler to implement DART in the range of overlap
    
    ##initialize chains (4 in total)
    sb_chain1 <- sb_forest(tau_rate_sb=tau_rate_dart, update_tau_sb = FALSE, update_tau_mean_sb = FALSE)
    sb_chain2 <- sb_forest(tau_rate_sb=tau_rate_dart, update_tau_sb = FALSE, update_tau_mean_sb = FALSE)
    sb_chain3 <- sb_forest(tau_rate_sb=tau_rate_dart, update_tau_sb = FALSE, update_tau_mean_sb = FALSE)
    sb_chain4 <- sb_forest(tau_rate_sb=tau_rate_dart, update_tau_sb = FALSE, update_tau_mean_sb = FALSE)
    
  }else if(bart_model=="smooth_bart"){
  	## initialize MakeForest sampler to implement smooth BART in the range of overlap
  	
  	##initialize chains (4 in total)
  	sb_chain1 <- sb_forest(alpha_sb=ncol(datov_scale[, !(colnames(datov_scale) == "Yobs")]), update_s_sb = FALSE, update_alpha_sb = FALSE)
  	sb_chain2 <- sb_forest(alpha_sb=ncol(datov_scale[, !(colnames(datov_scale) == "Yobs")]), update_s_sb = FALSE, update_alpha_sb = FALSE)
  	sb_chain3 <- sb_forest(alpha_sb=ncol(datov_scale[, !(colnames(datov_scale) == "Yobs")]), update_s_sb = FALSE, update_alpha_sb = FALSE)
  	sb_chain4 <- sb_forest(alpha_sb=ncol(datov_scale[, !(colnames(datov_scale) == "Yobs")]), update_s_sb = FALSE, update_alpha_sb = FALSE)
  }else{
    stop("Specify if bart_model from c('dbarts', 'sbart', 'dart', 'smooth_bart)")
  }

  
  
  
  #### initialize params and matrices ####
  
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
  delta_star_burn<-matrix(NA,nrow=nburn,ncol=nrow(datall))
  
  
  # for softbart_sampler
  if(softbart_sampler==T){
    s_temp_test <- matrix(NaN, ncol=n_forest, nrow=nrow(datovtest)) 
    s_temp_sigma <- matrix(NaN, ncol=n_forest, nrow=1)
    s_temp_pred <- matrix(NaN, ncol=n_forest, nrow=nrow(datovtest)) 
  }
  
  
  #### run the sampler ####
  for (i in 1:(nburn+nsamp)){
    
    ##### 1. run the dbarts sampler once and take a sample from the BART ppd ####

    if(softbart_sampler==FALSE){
      temp<-dbfit$run(numBurnIn=0,numSamples=1)  
      ppd_test<-rnorm(n=length(temp$test),mean=temp$test,sd=temp$sigma) 
    }else if(softbart_sampler==TRUE){
      # for softbart_sampler
    	# chain 1
      s_temp_test[,1]<- sb_chain1$do_gibbs(X=data.matrix(datov_scale[, !(colnames(datov_scale) == "Yobs")]), 
                                                    Y=data.matrix(datov_scale$Yobs),
                                                    X_test=data.matrix(datovtest_scale[,!(colnames(datovtest_scale) == "Yobs")]), 
                                                    i=1) %>%
          t()
      s_temp_sigma[,1] <- sb_chain1$get_sigma() # get_sigma() gives error standard deviation of the forest
      # chain 2
      s_temp_test[,2]<- sb_chain2$do_gibbs(X=data.matrix(datov_scale[, !(colnames(datov_scale) == "Yobs")]), 
      																		 Y=data.matrix(datov_scale$Yobs),
      																		 X_test=data.matrix(datovtest_scale[,!(colnames(datovtest_scale) == "Yobs")]), 
      																		 i=1) %>%
      	t()
      s_temp_sigma[,2] <- sb_chain2$get_sigma() # get_sigma() gives error standard deviation of the forest
      # chain 3
      s_temp_test[,3]<- sb_chain3$do_gibbs(X=data.matrix(datov_scale[, !(colnames(datov_scale) == "Yobs")]), 
      																		 Y=data.matrix(datov_scale$Yobs),
      																		 X_test=data.matrix(datovtest_scale[,!(colnames(datovtest_scale) == "Yobs")]), 
      																		 i=1) %>%
      	t()
      s_temp_sigma[,3] <- sb_chain3$get_sigma() # get_sigma() gives error standard deviation of the forest
      # chain 4
      s_temp_test[,4]<- sb_chain4$do_gibbs(X=data.matrix(datov_scale[, !(colnames(datov_scale) == "Yobs")]), 
      																		 Y=data.matrix(datov_scale$Yobs),
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
    }else{
      stop("Specify if softbart_sampler=T/F.")
    }
    
    ###### form the predicted individual causal effects in the RO from the BART predictions ####
    delta[which(datov$x==1),]<-datov$Yobs[which(datov$x==1)]-ppd_test[which(datov$x==1)]
    delta[which(datov$x==0),]<-ppd_test[which(datov$x==0)]-datov$Yobs[which(datov$x==0)]
    delta_sp<-delta[sp_ind,]
    
    ###### update Y1s #### 
    if(softbart_sampler==F){
      foosp<-temp$test[sp_ind]
    }else if(softbart_sampler==T){
      foosp<-s_temp_test[sp_ind]
    }
    Y1s[which(datov_sp$x==0)]<-foosp[which(datov_sp$x==0)]
    Y0s[which(datov_sp$x==1)]<-foosp[which(datov_sp$x==1)]
    
    ###### fit spline to the estimated causal effects in the RO from BART ####
    delta_star<-rep(NA,sum(RO==0))
    	
    bs_ps<-Hmisc::rcspline.eval(datov_sp$ps,
    															knots=quantile(datov_sp$ps,
    																						 probs=probs.rcs.ps),
    															inclx = T)
    
    
    ##### spline for treated units (E=1) in the RN ####
    if (nrow(datno1)>0){
      
      ##### 2. run the spline sampler once ####
    	bs_y1s<-Hmisc::rcspline.eval(Y1s,
    																 knots=quantile(Y1s,
    																 							 probs=probs.rcs.Y),
    																 inclx = T)
    
      if (ncol(datov_sp)<=3){
        X_spl<-as.matrix(cbind(1,bs_ps,bs_y1s))
      }
      else{
        X_spl<-as.matrix(cbind(1,bs_ps,bs_y1s,datov_sp[,4:ncol(datov_sp)]))
      }
      
      ###### sample the betas ## 
      Vbeta<-solve(Sigma0_inv+((1/sigma_spl1)*t(X_spl)%*%X_spl))
      Ebeta<-Vbeta%*%(Sigma0_inv%*%mu0+((1/sigma_spl1)*t(X_spl)%*%delta_sp))
      beta1<-matrix(MASS::mvrnorm(n=1,mu=Ebeta,Sigma=Vbeta),nrow=p_spl1,ncol=1)
      
      ###### sample the sigma_spls ## 
      a<-a0+(nrow(datov_sp)/2)
      b<-b0+((1/2)*sum((delta_sp-(X_spl%*%beta1))^2))
      sigma_spl1<-invgamma::rinvgamma(n=1,shape=a,scale=b)
      
      ##### 3. draw from the posterior predictive for the non-overlap region ####
      
      bs_ps_star<-Hmisc::rcspline.eval(datno1$ps,knots=quantile(datov_sp$ps,probs=probs.rcs.ps),inclx = T)
      bs_y1s_star<-Hmisc::rcspline.eval(datno1$Yobs,knots=quantile(Y1s,probs=probs.rcs.Y),inclx=T)
      
      if (ncol(datno)<=3){
        X_spl_star<-as.matrix(cbind(1,bs_ps_star,bs_y1s_star))
      }
      else{
        X_spl_star<-as.matrix(cbind(1,bs_ps_star,bs_y1s_star,datno1[,4:ncol(datno1)]))
      }
      Eppd<-X_spl_star%*%beta1
      if(softbart_sampler==FALSE){      
        delta_star[which(datno$x==1)]<-rnorm(n=nrow(datno1),mean=Eppd,sd=sqrt(sigma_spl1+ROdist[which(datno$x==1)]*10*(max(delta)-min(delta))))
      }
      if(softbart_sampler==TRUE){      
        delta_star[which(datno$x==1)]<-rnorm(n=nrow(datno1),mean=Eppd,sd=sqrt(sigma_spl1+ROdist[which(datno$x==1)]*var_infl*(max(delta)-min(delta))))
      }
    }
    
    #### spline for untreated units (E=0) in the RN ####
    if (nrow(datno0)>0){
      
      ###### 2. run the spline sampler once ####
    	
    	bs_y0s<-Hmisc::rcspline.eval(Y0s,knots=quantile(Y0s,probs=probs.rcs.Y),inclx=T)
    	
      if (ncol(datov_sp)<=3){
        X_spl<-as.matrix(cbind(1,bs_ps,bs_y0s))
      }
      else{
        X_spl<-as.matrix(cbind(1,bs_ps,bs_y0s,datov_sp[,4:ncol(datov_sp)]))
      }
      
      ##### sample the betas ####
      Vbeta<-solve(Sigma0_inv+((1/sigma_spl0)*t(X_spl)%*%X_spl))
      Ebeta<-Vbeta%*%(Sigma0_inv%*%mu0+((1/sigma_spl0)*t(X_spl)%*%delta_sp))
      beta0<-matrix(MASS::mvrnorm(n=1,mu=Ebeta,Sigma=Vbeta),nrow=p_spl0,ncol=1)
      
      ##### sample the sigma_spls ####
      a<-a0+(nrow(datov_sp)/2)
      b<-b0+((1/2)*sum((delta_sp-(X_spl%*%beta0))^2))
      sigma_spl0<-invgamma::rinvgamma(n=1,shape=a,scale=b)
      
      ##### 3. draw from the posterior predictive for the non-overlap region ####

      bs_ps_star<-Hmisc::rcspline.eval(datno0$ps,knots=quantile(datov_sp$ps,probs=probs.rcs.ps),inclx = T)
      bs_y0s_star<-Hmisc::rcspline.eval(datno0$Yobs,knots=quantile(Y0s,probs=probs.rcs.Y),inclx=T)
      
      
      if (ncol(datno)<=3){
        X_spl_star<-as.matrix(cbind(1,bs_ps_star,bs_y0s_star))
      }
      else{
        X_spl_star<-as.matrix(cbind(1,bs_ps_star,bs_y0s_star,datno0[,4:ncol(datno0)]))
      }
      Eppd<-X_spl_star%*%beta0
      if(softbart_sampler==FALSE){
        delta_star[which(datno$x==0)]<-rnorm(n=nrow(datno0),mean=Eppd,sd=sqrt(sigma_spl0+ROdist[which(datno$x==0)]*10*(max(delta)-min(delta))))
      }
      if(softbart_sampler==TRUE){
        delta_star[which(datno$x==0)]<-rnorm(n=nrow(datno0),mean=Eppd,sd=sqrt(sigma_spl0+ROdist[which(datno$x==0)]*var_infl*(max(delta)-min(delta))))
      }
    }
    
    
    #### 4. if we're past burn-in, save the output ####
    
    if (i<=nburn){
    	delta_star_burn[(i),which(RO==1)]<-c(delta)
    	delta_star_burn[(i),which(RO==0)]<-delta_star
    }
    
    if (i>nburn){
      delta_star_save[(i-nburn),which(RO==1)]<-c(delta)
      delta_star_save[(i-nburn),which(RO==0)]<-delta_star
    }
    
  }
  
  ## use the Bayesian bootstrap to get ACE posterior ##
  ace_pd<-apply(delta_star_save,1,aceBB)
  
  
  
  return(list(ace_pd=ace_pd,
  						delta_star_save=delta_star_save, 
              ice_mean=apply(delta_star_save,2,mean, na.rm=T), 
  						ice_median=apply(delta_star_save,2,median, na.rm=T), 
  						ice_sd=apply(delta_star_save,2,sd, na.rm=T),
              ice_lw=apply(delta_star_save,2,quantile,probs=.025, na.rm=T),
              ice_up=apply(delta_star_save,2,quantile,probs=.975, na.rm=T),
              share_NA=sum(is.na(delta_star_save)) / (dim(delta_star_save)[1]*dim(delta_star_save)[2])))
  
}

