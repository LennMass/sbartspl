####
# Estimate ATE effects 
# DGP based Nethery et al. (2019)
####


####
# Set simulation params ----
####

simnum_vec <- seq(1237, 1248, 1) # seed

# parameters
N = 2000 # 500 and 2000 used in paper
ncovs = 50 # 10, 25 and 50 used in paper
reps=10



####
# Start simulation loop over seeds----
####


for(j in 1:length(simnum_vec)){
	
	simnum <- simnum_vec[j]
	set.seed(simnum) # 10 seeds, 10 datasets per seed, gives 100 simulation datasets 
	
	filename_specs <- paste("n", N, 
													"_ncovs", ncovs, 
													"_seed", simnum,
													sep="")
	
	
	# save metrics
	method_names <- c("UGR", "UBART", "USBART", "TGR", "TBART", "TSBART", "XBCF", "BLR", "BARTSPL", "SBARTSPL")
	err_res <- matrix(NA, nrow = reps,  ncol=length(method_names))
	colnames(err_res) <- method_names
	
	cicoverage_res <- matrix(NA, nrow = reps,  ncol=length(method_names))
	colnames(cicoverage_res) <- method_names
	
	ciwidth_res <- matrix(NA, nrow = reps,  ncol=length(method_names))
	colnames(ciwidth_res) <- method_names
	
	
	####
	# Start simulation loop for 10 sets per seed----
	####
	
	
	
	for( i in 1:reps){
		
		####
		# Load and prepare dataset ----
		####
		
		simdat <- readRDS(paste("00_data/01_simulations/01_main/",
														filename_specs, 
														"_simid", i, 
														".rds",
														sep="")
		)
		
		# order categorical values to the end
		cat_vars <- simdat$untrimmed_dat[, c("X1", "X3", "X5", "X7", "X9")]
		dataset_rest <- simdat$untrimmed_dat[, !names(simdat$untrimmed_dat) %in% 
																				 	c("X1", "X3", "X5", "X7", "X9")] 
		test2<- cbind(dataset_rest, cat_vars)
		RO <- simdat$RO
		
		# remove unneccesary placeholder
		rm(cat_vars)
		rm(dataset_rest)
		
		# needed to create test set for BART, SoftBart
		test3<- test2
		test3$x<- 1-test3$x
		
		####
		# Start estimation ----
		####
		
		####
		## untrimmed BART ----
		####
		
		
		utrim_bartfit<-bartalone(xtr=test2[,2:ncol(test2)],
														 ytr=test2$Yobs,
														 xte=test3[,2:ncol(test2)])
		
		ubart <- list(ace_pd = utrim_bartfit$ace_pd,
									ate = utrim_bartfit$ace_pd %>% mean(),
									atel = utrim_bartfit$ace_pd %>% quantile(0.025),
									ateu = utrim_bartfit$ace_pd %>% quantile(0.975))
		
		err_res[i, "UBART"] <- ubart$ate - simdat$ate_true
		cicoverage_res[i, "UBART"] <- ifelse( (ubart$atel <= simdat$ate_true) && (ubart$ateu >= simdat$ate_true), 1, 0)
		ciwidth_res[i, "UBART"] <- ubart$ateu - ubart$atel
		
		
		
		####
		## untrimmed SoftBART ----
		####
		
		
		utrim_softbartfit<-softbartalone(xtr=test2[,2:ncol(test2)],
																		 ytr=test2$Yobs,
																		 xte=test3[,2:ncol(test2)])
		
		usoftbart <- list(ace_pd = utrim_softbartfit$ace_pd,
											ate= utrim_softbartfit$ace_pd %>% mean(),
											atel = utrim_softbartfit$ace_pd %>% quantile(0.025),
											ateu = utrim_softbartfit$ace_pd %>% quantile(0.975))
		
		err_res[i, "USBART"] <- usoftbart$ate - simdat$ate_true
		cicoverage_res[i, "USBART"] <- ifelse( (usoftbart$atel <= simdat$ate_true) && (usoftbart$ateu >= simdat$ate_true), 1, 0)
		ciwidth_res[i, "USBART"] <- usoftbart$ateu - usoftbart$atel
		
		####
		## trimmed BART ----
		####
		
		
		trim_bartfit<-bartalone(xtr=test2[which(RO==1), 2:ncol(test2)],
														ytr=test2$Yobs[which(RO==1)],
														xte=test3[which(RO==1), 2:ncol(test2)])
		
		tbart <- list(ace_pd = trim_bartfit$ace_pd,
									ate = trim_bartfit$ace_pd %>% mean(),
									atel = trim_bartfit$ace_pd %>% quantile(0.025),
									ateu = trim_bartfit$ace_pd %>% quantile(0.975))
		
		err_res[i, "TBART"] <- tbart$ate - simdat$ate_true
		cicoverage_res[i, "TBART"] <- ifelse( (tbart$atel <= simdat$ate_true) && (tbart$ateu >= simdat$ate_true), 1, 0)
		ciwidth_res[i, "TBART"] <- tbart$ateu - tbart$atel
		
		
		####
		## trimmed SoftBART ----
		####
		
		
		trim_softbartfit<-softbartalone(xtr=test2[which(RO==1), 2:ncol(test2)],
																		ytr=test2$Yobs[which(RO==1)],
																		xte=test3[which(RO==1), 2:ncol(test2)])
		
		tsoftbart <- list(ace_pd = trim_softbartfit$ace_pd,
											ate = trim_softbartfit$ace_pd %>% mean(),
											atel = trim_softbartfit$ace_pd %>% quantile(0.025),
											ateu = trim_softbartfit$ace_pd %>% quantile(0.975))
		
		err_res[i, "TSBART"] <- tsoftbart$ate - simdat$ate_true
		cicoverage_res[i, "TSBART"] <- ifelse( (tsoftbart$atel <= simdat$ate_true) && (tsoftbart$ateu >= simdat$ate_true), 1, 0)
		ciwidth_res[i, "TSBART"] <- tsoftbart$ateu - tsoftbart$atel
		
		
		####
		## BART+SPL ----
		####
		
		bartspl_model <-bartspl(datall=test2,RO=RO)
		
		bartspl.fit <- list(ace_pd = bartspl_model$ace_pd,
												ate = bartspl_model$ace_pd %>% mean(),
												atel = bartspl_model$ace_pd %>% quantile(0.025),
												ateu = bartspl_model$ace_pd %>% quantile(0.975))
		
		err_res[i, "BARTSPL"] <- bartspl.fit$ate - simdat$ate_true
		cicoverage_res[i, "BARTSPL"] <- ifelse( (bartspl.fit$atel <= simdat$ate_true) && (bartspl.fit$ateu >= simdat$ate_true), 1, 0)
		ciwidth_res[i, "BARTSPL"] <- bartspl.fit$ateu - bartspl.fit$atel
		
		####
		## SBART+SPL ----
		####
		
		sbartspl_model<-bartspl(datall=test2,
														RO=RO,
														nburn = 2500,
														nsamp = 1250,
														bart_model = "sbart",
														softbart_sampler = TRUE,
														probs.rcs.ps=c(.25,.5,0.75),
														probs.rcs.Y=c(.25,.5,0.75))
		
		sbartspl.fit <- list(ace_pd = sbartspl_model$ace_pd,
												 ate = sbartspl_model$ace_pd %>% mean(),
												 atel = sbartspl_model$ace_pd %>% quantile(0.025),
												 ateu = sbartspl_model$ace_pd %>% quantile(0.975))
		
		err_res[i, "SBARTSPL"] <- sbartspl.fit$ate - simdat$ate_true
		cicoverage_res[i, "SBARTSPL"] <- ifelse( (sbartspl.fit$atel <= simdat$ate_true) && (sbartspl.fit$ateu >= simdat$ate_true), 1, 0)
		ciwidth_res[i, "SBARTSPL"] <- sbartspl.fit$ateu - sbartspl.fit$atel
		
		####
		## trimmed GR ----
		####
		
		ce_tgr_cor<-gr(Y=test2$Yobs[which(RO==1)],
									 trt=test2$x[which(RO==1)],
									 ps=test2$ps_mis[which(RO==1)],
									 X=test2[which(RO==1),4:ncol(test2)],
									 M=500,qps=quantile(test2$ps_mis[which(RO==1)],probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))
		
		tgr <-list(ace_pd = ce_tgr_cor$ace,
							 ate = ce_tgr_cor$ace %>% mean(),
							 atel = ce_tgr_cor$ace %>% quantile(0.025),
							 ateu = ce_tgr_cor$ace %>% quantile(0.975))
		
		err_res[i, "TGR"] <- tgr$ate - simdat$ate_true
		cicoverage_res[i, "TGR"] <- ifelse( (tgr$atel <= simdat$ate_true) && (tgr$ateu >= simdat$ate_true), 1, 0)
		ciwidth_res[i, "TGR"] <- tgr$ateu - tgr$atel
		
		
		
		####
		## untrimmed GR ----
		####
		
		ugr_fit <-gr(Y=test2$Yobs,
								 trt=test2$x,
								 ps=test2$ps_mis,
								 X=test2[,4:ncol(test2)],
								 M=500,qps=quantile(test2$ps_mis,probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))
		
		ugr <- list(ace_pd = ugr_fit$ace %>% as.vector(),
								ate= ugr_fit$ace %>%  mean(),
								atel = ugr_fit$ace %>% quantile(.025),
								ateu = ugr_fit$ace %>% quantile(.975)
		)
		
		err_res[i, "UGR"] <- ugr$ate - simdat$ate_true
		cicoverage_res[i, "UGR"] <- ifelse( (ugr$atel <= simdat$ate_true) && (ugr$ateu >= simdat$ate_true), 1, 0)
		ciwidth_res[i, "UGR"] <- ugr$ateu - ugr$atel
		
		####
		## XBCF ----
		####
		
		xbcf.model <- XBCF::XBCF(y=test2$Yobs,
														 z=test2$x,
														 x_con=as.matrix(test2[,4:ncol(test2)]),
														 x_mod=as.matrix(test2[,4:ncol(test2)]),
														 pihat = test2$ps_mis,
														 pcat_con = 5, # five true binary covs
														 pcat_mod = 5, # five true binary covs
														 n_trees_mod = 20, num_sweeps = 100, Nmin = 20)
		xbcf.tau <- xbcf.model$tauhats.adjusted
		xbcf.ace_pd <- apply(xbcf.tau, 2, aceBB)
		xbcf <- list(ace_pd = xbcf.ace_pd,
								 ate = xbcf.ace_pd %>% mean(),
								 atel = xbcf.ace_pd %>% quantile(0.025),
								 ateu = xbcf.ace_pd %>% quantile(0.975))
		
		
		err_res[i, "XBCF"] <- xbcf$ate - simdat$ate_true
		cicoverage_res[i, "XBCF"] <- ifelse( (xbcf$atel <= simdat$ate_true) && (xbcf$ateu >= simdat$ate_true), 1, 0)
		ciwidth_res[i, "XBCF"] <- xbcf$ateu - xbcf$atel
		
		
		####
		## Bayesian Linear Regression (BLR) ----
		####
		
		# standardize covariate variables and outcome
		test2_scale <- test2
		cat_vars <- colnames(test2)[c(2, (ncol(test2)-4):ncol(test2))]
		cont_vars <- colnames(test2)[!(colnames(test2) %in% cat_vars)]
		test2_scale[cont_vars] <- scale(test2[cont_vars])
		test2_scale$Yobs <- scale(test2$Yobs)
		
		
		form_gam <- brms::bf(
			paste("Yobs ~ x",
						#paste(sprintf("s(%s, k = 4)", cont_vars[-1]), collapse = " + "),
						paste(sprintf(cont_vars[-1]), collapse = " + "),
						paste(sprintf(cat_vars[-1]), collapse = " + "),
						sep = " + "))
		
		
		priors_gam <- c(
			brms::prior(normal(0, 1.5), class = "b", coef = "x"),
			brms::prior(normal(0, 1), class = "b"),
			brms::prior(student_t(3, 0, 2.5), class = "sigma")
		)
		
		
		fit_ugam <- brms::brm(
			form_gam,
			data = test2_scale,
			prior = priors_gam,
			chains = 2,
			iter = 2000,
			control = list(adapt_delta = 0.9)
		)
		
		
		ugam_m1 <- brms::posterior_epred(fit_ugam, newdata = transform(test2_scale, x = 1))
		ugam_m0 <- brms::posterior_epred(fit_ugam, newdata = transform(test2_scale, x = 0))
		
		# ATE posterior draws, back-transformed on original scale
		ugam_ATE_post <- rowMeans(ugam_m1 - ugam_m0) * attr(test2_scale$Yobs, "scaled:scale")
		ugam <- list(ace_pd=ugam_ATE_post,
								 ate = ugam_ATE_post %>% mean(), 
								 atel = ugam_ATE_post %>% quantile(0.025), 
								 ateu = ugam_ATE_post %>% quantile(0.975) 
		)
		
		err_res[i, "BLR"] <- ugam$ate - simdat$ate_true
		cicoverage_res[i, "BLR"] <- ifelse( (ugam$atel <= simdat$ate_true) && (ugam$ateu >= simdat$ate_true), 1, 0)
		ciwidth_res[i, "BLR"] <- ugam$ateu - ugam$atel
		
		print(paste("Finished sim run i =", i, "in seed run j =", j, sep=" "))
		
	} # end of loop over sets per seed
	
	
	
	
	
	####
	# Save results ----
	####
	
	
	saveRDS(err_res, 
					paste("02_results/01_simulations/01_main/",
								"err_results_", 
								filename_specs,
								"_reps", reps, 
								".rds",
								sep="")
	)
	
	
	
	saveRDS(cicoverage_res, 
					paste("02_results/01_simulations/01_main/",
								"cicoverage_results_", 
								filename_specs,
								"_reps", reps, 
								".rds",
								sep="")
	)
	
	saveRDS(ciwidth_res, 
					paste("02_results/01_simulations/01_main/",
								"ciwidth_results_", 
								filename_specs,
								"_reps", reps, 
								".rds",
								sep="")
	)
	
	print(paste("Finished seed run j =", j, sep=" "))
	
	
} # end of loop over seeds

