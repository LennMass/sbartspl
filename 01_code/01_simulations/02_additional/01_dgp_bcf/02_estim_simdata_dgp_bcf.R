####
# Estimate ATE effects 
# DGP based Hahn et al. (2020), Krantsevich et al. (2022), Wang et al. (2024)
####


####
# Set simulation params ----
####

n = 500 
tau_setting="heterogeneous" # c("homogeneous", "heterogeneous")
p_cont_noise = 5 # c(5, 10)
p_cat_noise = 5 # c(5, 10)
reps <- 20

simnum <- 123
set.seed(simnum)


filename_specs <- paste("n", n, 
												"_tausetting", tau_setting, 
												"_pnoise", p_cont_noise+p_cat_noise,
												"_seed", simnum,
												sep="")


# save metrics
method_names <- c("UGR", "UBART", "USBART", "TGR", "TBART", "TSBART", "XBCF", "BARTSPL", "SBARTSPL")
err_res <- matrix(NA, nrow = reps,  ncol=length(method_names))
colnames(err_res) <- method_names

cicoverage_res <- matrix(NA, nrow = reps,  ncol=length(method_names))
colnames(cicoverage_res) <- method_names

ciwidth_res <- matrix(NA, nrow = reps,  ncol=length(method_names))
colnames(ciwidth_res) <- method_names


####
# Start simulation loop ----
####



for( i in 1:reps){
	
	####
	# Load and prepare dataset ----
	####
	
	simdat <- readRDS(paste("00_data/01_simulations/02_additional/01_dgp_bcf/",
																									filename_specs, 
																									"_simid", i, 
																									".rds",
																									sep="")
										)
	
	
	test2<- simdat$untrimmed_dat
	RO <- simdat$RO
	
	# needed to create test set for BART, SoftBart
	test3<- simdat$untrimmed_dat
	test3$z<- 1-test3$z
	
	####
	# Start estimation ----
	####
	
	####
	## untrimmed BART ----
	####
	
	
	utrim_bartfit<-bartalone(xtr=test2[,2:ncol(test2)],
													 ytr=test2$y,
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
																	 ytr=test2$y,
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
													ytr=test2$y[which(RO==1)],
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
																	ytr=test2$y[which(RO==1)],
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
													nsamp = 5000,
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
	
	ce_tgr_cor<-gr(Y=test2$y[which(RO==1)],
								 trt=test2$z[which(RO==1)],
								 ps=test2$ps[which(RO==1)],
								 X=test2[which(RO==1),4:ncol(test2)],
								 M=500,qps=quantile(test2$ps[which(RO==1)],probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))
	
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
	
	ugr_fit <-gr(Y=test2$y,
							 trt=test2$z,
							 ps=test2$ps,
							 X=test2[,4:ncol(test2)],
							 M=500,qps=quantile(test2$ps,probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))
	
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
	
	xbcf.model <- XBCF::XBCF(y=test2$y,
														z=test2$z,
														x_con=as.matrix(test2[,4:ncol(test2)]),
														x_mod=as.matrix(test2[,4:ncol(test2)]),
														pihat = test2$ps,
														pcat_con = simdat$p_cat_noise, 
														pcat_mod = simdat$p_cat_noise,
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
	
	print(paste("Finished sim run i =", i, sep=" "))
	
}


####
# Save results ----
####


saveRDS(err_res, 
				paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
							"err_results_", 
							filename_specs,
							"_reps", reps, 
							".rds",
							sep="")
)



saveRDS(cicoverage_res, 
				paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
							"cicoverage_results_", 
							filename_specs,
							"_reps", reps, 
							".rds",
							sep="")
)

saveRDS(ciwidth_res, 
				paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
							"ciwidth_results_", 
							filename_specs,
							"_reps", reps, 
							".rds",
							sep="")
)

	
	

