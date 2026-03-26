####
# Analysis of variance inflation parameter 
# Based on DGP BARTSPL
####


####
# Set simulation params ----
####

# parameters
N = 500 # sample size
ncovs = 10
reps=20
simnum <- 5
set.seed(simnum)

filename_specs <- paste("n", N, 
												"_ncovs", ncovs, 
												"_seed", simnum,
												sep="")


# save metrics
method_names <- c("BARTSPL", "BARTSPL_vp0", "SBARTSPL")
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
	
	simdat <- readRDS(paste("00_data/01_simulations/02_additional/02_varinf_param/",
													filename_specs, 
													"_simid", i, 
													".rds",
													sep="")
	)
	
	
	test2<- simdat$untrimmed_dat
	RO <- simdat$RO
	
	
	####
	# Start estimation ----
	####
	
	
	####
	## BART+SPL, var_infl_param = 10 (recommended default in Nethery et al. (2019)) ----
	####
	
	bartspl_model <-bartspl(datall=test2,RO=RO)
	
	bartspl.fit <- list(ace_pd = bartspl_model$ace_pd,
											ate = bartspl_model$ace_pd %>% mean(), 
											atel = bartspl_model$ace_pd %>% quantile(0.025), 
											ateu = bartspl_model$ace_pd %>% quantile(0.975))
	
	err_res[i, "BARTSPL"] <- bartspl.fit$ate - simdat$ate_true
	cicoverage_res[i, "BARTSPL"] <- ifelse( (bartspl.fit$atel <= simdat$ate_true) && (bartspl.fit$ateu >= simdat$ate_true), 1, 0)
	ciwidth_res[i, "BARTSPL"] <- bartspl.fit$ateu - bartspl.fit$atel
	
	print("Finished BART+SPL with var_infl_param = 10")
	
	####
	## BART+SPL, var_infl_param = 0, no variance inflation in RN ----
	####
	
	bartspl_model <-bartspl_vp(datall=test2,RO=RO, var_infl_param=0)
	
	bartspl.fit.vp0 <- list(ace_pd = bartspl_model$ace_pd,
													ate = bartspl_model$ace_pd %>% mean(), 
													atel = bartspl_model$ace_pd %>% quantile(0.025), 
													ateu = bartspl_model$ace_pd %>% quantile(0.975))
	
	err_res[i, "BARTSPL_vp0"] <- bartspl.fit.vp0$ate - simdat$ate_true
	cicoverage_res[i, "BARTSPL_vp0"] <- ifelse( (bartspl.fit.vp0$atel <= simdat$ate_true) && (bartspl.fit.vp0$ateu >= simdat$ate_true), 1, 0)
	ciwidth_res[i, "BARTSPL_vp0"] <- bartspl.fit.vp0$ateu - bartspl.fit.vp0$atel
	
	print("Finished BART+SPL with var_infl_param = 0")
	
	
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
	
	
	
	print(paste("Finished sim run i =", i, sep=" "))
	
}


####
# Save results ----
####


saveRDS(err_res, 
				paste("02_results/01_simulations/02_additional/02_varinf_param/",
							"err_results_", 
							filename_specs,
							"_reps", reps, 
							".rds",
							sep="")
)



saveRDS(cicoverage_res, 
				paste("02_results/01_simulations/02_additional/02_varinf_param/",
							"cicoverage_results_", 
							filename_specs,
							"_reps", reps, 
							".rds",
							sep="")
)

saveRDS(ciwidth_res, 
				paste("02_results/01_simulations/02_additional/02_varinf_param/",
							"ciwidth_results_", 
							filename_specs,
							"_reps", reps, 
							".rds",
							sep="")
)