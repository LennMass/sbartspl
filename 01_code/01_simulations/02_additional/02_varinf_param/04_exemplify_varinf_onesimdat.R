####
# Analysis of variance inflation parameter 
# Based on DGP BARTSPL
# Analysis of one simulation dataset
####


####
# Set simulation params ----
####

# parameters
N = 500 # sample size
ncovs = 10
sim_id=14
i=length(sim_id)
simnum <- 5
set.seed(simnum)

filename_specs <- paste("n", N, 
												"_ncovs", ncovs, 
												"_seed", simnum,
												sep="")

# save metrics
method_names <- c("BARTSPL", "BARTSPL_vp0", "SBARTSPL")
err_res <- matrix(NA, nrow = 1,  ncol=length(method_names))
colnames(err_res) <- method_names

cicoverage_res <- matrix(NA, nrow = 1,  ncol=length(method_names))
colnames(cicoverage_res) <- method_names

ciwidth_res <- matrix(NA, nrow = 1,  ncol=length(method_names))
colnames(ciwidth_res) <- method_names

####
# Load and prepare dataset ----
####

simdat <- readRDS(paste("00_data/01_simulations/02_additional/02_varinf_param/",
												filename_specs, 
												"_simid", sim_id, 
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

bartspl_model_vp0 <-bartspl_vp(datall=test2,RO=RO, var_infl_param=0)

bartspl.fit.vp0 <- list(ace_pd = bartspl_model_vp0$ace_pd,
												ate = bartspl_model_vp0$ace_pd %>% mean(), 
												atel = bartspl_model_vp0$ace_pd %>% quantile(0.025), 
												ateu = bartspl_model_vp0$ace_pd %>% quantile(0.975))

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



print(paste("Finished SBART+SPL"))




####
# Region of non-overlap
###

n_RO <- which(RO==0) %>% length()


ites_BARTSPL <- data.frame(id=seq(1:n_RO),
													 method = "BARTSPL",
													 true=simdat$ite_true[which(RO==0)],
													 ITE=bartspl_model$ice_mean[which(RO==0)],
													 lower=bartspl_model$ice_lw[which(RO==0)],
													 upper=bartspl_model$ice_up[which(RO==0)])

ites_BARTSPL_vp0 <- data.frame(id=seq(1:n_RO),
															 method = "BARTSPL_vp0",
															 true=simdat$ite_true[which(RO==0)],
															 ITE=bartspl_model_vp0$ice_mean[which(RO==0)],
															 lower=bartspl_model_vp0$ice_lw[which(RO==0)],
															 upper=bartspl_model_vp0$ice_up[which(RO==0)])


ites_SBARTSPL <- data.frame(id=seq(1:n_RO),
												 method = "SBARTSPL", 
												 true=simdat$ite_true[which(RO==0)],
												 ITE=sbartspl_model$ice_mean[which(RO==0)],
												 lower=sbartspl_model$ice_lw[which(RO==0)],
												 upper=sbartspl_model$ice_up[which(RO==0)])

ites <- rbind(ites_BARTSPL, ites_BARTSPL_vp0, ites_SBARTSPL) 



ite_nonoverlap <- ggplot(ites, aes(x=id, y=ITE, colour=method))+
	geom_point()+
	geom_errorbar(aes(ymin = lower, 
													 ymax = upper), 
											 width=0.2)+
	geom_hline(yintercept = simdat$ate_true, colour="black", alpha = 0.5, linetype="dashed")+
	geom_point(aes(y=true), colour="black", alpha = 0.5)+
	scale_fill_brewer(palette="Set2") + 
	scale_color_brewer(palette="Set2") +
	scale_shape_manual(values=c(1, 4))+
	facet_wrap(~method, nrow=3)+
	theme_bw(base_size=12)

ggsave(paste("ite_nonoverlap_", 
						 filename_specs, 
						 "_simid", sim_id,
						 ".pdf", 
						 sep=""),
			 path = "02_results/01_simulations/02_additional/02_varinf_param/", 
			 plot=ite_nonoverlap, width=6, height=4, units="in")


RMSE <- data.frame(BARTSPL = rmse(ites_BARTSPL$ITE, ites_BARTSPL$true) %>% round(digit=3), 
									 BARTSPL_vp0 = rmse(ites_BARTSPL_vp0$ITE, ites_BARTSPL_vp0$true) %>% round(digit=3),
									 SBARTSPL = rmse(ites_SBARTSPL$ITE, ites_SBARTSPL$true) %>% round(digit=3))


AbsBias <- data.frame(BARTSPL=mean(abs(ites_BARTSPL$ITE - ites_BARTSPL$true)) %>% round(digit=3),
											BARTSPL_vp0 = mean(abs(ites_BARTSPL_vp0$ITE - ites_BARTSPL_vp0$true)) %>% round(digit=3),
											SBARTSP=mean(abs(ites_SBARTSPL$ITE -  ites_SBARTSPL$true)) %>% round(digit=3))


coverage <- data.frame(BARTSPL=mean(ifelse(ites_BARTSPL$lower <= ites_BARTSPL$true & ites_BARTSPL$upper >= ites_BARTSPL$true, 1, 0)) %>% round(2), 
											 BARTSPL_vp0 = mean(ifelse(ites_BARTSPL_vp0$lower <= ites_BARTSPL_vp0$true & ites_BARTSPL_vp0$upper >= ites_BARTSPL_vp0$true, 1, 0)) %>% round(2),
											 SBARTSPL = mean(ifelse(ites_SBARTSPL$lower <= ites_SBARTSPL$true & ites_SBARTSPL$upper >= ites_SBARTSPL$true, 1, 0)) %>% round(2))

width <- data.frame(BARTSPL=mean(ites_BARTSPL$upper - ites_BARTSPL$lower),
										BARTSPL_vp0 = mean(ites_BARTSPL_vp0$upper - ites_BARTSPL_vp0$lower),
										SBARTSPL=mean(ites_SBARTSPL$upper - ites_SBARTSPL$lower))

RN_ite <- data.frame(RMSE = t(RMSE),
										 AbsBias = t(AbsBias), 
										 Cov = t(coverage), 
										 Width = t(width))

saveRDS(RN_ite, paste("02_results/01_simulations/02_additional/02_varinf_param/", 
											"RN_ite_results", 
											".rds", 
											sep=""))

####
# Region of overlap
###

n_RO <- which(RO==1) %>% length()


ites_BARTSPL <- data.frame(id=seq(1:n_RO),
													 method = "BARTSPL",
													 true=simdat$ite_true[which(RO==1)],
													 ITE=bartspl_model$ice_mean[which(RO==1)],
													 lower=bartspl_model$ice_lw[which(RO==1)],
													 upper=bartspl_model$ice_up[which(RO==1)])

ites_BARTSPL_vp0 <- data.frame(id=seq(1:n_RO),
															 method = "BARTSPL_vp0",
															 true=simdat$ite_true[which(RO==1)],
															 ITE=bartspl_model_vp0$ice_mean[which(RO==1)],
															 lower=bartspl_model_vp0$ice_lw[which(RO==1)],
															 upper=bartspl_model_vp0$ice_up[which(RO==1)])


ites_SBARTSPL <- data.frame(id=seq(1:n_RO),
														method = "SBARTSPL", 
														true=simdat$ite_true[which(RO==1)],
														ITE=sbartspl_model$ice_mean[which(RO==1)],
														lower=sbartspl_model$ice_lw[which(RO==1)],
														upper=sbartspl_model$ice_up[which(RO==1)])

ites <- rbind(ites_BARTSPL, ites_BARTSPL_vp0, ites_SBARTSPL) 



ite_overlap <- ggplot(ites, aes(x=id, y=ITE, colour=method))+
	geom_point()+
	geom_errorbar(aes(ymin = lower, 
										ymax = upper), 
								width=0.2)+
	geom_hline(yintercept = simdat$ate_true, colour="black", alpha = 1, linetype="dashed")+
	geom_point(aes(y=true), colour="black", alpha = 0.5)+
	scale_fill_brewer(palette="Set2") + 
	scale_color_brewer(palette="Set2") +
	scale_shape_manual(values=c(1, 4))+
	facet_wrap(~method, nrow=3)+
	theme_bw(base_size=12)


ggsave(paste("ite_overlap_", 
						 filename_specs, 
						 "_simid", sim_id,
						 ".pdf", 
						 sep=""),
			 path = "02_results/01_simulations/02_additional/02_varinf_param/", 
			 plot=ite_overlap, width=6, height=4, units="in")


RMSE <- data.frame(BARTSPL = rmse(ites_BARTSPL$ITE, ites_BARTSPL$true) %>% round(digit=3), 
									 BARTSPL_vp0 = rmse(ites_BARTSPL_vp0$ITE, ites_BARTSPL_vp0$true) %>% round(digit=3),
									 SBARTSPL = rmse(ites_SBARTSPL$ITE, ites_SBARTSPL$true) %>% round(digit=3))


AbsBias <- data.frame(BARTSPL=mean(abs(ites_BARTSPL$ITE - ites_BARTSPL$true)) %>% round(digit=3),
											BARTSPL_vp0 = mean(abs(ites_BARTSPL_vp0$ITE - ites_BARTSPL_vp0$true)) %>% round(digit=3),
											SBARTSP=mean(abs(ites_SBARTSPL$ITE -  ites_SBARTSPL$true)) %>% round(digit=3))
	
	
coverage <- data.frame(BARTSPL=mean(ifelse(ites_BARTSPL$lower <= ites_BARTSPL$true & ites_BARTSPL$upper >= ites_BARTSPL$true, 1, 0)) %>% round(2), 
											 BARTSPL_vp0 = mean(ifelse(ites_BARTSPL_vp0$lower <= ites_BARTSPL_vp0$true & ites_BARTSPL_vp0$upper >= ites_BARTSPL_vp0$true, 1, 0)) %>% round(2),
											 SBARTSPL = mean(ifelse(ites_SBARTSPL$lower <= ites_SBARTSPL$true & ites_SBARTSPL$upper >= ites_SBARTSPL$true, 1, 0)) %>% round(2))

width <- data.frame(BARTSPL=mean(ites_BARTSPL$upper - ites_BARTSPL$lower),
										BARTSPL_vp0 = mean(ites_BARTSPL_vp0$upper - ites_BARTSPL_vp0$lower),
										SBARTSPL=mean(ites_SBARTSPL$upper - ites_SBARTSPL$lower))

RO_ite <- data.frame(RMSE = t(RMSE),
										 AbsBias = t(AbsBias), 
										 Cov = t(coverage), 
										 Width = t(width))

saveRDS(RO_ite, paste("02_results/01_simulations/02_additional/02_varinf_param/", 
											"RO_ite_results", 
											".rds", 
											sep=""))









