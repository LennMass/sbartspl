####
# Thyroid cancer and NG compressor station analysis ----
####

# load packages
library(Hmisc)
library(MASS)
library(stats)
library(splines)
library(MCMCpack)
library(BayesTree)
library(dbarts)
library(sf) # alternative for rgdal-package, as this is depreciated. We use sf::read_f() instead of rgdal::readOGR().
library(SoftBart) 
library(tidyverse) 

# Load helper functions
source('01_code/00_functions/01_helper/pw_overlap.R')
source('01_code/00_functions/01_helper/aceBB.R')
source('01_code/00_functions/01_helper/normalize.R')

# Load self-written methods
source('01_code/00_functions/02_methods/bartspl_appl.R')
source('01_code/00_functions/02_methods/bartspl.R')
source('01_code/00_functions/02_methods/bartalone.R')
source('01_code/00_functions/02_methods/softbartalone.R')
source('01_code/00_functions/02_methods/gr.R')






####
# outcome: thyroid cancer mortality rate ----
####

set.seed(4)

## create an empty table to store thyroid cancer analysis results ##
results_tab<-NULL

# read in dataset
test2 <- readRDS("00_data/02_application/02_preprocess/thycancer.rds")

## identify the RO and RN ##
RO<-pw_overlap(ps=test2$ps,E=test2$expose,a=.1*(max(test2$ps)-min(test2$ps)),b=10)

####
# Estimate effect of exposure on thyroid mortality rates ----
####

####
## untrimmed BART ----
####

test3<-test2
test3$expose<-1-test3$expose
utrim_bartfit<-bartalone(xtr=test2[,c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],
												 ytr=test2[,2],xte=test3[,c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))])

ubart <- list(ace_pd = utrim_bartfit$ace_pd,
							ate = utrim_bartfit$ace_pd %>% mean(), 
							atel = utrim_bartfit$ace_pd %>% quantile(0.025),
							ateu = utrim_bartfit$ace_pd %>% quantile(0.975))

####
## untrimmed SoftBART ----
####

test3<-test2
test3$expose<-1-test3$expose
utrim_softbartfit<-softbartalone(xtr=test2[,c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],
																 ytr=test2[,2],
																 xte=test3[,c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))])

usoftbart <- list(ace_pd = utrim_softbartfit$ace_pd,
									ate= utrim_softbartfit$ace_pd %>% mean(),
									atel = utrim_softbartfit$ace_pd %>% quantile(0.025), 
									ateu = utrim_softbartfit$ace_pd %>% quantile(0.975))

####
## trimmed BART ----
####

test3<-test2
test3$expose<-1-test3$expose
trim_bartfit<-bartalone(xtr=test2[which(RO==1),c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],
												ytr=test2[which(RO==1),2],xte=test3[which(RO==1),c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))])

tbart <- list(ace_pd = trim_bartfit$ace_pd,
							ate = trim_bartfit$ace_pd %>% mean(), 
							atel = trim_bartfit$ace_pd %>% quantile(0.025),
							ateu = trim_bartfit$ace_pd %>% quantile(0.975))


####
## trimmed SoftBART ----
####

test3<-test2
test3$expose<-1-test3$expose
trim_softbartfit<-softbartalone(xtr=test2[which(RO==1),c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],
																ytr=test2[which(RO==1),2],xte=test3[which(RO==1),c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))])

tsoftbart <- list(ace_pd = trim_softbartfit$ace_pd,
									ate = trim_softbartfit$ace_pd %>% mean(), 
									atel = trim_softbartfit$ace_pd %>% quantile(0.025),
									ateu = trim_softbartfit$ace_pd %>% quantile(0.975))

####
## BART+SPL ----
####

bartspl_model <-bartspl_appl(datall=test2[,c(2,ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],RO=RO)

bartspl.fit <- list(ace_pd = bartspl_model$ace_pd,
										ate = bartspl_model$ace_pd %>% mean(), 
										atel = bartspl_model$ace_pd %>% quantile(0.025), 
										ateu = bartspl_model$ace_pd %>% quantile(0.975))

####
## SBART+SPL ----
####

sbartspl_model<-bartspl(datall=test2[,c(2,ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],
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

####
## trimmed GR ----
####

ce_tgr_cor<-gr(Y=test2[which(RO==1),2],
							 trt=test2$expose[which(RO==1)],
							 ps=test2$ps[which(RO==1)],
							 X=test2[which(RO==1),3:(ncol(test2)-2)],
							 M=500,qps=quantile(test2$ps[which(RO==1)],probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))

tgr <-list(ace_pd = ce_tgr_cor$ace,
					 ate = ce_tgr_cor$ace %>% mean(), 
					 atel = ce_tgr_cor$ace %>% quantile(0.025),
					 ateu = ce_tgr_cor$ace %>% quantile(0.975))




####
## untrimmed GR ----
####

ugr_fit <-gr(Y=test2[,2],
						 trt=test2$expose,
						 ps=test2$ps,
						 X=test2[,3:(ncol(test2)-2)],
						 M=500,qps=quantile(test2$ps,probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))

ugr <- list(ace_pd = ugr_fit$ace %>% as.vector(), 
						ate= ugr_fit$ace %>%  mean(),
						atel = ugr_fit$ace %>% quantile(.025), 
						ateu = ugr_fit$ace %>% quantile(.975)
)


####
## XBCF ----
####

xbcf.model <- XBCF:: XBCF(y=test2[,2],
													z=test2$expose,
													x_con=as.matrix(test2[,c(3:(ncol(test2)-2))]),
													x_mod=as.matrix(test2[,c(3:(ncol(test2)-2))]),
													pihat = test2$ps,
													pcat_con = 0, 
													pcat_mod = 0,
													n_trees_mod = 20, num_sweeps = 100, Nmin = 20)
xbcf.tau <- xbcf.model$tauhats.adjusted
xbcf.ace_pd <- apply(xbcf.tau, 2, aceBB)
xbcf <- list(ace_pd = xbcf.ace_pd, 
						 ate = xbcf.ace_pd %>% mean(), 
						 atel = xbcf.ace_pd %>% quantile(0.025), 
						 ateu = xbcf.ace_pd %>% quantile(0.975))


####
## GAM, untrimmed ----
####

# standardize covariate variables and outcome
test2_scale <- test2
cont_vars <- colnames(test2)[-c(2, ncol(test2)-1)]
test2_scale[cont_vars] <- scale(test2[cont_vars])
test2_scale$out <- scale(test2$out)


form <- brms::bf(
	paste("out ~ expose",
				paste(sprintf("s(%s, k = 3)", cont_vars), collapse = " + "),
				sep = " + "))

ugam_fit <- fit <- brms::brm(
	form,
	data = test2_scale,
	chains = 4,
	cores = 4, 
	control = list(adapt_delta = 0.995)
)


ugam_m1 <- brms::posterior_epred(ugam_fit, newdata = transform(test2_scale, expose = 1))
ugam_m0 <- brms::posterior_epred(ugam_fit, newdata = transform(test2_scale, expose = 0))

# ATE posterior draws, back-transformed on original scale
ugam_ATE_post <- rowMeans(ugam_m1 - ugam_m0) * attr(test2_scale$out, "scaled:scale")
ugam <- list(ace_pd=ugam_ATE_post,
						 ate = ugam_ATE_post %>% mean(), 
						 atel = ugam_ATE_post %>% quantile(0.025), 
						 ateu = ugam_ATE_post %>% quantile(0.975) 
)


####
## GAM, trimmed ----
####


tgam_fit <- fit <- brms::brm(
	form,
	data = test2_scale[which(RO==1),],
	chains = 4,
	cores = 4, 
	control = list(adapt_delta = 0.995)
)


tgam_m1 <- brms::posterior_epred(tgam_fit, newdata = transform(test2_scale[which(RO==1),], expose = 1))
tgam_m0 <- brms::posterior_epred(tgam_fit, newdata = transform(test2_scale[which(RO==1),], expose = 0))

# ATE posterior draws, back-transformed on original scale
tgam_ATE_post <- rowMeans(tgam_m1 - tgam_m0) * attr(test2_scale$out, "scaled:scale")
tgam <- list(ace_pd=tgam_ATE_post, 
						 ate= tgam_ATE_post %>% mean(), 
						 atel = tgam_ATE_post %>% quantile(0.025),
						 ateu = tgam_ATE_post %>% quantile(0.975))


####
## Save results ----
####

results_tab<-rbind(results_tab,
									 c(ugr$ate, ugr$atel, ugr$ateu, ugr$ateu - ugr$atel),
									 c(ubart$ate, ubart$atel, ubart$ateu, ubart$ateu - ubart$atel),
									 c(usoftbart$ate, usoftbart$atel, usoftbart$ateu, usoftbart$ateu - usoftbart$atel),
									 c(ugam$ate, ugam$atel, ugam$ateu, ugam$ateu - ugam$atel),
									 c(tbart$ate, tbart$atel, tbart$ateu, tbart$ateu - tbart$atel),
									 c(tsoftbart$ate, tsoftbart$atel, tsoftbart$ateu, tsoftbart$ateu - tsoftbart$atel),
									 c(tgr$ate,tgr$atel, tgr$ateu, tgr$ateu - tgr$atel),
									 c(tgam$ate, tgam$atel, tgam$ateu, tgam$ateu - tgam$atel),
									 c(xbcf$ate, xbcf$atel, xbcf$ateu, xbcf$ateu - xbcf$atel),
									 c(bartspl.fit$ate, bartspl.fit$atel, bartspl.fit$ateu, bartspl.fit$ateu - bartspl.fit$atel),
									 c(sbartspl.fit$ate, sbartspl.fit$atel, sbartspl.fit$ateu, sbartspl.fit$ateu - sbartspl.fit$atel)
)
methods <- c("UGR", "UBART", "USBART", "UGAM", "TGR", "TBART", "TSBART", "TGAM", "XBCF", "BARTSPL", "SBARTSPL")
colnames(results_tab) <- c( "mean", "lower", "upper", "width")
results_tab <- data.frame(round(results_tab, digits=3))
results_tab <- cbind(method=methods, results_tab)

saveRDS(results_tab, "02_results/02_application/results_thycancer.rds")





####
# outcome: CHANGE in thyroid cancer mortality rate ----
####

set.seed(4)

## create an empty table to store change in thyroid cancer analysis results ##
results_tab<-NULL

# read in dataset
test2 <- readRDS("00_data/02_application/02_preprocess/thycancer_change.rds")

# Identify the RO and RN 
RO<-pw_overlap(ps=test2$ps,E=test2$expose,a=.1*(max(test2$ps)-min(test2$ps)),b=10)

####
# Estimate effect of exposure on CHANGE in thyroid cancer mortality rates ----
####

####
## untrimmed BART ----
####

test3<-test2
test3$expose<-1-test3$expose
utrim_bartfit<-bartalone(xtr=test2[,c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],
												 ytr=test2[,2],xte=test3[,c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))])

ubart <- list(ace_pd = utrim_bartfit$ace_pd,
							ate = utrim_bartfit$ace_pd %>% mean(), 
							atel = utrim_bartfit$ace_pd %>% quantile(0.025),
							ateu = utrim_bartfit$ace_pd %>% quantile(0.975))

####
## untrimmed SoftBART ----
####

test3<-test2
test3$expose<-1-test3$expose
utrim_softbartfit<-softbartalone(xtr=test2[,c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],
																 ytr=test2[,2],
																 xte=test3[,c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))])

usoftbart <- list(ace_pd = utrim_softbartfit$ace_pd,
									ate= utrim_softbartfit$ace_pd %>% mean(),
									atel = utrim_softbartfit$ace_pd %>% quantile(0.025), 
									ateu = utrim_softbartfit$ace_pd %>% quantile(0.975))

####
## trimmed BART ----
####

test3<-test2
test3$expose<-1-test3$expose
trim_bartfit<-bartalone(xtr=test2[which(RO==1),c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],
												ytr=test2[which(RO==1),2],xte=test3[which(RO==1),c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))])

tbart <- list(ace_pd = trim_bartfit$ace_pd,
							ate = trim_bartfit$ace_pd %>% mean(), 
							atel = trim_bartfit$ace_pd %>% quantile(0.025),
							ateu = trim_bartfit$ace_pd %>% quantile(0.975))


####
## trimmed SoftBART ----
####

test3<-test2
test3$expose<-1-test3$expose
trim_softbartfit<-softbartalone(xtr=test2[which(RO==1),c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],
																ytr=test2[which(RO==1),2],xte=test3[which(RO==1),c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))])

tsoftbart <- list(ace_pd = trim_softbartfit$ace_pd,
									ate = trim_softbartfit$ace_pd %>% mean(), 
									atel = trim_softbartfit$ace_pd %>% quantile(0.025),
									ateu = trim_softbartfit$ace_pd %>% quantile(0.975))

####
## BART+SPL ----
####

bartspl_model <-bartspl_appl(datall=test2[,c(2,ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],RO=RO)

bartspl.fit <- list(ace_pd = bartspl_model$ace_pd,
										ate = bartspl_model$ace_pd %>% mean(), 
										atel = bartspl_model$ace_pd %>% quantile(0.025), 
										ateu = bartspl_model$ace_pd %>% quantile(0.975))

####
## SBART+SPL ----
####

sbartspl_model<-bartspl(datall=test2[,c(2,ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],
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

####
## trimmed GR ----
####

ce_tgr_cor<-gr(Y=test2[which(RO==1),2],
							 trt=test2$expose[which(RO==1)],
							 ps=test2$ps[which(RO==1)],
							 X=test2[which(RO==1),3:(ncol(test2)-2)],
							 M=500,qps=quantile(test2$ps[which(RO==1)],probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))

tgr <-list(ace_pd = ce_tgr_cor$ace,
					 ate = ce_tgr_cor$ace %>% mean(), 
					 atel = ce_tgr_cor$ace %>% quantile(0.025),
					 ateu = ce_tgr_cor$ace %>% quantile(0.975))




####
## untrimmed GR ----
####

ugr_fit <-gr(Y=test2[,2],
						 trt=test2$expose,
						 ps=test2$ps,
						 X=test2[,3:(ncol(test2)-2)],
						 M=500,qps=quantile(test2$ps,probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))

ugr <- list(ace_pd = ugr_fit$ace %>% as.vector(), 
						ate= ugr_fit$ace %>%  mean(),
						atel = ugr_fit$ace %>% quantile(.025), 
						ateu = ugr_fit$ace %>% quantile(.975)
)


####
## XBCF ----
####

xbcf.model <- XBCF:: XBCF(y=test2[,2],
													z=test2$expose,
													x_con=as.matrix(test2[,c(3:(ncol(test2)-2))]),
													x_mod=as.matrix(test2[,c(3:(ncol(test2)-2))]),
													pihat = test2$ps,
													pcat_con = 0, 
													pcat_mod = 0,
													n_trees_mod = 20, num_sweeps = 100, Nmin = 20)
xbcf.tau <- xbcf.model$tauhats.adjusted
xbcf.ace_pd <- apply(xbcf.tau, 2, aceBB)
xbcf <- list(ace_pd = xbcf.ace_pd, 
						 ate = xbcf.ace_pd %>% mean(), 
						 atel = xbcf.ace_pd %>% quantile(0.025), 
						 ateu = xbcf.ace_pd %>% quantile(0.975))


####
## GAM, untrimmed ----
####

# standardize covariate variables and outcome
test2_scale <- test2
cont_vars <- colnames(test2)[-c(2, ncol(test2)-1)]
test2_scale[cont_vars] <- scale(test2[cont_vars])
test2_scale$out <- scale(test2$out)


form <- brms::bf(
	paste("out ~ expose",
				paste(sprintf("s(%s, k = 3)", cont_vars), collapse = " + "),
				sep = " + "))

ugam_fit <- fit <- brms::brm(
	form,
	data = test2_scale,
	chains = 4,
	cores = 4, 
	control = list(adapt_delta = 0.995)
)


ugam_m1 <- brms::posterior_epred(ugam_fit, newdata = transform(test2_scale, expose = 1))
ugam_m0 <- brms::posterior_epred(ugam_fit, newdata = transform(test2_scale, expose = 0))

# ATE posterior draws, back-transformed on original scale
ugam_ATE_post <- rowMeans(ugam_m1 - ugam_m0) * attr(test2_scale$out, "scaled:scale")
ugam <- list(ace_pd=ugam_ATE_post,
						 ate = ugam_ATE_post %>% mean(), 
						 atel = ugam_ATE_post %>% quantile(0.025), 
						 ateu = ugam_ATE_post %>% quantile(0.975) 
)


####
## GAM, trimmed ----
####


tgam_fit <- fit <- brms::brm(
	form,
	data = test2_scale[which(RO==1),],
	chains = 4,
	cores = 4, 
	control = list(adapt_delta = 0.995)
)


tgam_m1 <- brms::posterior_epred(tgam_fit, newdata = transform(test2_scale[which(RO==1),], expose = 1))
tgam_m0 <- brms::posterior_epred(tgam_fit, newdata = transform(test2_scale[which(RO==1),], expose = 0))

# ATE posterior draws, back-transformed on original scale
tgam_ATE_post <- rowMeans(tgam_m1 - tgam_m0) * attr(test2_scale$out, "scaled:scale")
tgam <- list(ace_pd=tgam_ATE_post, 
						 ate= tgam_ATE_post %>% mean(), 
						 atel = tgam_ATE_post %>% quantile(0.025),
						 ateu = tgam_ATE_post %>% quantile(0.975))


####
## Save results ----
####

results_tab<-rbind(results_tab,
									 c(ugr$ate, ugr$atel, ugr$ateu, ugr$ateu - ugr$atel),
									 c(ubart$ate, ubart$atel, ubart$ateu, ubart$ateu - ubart$atel),
									 c(usoftbart$ate, usoftbart$atel, usoftbart$ateu, usoftbart$ateu - usoftbart$atel),
									 c(ugam$ate, ugam$atel, ugam$ateu, ugam$ateu - ugam$atel),
									 c(tbart$ate, tbart$atel, tbart$ateu, tbart$ateu - tbart$atel),
									 c(tsoftbart$ate, tsoftbart$atel, tsoftbart$ateu, tsoftbart$ateu - tsoftbart$atel),
									 c(tgr$ate,tgr$atel, tgr$ateu, tgr$ateu - tgr$atel),
									 c(tgam$ate, tgam$atel, tgam$ateu, tgam$ateu - tgam$atel),
									 c(xbcf$ate, xbcf$atel, xbcf$ateu, xbcf$ateu - xbcf$atel),
									 c(bartspl.fit$ate, bartspl.fit$atel, bartspl.fit$ateu, bartspl.fit$ateu - bartspl.fit$atel),
									 c(sbartspl.fit$ate, sbartspl.fit$atel, sbartspl.fit$ateu, sbartspl.fit$ateu - sbartspl.fit$atel)
)
methods <- c("UGR", "UBART", "USBART", "UGAM", "TGR", "TBART", "TSBART", "TGAM", "XBCF", "BARTSPL", "SBARTSPL")
colnames(results_tab) <- c( "mean", "lower", "upper", "width")
results_tab <- data.frame(round(results_tab, digits=3))
results_tab <- cbind(method=methods, results_tab)

saveRDS(results_tab, "02_results/02_application/results_thycancer_change.rds")

