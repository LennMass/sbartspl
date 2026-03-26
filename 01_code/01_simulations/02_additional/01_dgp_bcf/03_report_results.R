####
# Table 4 results ATE ----
####


####
## First panel, left ----
####

n = 500 
tau_setting="homogeneous" # c("homogeneous", "heterogeneous")
p_cont_noise = 5 # c(10, 20)
p_cat_noise = 5 # c(10, 20)
reps <- 20
simnum <- 123


filename_specs <- paste("n", n, 
												"_tausetting", tau_setting, 
												"_pnoise", p_cont_noise+p_cat_noise,
												"_seed", simnum,
												sep="")


err <- readRDS(paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
										 "err_results_", 
										 filename_specs,
										 "_reps", reps, 
										 ".rds",
										 sep=""))



cicoverage <- readRDS(paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
														"cicoverage_results_", 
														filename_specs,
														"_reps", reps, 
														".rds",
														sep=""))

ciwidth <-  readRDS(paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
													"ciwidth_results_", 
													filename_specs,
													"_reps", reps, 
													".rds",
													sep=""))


RMSE <- err %>% 
	apply(2, function(x) sqrt( mean( ( x )^2 ) )) %>%
	t() %>%
	t()

ABSBIAS <- err %>% 
	apply(2, function(x) mean( abs( x ) )) %>%
	t() %>%
	t()

COVERAGE <- cicoverage %>%
	colMeans() %>% 
	t() %>%
	t()

WIDTH <- ciwidth %>%
	colMeans() %>% 
	t() %>%
	t()


ate_table <- data.frame(RMSE = RMSE, AbsBias = ABSBIAS, CICoverage = COVERAGE, CIWidth = WIDTH) 

xtable(ate_table,digits = 3)

####
## First panel, right ----
####

n = 500 
tau_setting="heterogeneous" # c("homogeneous", "heterogeneous")
p_cont_noise = 5 # c(10, 20)
p_cat_noise = 5 # c(10, 20)
reps <- 20
simnum <- 123


filename_specs <- paste("n", n, 
												"_tausetting", tau_setting, 
												"_pnoise", p_cont_noise+p_cat_noise,
												"_seed", simnum,
												sep="")


err <- readRDS(paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
										 "err_results_", 
										 filename_specs,
										 "_reps", reps, 
										 ".rds",
										 sep=""))



cicoverage <- readRDS(paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
														"cicoverage_results_", 
														filename_specs,
														"_reps", reps, 
														".rds",
														sep=""))

ciwidth <-  readRDS(paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
													"ciwidth_results_", 
													filename_specs,
													"_reps", reps, 
													".rds",
													sep=""))


RMSE <- err %>% 
	apply(2, function(x) sqrt( mean( ( x )^2 ) )) %>%
	t() %>%
	t()

ABSBIAS <- err %>% 
	apply(2, function(x) mean( abs( x ) )) %>%
	t() %>%
	t()

COVERAGE <- cicoverage %>%
	colMeans() %>% 
	t() %>%
	t()

WIDTH <- ciwidth %>%
	colMeans() %>% 
	t() %>%
	t()


ate_table <- data.frame(RMSE = RMSE, AbsBias = ABSBIAS, CICoverage = COVERAGE, CIWidth = WIDTH) 

xtable(ate_table,digits = 3)


####
## Second panel, left ----
####

# parameters
n = 500 
tau_setting="homogeneous" # c("homogeneous", "heterogeneous")
p_cont_noise = 10 
p_cat_noise = 10 
reps <- 20
simnum <- 123


filename_specs <- paste("n", n, 
												"_tausetting", tau_setting, 
												"_pnoise", p_cont_noise+p_cat_noise,
												"_seed", simnum,
												sep="")


err <- readRDS(paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
										 "err_results_", 
										 filename_specs,
										 "_reps", reps, 
										 ".rds",
										 sep=""))



cicoverage <- readRDS(paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
														"cicoverage_results_", 
														filename_specs,
														"_reps", reps, 
														".rds",
														sep=""))

ciwidth <-  readRDS(paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
													"ciwidth_results_", 
													filename_specs,
													"_reps", reps, 
													".rds",
													sep=""))


RMSE <- err %>% 
	apply(2, function(x) sqrt( mean( ( x )^2 , na.rm=T) )) %>%
	t() %>%
	t()

ABSBIAS <- err %>% 
	apply(2, function(x) mean( abs( x ) , na.rm=T)) %>%
	t() %>%
	t()

COVERAGE <- cicoverage %>%
	colMeans(na.rm=T) %>% 
	t() %>%
	t()

WIDTH <- ciwidth %>%
	colMeans(na.rm=T) %>% 
	t() %>%
	t()


ate_table <- data.frame(RMSE = RMSE, AbsBias = ABSBIAS, CICoverage = COVERAGE, CIWidth = WIDTH) 

xtable(ate_table,digits = 3)


####
## Second panel, right ----
####

# parameters
n = 500 
tau_setting="heterogeneous" # c("homogeneous", "heterogeneous")
p_cont_noise = 10 
p_cat_noise = 10 
reps <- 20
simnum <- 123


filename_specs <- paste("n", n, 
												"_tausetting", tau_setting, 
												"_pnoise", p_cont_noise+p_cat_noise,
												"_seed", simnum,
												sep="")


err <- readRDS(paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
										 "err_results_", 
										 filename_specs,
										 "_reps", reps, 
										 ".rds",
										 sep=""))



cicoverage <- readRDS(paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
														"cicoverage_results_", 
														filename_specs,
														"_reps", reps, 
														".rds",
														sep=""))

ciwidth <-  readRDS(paste("02_results/01_simulations/02_additional/01_dgp_bcf/",
													"ciwidth_results_", 
													filename_specs,
													"_reps", reps, 
													".rds",
													sep=""))


RMSE <- err %>% 
	apply(2, function(x) sqrt( mean( ( x )^2 , na.rm=T) )) %>%
	t() %>%
	t()

ABSBIAS <- err %>% 
	apply(2, function(x) mean( abs( x ) , na.rm=T)) %>%
	t() %>%
	t()

COVERAGE <- cicoverage %>%
	colMeans(na.rm=T) %>% 
	t() %>%
	t()

WIDTH <- ciwidth %>%
	colMeans(na.rm=T) %>% 
	t() %>%
	t()


ate_table <- data.frame(RMSE = RMSE, AbsBias = ABSBIAS, CICoverage = COVERAGE, CIWidth = WIDTH) 

xtable(ate_table,digits = 3)







