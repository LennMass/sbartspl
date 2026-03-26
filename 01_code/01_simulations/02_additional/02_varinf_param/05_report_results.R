
####
# Table 5 results ATE ----
####


####
## First panel ----
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

err <- readRDS(paste("02_results/01_simulations/02_additional/02_varinf_param/",
							"err_results_", 
							filename_specs,
							"_reps", reps, 
							".rds",
							sep=""))



cicoverage <- readRDS(paste("02_results/01_simulations/02_additional/02_varinf_param/",
							"cicoverage_results_", 
							filename_specs,
							"_reps", reps, 
							".rds",
							sep=""))

ciwidth <-  readRDS(paste("02_results/01_simulations/02_additional/02_varinf_param/",
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
	
xtable(ate_table,digits = 4)


####
## Second panel ----
####

# parameters
n = 500 
tau_setting="homogeneous" # c("homogeneous", "heterogeneous")
p_cont_noise = 10 
p_cat_noise = 10 
reps <- 20
simnum <- 5
set.seed(simnum)


filename_specs <- paste("n", n, 
												"_tausetting", tau_setting, 
												"_pnoise", p_cont_noise+p_cat_noise,
												"_seed", simnum,
												sep="")


err <- readRDS(paste("02_results/01_simulations/02_additional/02_varinf_param/",
										 "err_results_", 
										 filename_specs,
										 "_reps", reps, 
										 ".rds",
										 sep=""))



cicoverage <- readRDS(paste("02_results/01_simulations/02_additional/02_varinf_param/",
														"cicoverage_results_", 
														filename_specs,
														"_reps", reps, 
														".rds",
														sep=""))

ciwidth <-  readRDS(paste("02_results/01_simulations/02_additional/02_varinf_param/",
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

xtable(ate_table,digits = 4)


####
## Third panel ----
####

# parameters
n = 500 
tau_setting="heterogeneous" # c("homogeneous", "heterogeneous")
p_cont_noise = 10 
p_cat_noise = 10 
reps <- 20
simnum <- 5
set.seed(simnum)


filename_specs <- paste("n", n, 
												"_tausetting", tau_setting, 
												"_pnoise", p_cont_noise+p_cat_noise,
												"_seed", simnum,
												sep="")


err <- readRDS(paste("02_results/01_simulations/02_additional/02_varinf_param/",
										 "err_results_", 
										 filename_specs,
										 "_reps", reps, 
										 ".rds",
										 sep=""))



cicoverage <- readRDS(paste("02_results/01_simulations/02_additional/02_varinf_param/",
														"cicoverage_results_", 
														filename_specs,
														"_reps", reps, 
														".rds",
														sep=""))

ciwidth <-  readRDS(paste("02_results/01_simulations/02_additional/02_varinf_param/",
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

xtable(ate_table,digits = 4)





