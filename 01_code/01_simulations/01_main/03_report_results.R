####
# Table 1 results ATE ----
####



# parameters
N = 2000 # sample size
ncovs = 50
reps=10
simnum_vec <- seq(1237, 1248, 1) # seed

err <- NULL
cicoverage <- NULL
ciwidth <- NULL

for(j in 1:length(simnum_vec)){
	
	simnum <- simnum_vec[j]

	filename_specs <- paste("n", N, 
													"_ncovs", ncovs, 
													"_seed", simnum,
													sep="")
	
	err_current <- readRDS(paste("02_results/01_simulations/01_main/",
											 "err_results_", 
											 filename_specs,
											 "_reps", reps, 
											 ".rds",
											 sep=""))
	err <- rbind(err, err_current)
	
	
	
	cicoverage_current <- readRDS(paste("02_results/01_simulations/01_main/",
															"cicoverage_results_", 
															filename_specs,
															"_reps", reps, 
															".rds",
															sep=""))
	cicoverage <- rbind(cicoverage, cicoverage_current)
	
	
	ciwidth_current <-  readRDS(paste("02_results/01_simulations/01_main/",
														"ciwidth_results_", 
														filename_specs,
														"_reps", reps, 
														".rds",
														sep=""))
	ciwidth <- rbind(ciwidth, ciwidth_current)
	
}


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


