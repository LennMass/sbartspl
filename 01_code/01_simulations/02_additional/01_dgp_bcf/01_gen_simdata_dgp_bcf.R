simnum <- 123 # seed
set.seed(simnum)

# parameters
n = 500 # sample size
tau_setting="homogeneous" # c("homogeneous", "heterogeneous") 
p_cont_noise = 5 # c(5, 10)
p_cat_noise = 5 # c(5, 10)
reps=20

filename_specs <- paste("n", n, 
									"_tausetting", tau_setting, 
									"_pnoise", p_cont_noise+p_cat_noise,
									"_seed", simnum, 
									sep="")

for(i in 1:reps){
	
	create_dataset <- dgp_bcf(n=n, 
														p_cont_noise = p_cont_noise,
														p_cat_noise = p_cat_noise, 
														tau_setting = tau_setting)
	saveRDS(create_dataset, paste("00_data/01_simulations/02_additional/01_dgp_bcf/",
												filename_specs, 
												"_simid", i, 
												".rds",
												sep="")
												)
}
