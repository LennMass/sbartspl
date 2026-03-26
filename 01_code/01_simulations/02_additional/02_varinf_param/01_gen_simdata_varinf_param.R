####
# DGP based on BCF
####

simnum <- 5 # seed
set.seed(simnum)

# parameters
n = 500 # sample size
tau_setting="homogeneous" # c("homogeneous", "heterogeneous") 
p_cont_noise = 10
p_cat_noise = 10 
reps=20

filename_specs <- paste("n", n, 
												"_tausetting", tau_setting, 
												"_pnoise", p_cont_noise+p_cat_noise,
												"_seed", simnum,
												sep="")

for(i in 1:reps){
	
	create_dataset <- dgp_bcf(n=n, 
														tau_setting=tau_setting,
														p_cat_noise=p_cat_noise, 
														p_cont_noise=p_cont_noise)
	saveRDS(create_dataset, paste("00_data/01_simulations/02_additional/02_varinf_param/",
																filename_specs, 
																"_simid", i, 
																".rds",
																sep="")
	)
}


####
# DGP based on BARTSPL
####

simnum <- 5 # seed
set.seed(simnum)

# parameters
N = 500 # sample size
ncovs = 10
reps=20

filename_specs <- paste("n", N, 
												"_ncovs", ncovs, 
												"_seed", simnum,
												sep="")


for(i in 1:reps){
	
	create_dataset <- dgp_bartspl(N=N, 
																ncovs = ncovs)
	saveRDS(create_dataset, paste("00_data/01_simulations/02_additional/02_varinf_param/",
																filename_specs, 
																"_simid", i, 
																".rds",
																sep="")
	)
	
}



