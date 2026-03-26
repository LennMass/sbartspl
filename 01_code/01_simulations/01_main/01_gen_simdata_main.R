####
# DGP based on BARTSPL
####


####
# Generate 100 datasets with 10 different seeds
####

simnum_vec <- seq(1237, 1248, 1) # seed

# parameters
N = 2000 # c(500, 2000)
ncovs = 25 # c(10, 25, 50)
reps=10



for(j in 1:length(simnum_vec)){
	
	simnum <- simnum_vec[j]
	set.seed(simnum) # 10 seeds, 10 datasets per seed, gives 100 simulation datasets 
	
	filename_specs <- paste("n", N, 
													"_ncovs", ncovs, 
													"_seed", simnum,
													sep="")
	
	for(i in 1:reps){
		
		create_dataset <- dgp_bartspl(N=N, 
																	ncovs = ncovs)
		saveRDS(create_dataset, paste("00_data/01_simulations/01_main/",
																	filename_specs, 
																	"_simid", i, 
																	".rds",
																	sep="")
		)
		
	}
}






