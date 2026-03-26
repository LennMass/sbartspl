library(usmap)
library(tidyverse)
library(cowplot)



##### Mortality rates ####

###### Leukemia ####

dat <- readRDS("00_data/02_application/02_preprocess/leukemia.rds") %>%
	dplyr::select(fips, out, ps, expose) %>%
	dplyr::transmute(fips=fips, 
									 values=exp(out), # to get deaths per 100 000 population (rescale from log),
									 ps=ps, 
									 expose=expose
									 )

# mid western focus
p_leukemia <- usmap::plot_usmap(data=dat, 
									regions = "counties", 
									include = c("AZ", "AR", "CO", "IL", "IA", "KS", "LA", "MN", "MS", "MO", 
															"MT", "NE", "NM", "ND", "OK", "SD", "TX", "UT", "WI", "WY")) +
	scale_fill_continuous(name = "Leukemia mortality (2014)", na.value = "transparent")+ 
	theme(legend.position = "right")


###### Leukemia change ####

dat <- readRDS("00_data/02_application/02_preprocess/leukemia_change.rds") %>%
	dplyr::select(fips, out) %>%
	dplyr::transmute(fips=fips, 
									 values=out # percentage points
	)

# mid western focus
p_leukemiachange <- usmap::plot_usmap(data=dat, 
									regions = "counties", 
									include = c("AZ", "AR", "CO", "IL", "IA", "KS", "LA", "MN", "MS", "MO", 
															"MT", "NE", "NM", "ND", "OK", "SD", "TX", "UT", "WI", "WY")) +
	scale_fill_continuous(name = "Leukemia change (1980 to 2014)", na.value = "transparent")+ 
	theme(legend.position = "right")


###### Thyroid ####

dat <- readRDS("00_data/02_application/02_preprocess/thycancer.rds") %>%
	dplyr::select(fips, out) %>%
	dplyr::transmute(fips=fips, 
									 values=exp(out) # percentage points
	)

# mid western focus
p_thy <-  usmap::plot_usmap(data=dat, 
									regions = "counties", 
									include = c("AZ", "AR", "CO", "IL", "IA", "KS", "LA", "MN", "MS", "MO", 
															"MT", "NE", "NM", "ND", "OK", "SD", "TX", "UT", "WI", "WY")) +
	scale_fill_continuous(name = "Thyorid mortality (2014)", na.value = "transparent")+ 
	theme(legend.position = "right")



###### Thyroid change ####

dat <- readRDS("00_data/02_application/02_preprocess/thycancer_change.rds") %>%
	dplyr::select(fips, out) %>%
	dplyr::transmute(fips=fips, 
									 values=out # percentage points
	)

# mid western focus
p_thychange <- usmap::plot_usmap(data=dat, 
									regions = "counties", 
									include = c("AZ", "AR", "CO", "IL", "IA", "KS", "LA", "MN", "MS", "MO", 
															"MT", "NE", "NM", "ND", "OK", "SD", "TX", "UT", "WI", "WY")) +
	scale_fill_continuous(name = "Thyroid change (1980 to 2014)",na.value = "transparent") + 
	theme(legend.position = "right")



mort_2014 <- plot_grid(p_leukemia, p_thy, nrow=2)

change_19802014 <- plot_grid(p_leukemiachange, p_thychange, nrow=2)




#### Region of overlap/non-overlap ####

RO<-pw_overlap(ps=dat$ps,E=dat$expose,a=.1*(max(dat$ps)-min(dat$ps)),b=10)

dat <- data.frame(fips=dat$fips, values=as.factor(RO)) 

overlap_plot <- usmap::plot_usmap(data=dat, 
									regions = "counties", 
									include = c("AZ", "AR", "CO", "IL", "IA", "KS", "LA", "MN", "MS", "MO", 
															"MT", "NE", "NM", "ND", "OK", "SD", "TX", "UT", "WI", "WY"), ) +
	scale_fill_manual(name = "Region",
										values = c("0" = "orange", "1" = "#2166ac"),
										na.value = "transparent",
										labels = c("0" = "Non-overlap", "1" = "Overlap")) +
	theme(legend.position = "right")

