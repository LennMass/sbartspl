

####
# Load packages ---- 
####

# used in application
# load packages
library(here)
library(Hmisc)
library(MASS)
library(stats)
library(splines)
library(MCMCpack)
library(BayesTree)
library(dbarts)
library(sf) # alternative for rgdal-package, as this is depreciated. We use sf::read_f() instead of rgdal::readOGR().
library(xtable)
library(SoftBart) 
library(tidyverse)
library(brms)
library(XBCF)
library(invgamma)


# source self-written functions
lapply(
	c(list.files(here::here("01_code/00_functions/01_helper"), pattern = "\\.R$", full.names = TRUE),
		list.files(here::here("01_code/00_functions/02_methods/"), pattern = "\\.R$", full.names = TRUE)
		),
	source
)


