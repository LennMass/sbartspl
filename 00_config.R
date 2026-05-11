####
# Setup: Load (and install if needed) packages ----
####

# {pacman} bootstraps everything else: p_load() installs from CRAN if missing,
# then attaches; p_load_gh() does the same for GitHub-only packages.
if (!requireNamespace("pacman", quietly = TRUE)) {
	install.packages("pacman", repos = "https://cloud.r-project.org")
}

# ---- CRAN packages  ----
# Note: {sf} replaces the deprecated {rgdal} in Nethery et al. (2019). We use sf::st_read() instead of
# rgdal::readOGR(). {stats} and {splines} ship with base R, no install needed.
pacman::p_load(
	here,        # project-relative file paths
	Hmisc,       # misc data analysis utilities
	MASS,        # classical stats (e.g. mvrnorm)
	stats,       # base R stats
	splines,     # spline basis functions
	MCMCpack,    # MCMC samplers
	BayesTree,   # original BART implementation
	dbarts,      # discrete BART (faster BART)
	sf,          # spatial data; replacement for {rgdal}
	xtable,      # LaTeX/HTML tables
	SoftBart,    # Soft BART
	tidyverse,   # data wrangling + ggplot2
	brms,        # Bayesian regression via Stan
	invgamma     # inverse-gamma distribution
)

# ---- GitHub-only packages ----
# XBCF is not on CRAN; p_load_gh installs from the maintainer's repo if missing.
pacman::p_load_gh("JingyuHe/XBCF")

####
# Source self-written functions ----
####
# Pulls in every .R file under 01_code/00_functions/{01_helper, 02_methods}
# so all helpers and method implementations are available downstream.
invisible(lapply(
	unlist(lapply(
		c(
			here::here("01_code", "00_functions", "01_helper"),
			here::here("01_code", "00_functions", "02_methods")
		),
		list.files, pattern = "\\.R$", full.names = TRUE
	)),
	source
))
