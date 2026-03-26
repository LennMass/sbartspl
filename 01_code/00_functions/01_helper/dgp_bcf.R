

dgp_bcf <-function(n = 500, # sample size
									 a = 0.1, # a parameter in non-overlap definition#
									 b = 7, # b parameter in non-overlap definition#
									 kappa = 0.5,
									 rho = 1,
									 tau_setting="homogeneous", 
									 p_cont_noise = 10, # 10, 20
									 p_cat_noise = 10 # 10, 20
									 ){
	
	g <- function(x) {
		ret <- rep(0, length(x))
		ret[x == 1] <- 2
		ret[x == 2] <- -1
		ret[x == 3] <- -4
		return(ret)
	}
	
	h <- function(x, p_cont_true, p_cont_noise) {
		eta <- 
			2 * (x[,1]^2 - 0.5)^2 -      # quadratic peak/valley along x1
			1.5 * (x[,2]^2 - 0.2)^2 +    # another peak/valley along x2
			0.8 * sin(pi * x[,3]) +       # oscillatory pattern along x3
			pnorm(x[,p_cont_true+p_cont_noise+1] - x[,p_cont_true+p_cont_noise+2])
		
		ps <- plogis(eta)               # map to (0,1)
		
		
		return(ps)
	}
	
	pi_func <- function(x, p_cont_true, p_cont_noise) {
		res <- h(x, p_cont_true, p_cont_noise)
		res <- pmin(pmax(res, 0.01), 0.99)
		return(res)
	}
	
	# non-linear
	mu <- function(x, p_cont_true, p_cont_noise) {
		1 + g(x[, p_cont_true+p_cont_noise+1]) + x[, 1] * x[, 3]
	}
	
	
		
	if (tau_setting == "homogeneous") {
		
		# homogenerous
		tau <- function(x, p_cont_true, p_cont_noise) {
			rep(3, nrow(x))
		}
	} else if (tau_setting == "heterogeneous") {
		
		# heterogeneous
		tau <- function(x, p_cont_true, p_cont_noise) {
			1 + 2 * x[, 2] * x[, p_cont_true+p_cont_noise+2]
		}
	} else {
		stop("Choose correct 'tau_setting' setting.")
	}
		
		
	generate_x <- function(n=500, p_cont = 3, p_cat = 4, cat_levels = list(c(0,1), 1:3)) {
		# p_cont: number of continuous covariates
		# p_cat: number of categorical covariates (at the end)
		# cat_levels: list of possible values for each categorical covariate
		
		# 1. Create continuous covariates
		x_cont <- matrix(rnorm(n * p_cont), n, p_cont)
		
		# 2. Create categorical covariates
		x_cat <- matrix(0, n, p_cat)
		for (j in 1:p_cat) {
			select_cat_level <- sample(seq(1:length(cat_levels)), 1)
			x_cat[, j] <- sample(cat_levels[[select_cat_level]], n, replace = TRUE)
		}
		
		# 3. Combine into full matrix: continuous first, categorical at the end
		x <- cbind(x_cont, x_cat)
		
		return(x)
	}
	
	max_iter <- 1000
	iter <- 0
	
	repeat {
		iter <- iter + 1
		
		
		p_cont_true = 3 # true continuous vars
		p_cat_true = 2  # true categorical vars
		x <- generate_x(n=n,
										p_cont = p_cont_true + p_cont_noise,
										p_cat  = p_cat_true + p_cat_noise)
		
		
		ce_ps <- pi_func(x, p_cont_true, p_cont_noise)
		z <- rbinom(n, 1, ce_ps) 
		
		# draw outcome
		f_xz <- mu(x, p_cont_true, p_cont_noise) + rho * tau(x, p_cont_true, p_cont_noise) * z
		sigma <- kappa * sd(rho * f_xz)
		y1 <- mu(x, p_cont_true, p_cont_noise) + rho * tau(x, p_cont_true, p_cont_noise)
		y0 <- mu(x, p_cont_true, p_cont_noise)
		y <- f_xz + sigma * rnorm(n)
		
		# calculate the true average treatment effect (ATE)
		ate_true <- mean(tau(x, p_cont_true, p_cont_noise)) 
		
		ce_true <- y1 - y0
		ce_trt <- z
		
		sink("/dev/null")
		fitz <- nnet::nnet(z ~ ., data = cbind(z, x), size = 3, rang = 0.1, maxit = 1000, abstol = 1.0e-8, decay = 5e-2)
		sink() # close the stream
		ps <- fitz$fitted.values
		
		## GR check ##
		ce_gr<-gr(Y=y,
							trt=z,
							ps=ps,
							X=x,
							M=500,
							qps=quantile(ps,probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))
		check_GR<-ce_gr[[1]]
		
		
		## find RO and RN ##
		RO_cor <- pw_overlap(ps = ps, E = z, a = a, b = b)
	
		order_ps_cor <- ps[order(ps)]
		order_RO_cor <- RO_cor[order(ps)]
		ll_cor <- max(which(order_RO_cor == 1))
		
		## ignore small non-overlap in the left tail (add all left tail to the RO) ##
		temp1 <- order_ps_cor[1:ll_cor]
		ps_addRO <- temp1[which(order_RO_cor[1:ll_cor] == 0)]
		RO_cor[which((ps %in% ps_addRO) == 1)] <- 1
		
		RO_mean <- mean(RO_cor)
		
		
		
		# safeguards to generate relevant overlap cases
		if( (RO_mean >= 0.70) && (RO_mean <=0.90) && check_GR==1) {
			
			break
			
		}
		
		if (iter >= max_iter) {
			stop("Timeout: condition not met after ", max_iter, " iterations")
		}
		
	}
	
	
	

	## create untrimmed and trimmed dataset ##
	datall_cor <- data.frame(y, z, ps, x)
	dattr_cor <- datall_cor[which(RO_cor == 1), ]
	
	
	
	return(list(untrimmed_dat = datall_cor,
							trimmed_dat = dattr_cor,
							ate_true=ate_true,
							RO=RO_cor,
							RO_share=RO_mean, 
							p_cont_noise = p_cont_noise, 
							p_cat_noise = p_cat_noise
							))
	
}






