tune_ncomp <- function(excl_sub, s.train.vects, s.test.vects, d.train.vects, d.test.vects, c.upper_vect, c.test_upper_vect, behav) {

	if (!ptest) cat("**Tuning ncomp: subject", excl_sub, "\n")

	# best ncomp for static
	if (opt=="s") {
		if (ptest) ns <- 1:6 else ns <- 1:(n_train_sub-2)
		rs <- rep(0, length(ns))
		for (n in ns) {
			
			lr <- loo_v(s.train.vects[,-excl_sub], s.test.vects[,-excl_sub], behav[-excl_sub], ncomp=n)
			rs[n] <- lr

		}

		s.best_ncomp <- ns[which.max(rs)]
		s.best_r <- max(rs)
		if (!ptest) cat("Subject", excl_sub, "Performance across ncomps:", rs, "; Optimal s.ncomp =", s.best_ncomp, ", r =", s.best_r, ".\n\n")
		return(list("s"=s.best_ncomp))

	}

	# best ncomp for dynamic
	if (opt=="d") {
		if (ptest) ns <- 1:6 else ns <- 1:(n_train_sub-2)
		rs <- rep(0, length(ns))
		for (n in ns) {

			lr <- loo_v(d.train.vects[,-excl_sub], d.test.vects[,-excl_sub], behav[-excl_sub], ncomp=n)
			rs[n] <- lr
		
		}	

		d.best_ncomp <- ns[which.max(rs)]
		d.best_r <- max(rs)
		if (!ptest) cat("Subject", excl_sub, "Performance across ncomps:", rs, "; Optimal d.ncomp =", d.best_ncomp, ", r =", d.best_r, ".\n\n")

		return(list("d"=d.best_ncomp))

	}

	# best ncomp for combined
	if (opt=="c") {
		if (ptest) ns <- 1:6 else ns <- 1:(n_train_sub-2)
		rs <- rep(0, length(ns))
		for (n in ns) {
			
			lr <- loo_v(c.upper_vect, c.test_upper_vect, behav[-excl_sub], ncomp=n)
			rs[n] <- lr

		}

		c.best_ncomp <- ns[which.max(rs)]
		c.best_r <- max(rs)
		if (!ptest) cat("Subject", excl_sub, "Performance across ncomps:", rs, "; Optimal c.ncomp =", c.best_ncomp, ", r =", c.best_r, ".\n\n")

		return(list("c"=c.best_ncomp))

	}

}	
