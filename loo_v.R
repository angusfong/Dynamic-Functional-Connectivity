loo_v <- function(train.vect, test.vect, behav.loo, ncomp=5) {
# runs CPM from train on test, returns r between predicted and observed
# assumes vectorized, and upper triangle of the covariance matrix

	cor_method <- "spearman"
    n_subl <- dim(train.vect)[2]
		
    loo.preds <- rep(0, n_subl)
    for (excl_subl in 1:n_subl) {

	    upper_vect <- train.vect[,-excl_subl]

		test_upper_vect <- test.vect[,excl_subl]

		# note that behav.loo already leaves out the OOS point
		m_pls <- plsr(behav.loo[-excl_subl] ~ t(upper_vect), ncomp=ncomp, method=pls_method)
		loo.preds[excl_subl] <- predict(m_pls, ncomp=ncomp, newdata=t(test_upper_vect))

	}
    

    corr <- cor(loo.preds, behav.loo, method=cor_method)
    return(corr)
	
	    

}

