comb_predict <- function(opt, s.train.vects, s.test.vects, d.train.vects, d.test.vects, behav, excl_sub, ncomps=1) {
# glm with static pos, neg, glm, net; dynamic neg    

	n_sub <- dim(s.train.mats)[3]
	n_train_sub <- n_sub-1

	s.upper_vect <- s.train.vects[,-excl_sub]
	d.upper_vect <- d.train.vects[,-excl_sub]
	s.other_upper_vect <- s.test.vects[,-excl_sub] #other means all subject, but test type (same if RR/TT)
	d.other_upper_vect <- d.test.vects[,-excl_sub]
	s.test_upper_vect <- s.test.vects[,excl_sub]
	d.test_upper_vect <- d.test.vects[,excl_sub]

	# make train and test vectors for d, s, and c (concatenated)
	c.upper_vect <- rbind(s.upper_vect, d.upper_vect)
	c.other_upper_vect <- rbind(s.other_upper_vect, d.other_upper_vect)
	c.test_upper_vect <- c(s.test_upper_vect, d.test_upper_vect)

	n_edge <- nrow(s.upper_vect)

	# tune ncomp
	if (ncomp_free) best_ncomp <- tune_ncomp(excl_sub, s.train.vects, s.test.vects, d.train.vects, 
		d.test.vects, c.upper_vect, c.other_upper_vect, behav) else best_ncomp = list("s"=ncomps, "d"=ncomps, "c"=ncomps)

	if (opt=="s") {
		s.m_pls <- plsr(behav[-excl_sub] ~ t(s.upper_vect), ncomp=best_ncomp$s, method=pls_method)
		return(predict(s.m_pls, ncomp=best_ncomp$s, newdata=t(s.test_upper_vect)))
	} else
	if (opt=="d") {
		d.m_pls <- plsr(behav[-excl_sub] ~ t(d.upper_vect), ncomp=best_ncomp$d, method=pls_method)
		return(predict(d.m_pls, ncomp=best_ncomp$d, newdata=t(d.test_upper_vect)))
	} else
	if (opt=="c") {
		c.m_pls <- plsr(behav[-excl_sub] ~ t(c.upper_vect), ncomp=best_ncomp$c, method=pls_method)
		return(predict(c.m_pls, ncomp=best_ncomp$c, newdata=t(c.test_upper_vect)))
	}	
	
}
