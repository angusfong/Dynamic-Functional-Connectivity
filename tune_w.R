tune_w <- function(excl_sub, filt="", taper="", measure="", algo="plsr") {

	ws <- c("10","20", "30","40","50","60")
	rs <- rep(0, length(ws))

	for (w in 1:length(ws)) {

		# get training (now n-1)

		load(paste(ifelse(train_on=="task","t25_","r25_"),ws[w],filt, taper,measure,hm_resids,".RData", sep=""))
		this.w.train <- connectomes[,,-excl_sub]
		this.w.train.vect <- apply(this.w.train,3,function(mat) as.vector(mat))[upper,]

		load(paste(ifelse(test_on=="task","t25_","r25_"),ws[w],filt, taper,measure,hm_resids,".RData", sep=""))
		this.w.test <- connectomes[,,-excl_sub]
		this.w.test.vect <- apply(this.w.test,3,function(mat) as.vector(mat))[upper,]
		

		lr <- loo_v(this.w.train.vect, this.w.test.vect, behav[-excl_sub], ncomp=1)
		rs[w] <- lr

		#cat("Iteration", excl_sub, ": w =", ws[w], ": r =",lr[1],"\n")
		
	}

	best_w <- ws[which.max(as.numeric(rs))]
	
	cat("Iteration", excl_sub, ": best w =", best_w,"\n")
	return(as.numeric(best_w))
}	
