paper_dir <- "/gpfs/milgram/project/chun/hf246/Predictions/paper_code"
setwd(paper_dir)

args = commandArgs(trailingOnly=TRUE)
ptest <- FALSE #ncomp only 1:6; doesn't cat opt

opt <- args[1] # c/s/d/c1
train_on <- args[2]
test_on <- args[3]

#params
filt <- "" # from ""(no filt)/"_hpf"/"_hpf_lpf"
mot <- FALSE # t/f only
taper <- "" # from ""(no taper) / "_tuk0.x"
hm_resids <- "" # or ""

## ----setup, include=FALSE------------------------------------------------
library(R.matlab)
library(foreach)
library(doParallel)
library(pls)

registerDoParallel(25)

source("loo_v.R")
source("tune_w.R")
source("tune_ncomp.R")
source("comb_predict.R")
## ------------------------------------------------------------------------

#presets
algo <- "plsr"
cor_method <- "spearman"
pls_method <- "kernel"
ncomp_free <- TRUE

setwd(paper_dir)

cat("**Model:", opt, "\n")
cat("**Training on", train_on, ", testing on", test_on, "\n")
if (mot) cat("\nEmploying feature-level motion control!!\n\n")
cat("Using", algo, "with", filt, taper, hm_resids,"params!\n")

n_node <- 268
thresh <- 0.01
mat_env <- readMat("Rosenberg2016_attn_data_share.mat")
behav <- as.vector(mat_env$dprime)

ws <- c("10","20","30","40","50","60")

## ------------------------------------------------------------------------
aa <- matrix(data=NA, n_node, n_node)
upper <- which(upper.tri(aa))
n_edge <- length(upper)
## ------------------------------------------------------------------------
#load static matrices

if (hm_resids=="_hmresiduals") {
	load(paste("static_",train_on,"_corrected.RData",sep=""))
	s.train.mats <- get(paste("static_",train_on,"_corrected",sep=""))
	load(paste("static_",test_on,"_corrected.RData",sep=""))
	s.test.mats <- get(paste("static_",test_on,"_corrected",sep=""))
} else {
	s.train.mats <- mat_env[paste("gradCPT.", train_on, ".mats", sep="")][[1]]
	s.test.mats <- mat_env[paste("gradCPT.", test_on, ".mats", sep="")][[1]]
}


s.train.vects <- apply(s.train.mats,3,function(mat) as.vector(mat))[upper,]
s.test.vects <- apply(s.test.mats,3,function(mat) as.vector(mat))[upper,]

n_sub <- dim(s.train.mats)[3]
cat(n_sub, "subjects in total \n")
n_train_sub <- n_sub-1
## ------------------------------------------------------------------------

if (mot) { 
#then no tune w to get mats, already 20. just run this part w/ variables to see # edges excluded
#feeds processed matrices to prediction

	load(paste(ifelse(train_on=="task","t25_","r25_"),"20",filt,taper,hm_resids,".RData", sep=""))
	d.train.mats <- connectomes

	load(paste(ifelse(test_on=="task","t25_","r25_"),"20",filt,taper,hm_resids,".RData", sep=""))
	d.test.mats <- connectomes

	s.train.vects <- apply(s.train.mats,3,function(mat) as.vector(mat))[upper,]
	s.test.vects <- apply(s.test.mats,3,function(mat) as.vector(mat))[upper,]
	d.train.vects <- apply(d.train.mats,3,function(mat) as.vector(mat))[upper,]
	d.test.vects <- apply(d.test.mats,3,function(mat) as.vector(mat))[upper,]

	bad_s_edges <- NULL
	bad_d_edges <- NULL

	m <- t(readMat("gradCPT_motion.mat")[[1]])
	for (i in 1:nrow(s.train.vects)) {
		for (j in 1:nrow(m)) {
			if (cor.test(s.train.vects[i,], m[j,], method=cor_method)$p.value < 0.05 | cor.test(s.test.vects[i,], m[j,], method=cor_method)$p.value < 0.05) {
				bad_s_edges <- c(bad_s_edges,i)
				break
			}
		}
	}
	s.train.vects <- s.train.vects[-bad_s_edges,]
	s.test.vects <- s.test.vects[-bad_s_edges,]
	for (i in 1:nrow(d.train.vects)) {
		for (j in 1:nrow(m)) {
			if (cor.test(d.train.vects[i,], m[j,], method=cor_method)$p.value < 0.05 | cor.test(d.test.vects[i,], m[j,], method=cor_method)$p.value < 0.05) {
				bad_d_edges <- c(bad_d_edges,i)
				break
			}
		}
	}
	d.train.vects <- d.train.vects[-bad_d_edges,]
	d.test.vects <- d.test.vects[-bad_d_edges,]

	cat(length(bad_s_edges), "static and", length(bad_d_edges), "dynamic edges were excluded for motion.\n")
}



system.time({
  lps <- foreach (excl_sub=1:n_sub, .combine='rbind') %dopar%
    do.call(
      function(excl_sub) {

	#learn optimal w for this subject based on RR/TT/RT/TR

	# best_w <- pre_ws[excl_sub]

	if (opt != 's') best_w <- tune_w(excl_sub, filt=filt, taper=taper) else best_w <- 20 # dummy, doesn't matter

	if (!mot) {
		load(paste(ifelse(train_on=="task","t25_","r25_"),best_w,filt,taper,hm_resids,".RData", sep=""))
		d.train.mats <- connectomes

		load(paste(ifelse(test_on=="task","t25_","r25_"),best_w,filt,taper,hm_resids,".RData", sep=""))
		d.test.mats <- connectomes


		d.train.vects <- apply(d.train.mats,3,function(mat) as.vector(mat))[upper,]
		d.test.vects <- apply(d.test.mats,3,function(mat) as.vector(mat))[upper,]
	}
	
	return(comb_predict(opt=opt, s.train.vects, s.test.vects, d.train.vects, d.test.vects, behav, excl_sub=excl_sub, ncomps=ncomps))

      }, list(excl_sub)
    )
})

# plots and model evaluation
preds <- unlist(lps)
print(preds)
corr <- cor(preds, behav, method=cor_method)
#p <- cor.test(preds, behav, method=cor_method)$p.value	
cat("Final r from PLSR =", corr, "\n")
#model_f <- paste(opt, "ModelPLSR",ncomps, train_on, "_", test_on, filt,taper,filt,"plots.pdf", sep="")
data_f <- paste(opt, train_on, test_on, "preds.RData", sep="")
save("preds", "behav", file=data_f)
#pdf(file=model_f)
#plot(preds ~ behav, xlab="Observed dprime", ylab="Predicted dprime", 
#	main=paste("PLSR Combined Model:", train_on, "on", test_on, ": r =", round(corr,2), "p =", round(p, 2)))
#abline(lm(preds ~ behav), col="red")
#dev.off()

