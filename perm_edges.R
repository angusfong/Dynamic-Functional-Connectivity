# One tailed permutation test to find significant edges
nperm <- 10^5
thres <- 0.05
set.seed(123)

library(R.matlab)
library(pls)


setwd("/nexsan/home01/chun/hf246/Lab/Predictions")

opt <- "c" # "c" / "s" / "d"
cat(opt, "model...\n")

ncomp=2
cat(ncomp, "components...\n")

taper <- "_tuk0.5" # from "" / "_tuk0.2" / "_tuk0.5"
cat("Using", taper, "taper!\n")

# External validation
train_on <- "task" # internal sample
test_on <- "rest" # external sample

cat(train_on, "on", test_on, "...\n")

best_w = 20

# 0- the train matrices
setwd("/nexsan/home01/chun/hf246/Lab/Predictions")
load(paste(ifelse(train_on=="task","t25_","r25_"),best_w,taper,".RData", sep=""))
d.train.mats <- aa

mat_env <- readMat("Rosenberg2016_attn_data_share.mat")
s.train.mats <- mat_env[paste("gradCPT.", train_on, ".mats", sep="")][[1]]

# 1- the train behavior
train.behav <- as.vector(mat_env$dprime)
n_train_sub <- length(train.behav)

# 2- train a full PLSR model

upper <- which(upper.tri(s.train.mats[,,1]))
s.train.vect <- apply(s.train.mats, 3, function(mat) as.vector(mat))[upper,]
d.train.vect <- apply(d.train.mats, 3, function(mat) as.vector(mat))[upper,]


if (opt=="c") {
	train_upper_vect <- t(rbind(s.train.vect, d.train.vect)) # n_sub by n_edge
} else
if (opt=="s") {
	train_upper_vect <- t(s.train.vect) # n_sub by n_edge
} else
if (opt=="d") {
	train_upper_vect <- t(d.train.vect) # n_sub by n_edge
}

# PLSR
m_pls <- plsr(train.behav ~ train_upper_vect, ncomp=ncomp)
stats <- m_pls$coefficients[,,ncomp]
stats_abs <- abs(stats)
signs <- sign(stats)
n_edge <- length(stats)
count_edges <- rep(0, n_edge)

for (perm in 1:nperm) {
	if (perm %% 100 == 0) cat("Permutation", perm, "...\n")
	perm.behav <- sample(train.behav)
	m_perm <- plsr(perm.behav ~ train_upper_vect, ncomp=ncomp)
	perm_edges <- m_perm$coefficients[,,ncomp]
	#count_edges is here # of permutations with greater coefficient (same sign) than stat for each edge
	count_edges <- count_edges + as.numeric(abs(perm_edges) > stats_abs & sign(perm_edges) == sign(stats))
}

#sig_edges is a mask of significant edgess
sig_edges <- count_edges / nperm < thres

#save edges
save(sig_edges, stats, signs, file=paste(opt, train_on, test_on, "_edges.RData", sep=""))