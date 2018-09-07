task <- "SST"
opt <- "d" # c/s/d/c1
train_on <- "task"
test_on <- "task"

stride4 <- FALSE


if (train_on == "rest") ncomp <- 1 else {
if (test_on=="rest") {
	if (opt=="d") ncomp <- 4 else ncomp <- 2
} else {
	if (opt=="c") ncomp <- 1 else ncomp <- 3}}
ncomp

library(R.matlab)
library(pls)

setwd("/gpfs/milgram/project/chun/hf246/Predictions")

# CoV or Var?
cor_method <- "spearman"

cat(opt, "model...\n")
full_subs <- TRUE

cat(task, "dataset...\n")
algo <- "plsr"
cat(algo, ncomp, "components...\n")
measure <- "" # from ""(var)/"_sd"/"_cv"
filt <- "" # from "_hpf", "_hpf_lpf"
mot <- "" # from ""(no excl)/"_mot_excl"
mot_thresh <- 0.05
taper <- "_tuk0.5" # from "" / "_tuk0.2" / "_tuk0.5"
cat("Using", measure, filt, taper, "!\n")

# External validation


cat(train_on, "on", test_on, "...\n")

if (algo=="lr") {# var --> 29,20,20,28 (if no hpf)
    if (train_on=="rest" & test_on=="rest") best_w=29
    if (train_on=="rest" & test_on=="task") best_w=20
    if (train_on=="task" & test_on=="rest") best_w=20
    if (train_on=="task" & test_on=="task") best_w=28
} else

if (filt == "_hpf") {# var --> 21,22,21,26 (if hpf)
    if (train_on=="rest" & test_on=="rest") best_w=21
    if (train_on=="rest" & test_on=="task") best_w=22
    if (train_on=="task" & test_on=="rest") best_w=21
    if (train_on=="task" & test_on=="task") best_w=26
}

if (algo == "plsr") { #(if no hpf)
    if (train_on=="rest" & test_on=="rest") best_w=20
    if (train_on=="rest" & test_on=="task") best_w=20
    if (train_on=="task" & test_on=="rest") best_w=20
    if (train_on=="task" & test_on=="task") best_w=20
} 

# 0- the train matrices
setwd("/gpfs/milgram/project/chun/hf246/Predictions")
load(paste(ifelse(train_on=="task","t25_","r25_"),best_w,taper,measure,filt,
	ifelse(stride4,"_stride4",""),".RData", sep=""))
d.train.mats <- aa

mat_env <- readMat("Rosenberg2016_attn_data_share.mat")
s.train.mats <- mat_env[paste("gradCPT.", train_on, ".mats", sep="")][[1]]

# 1- the train behavior
train.behav <- as.vector(mat_env$dprime)
n_train_sub <- length(train.behav)

# 2- the test matrices
setwd("/gpfs/milgram/project/chun/hf246/Predictions")
load(paste(ifelse(train_on=="task", "t", "r"),ifelse(test_on=="task", "t", "r"),
	task,algo,taper,filt,ifelse(stride4,"_stride4",""),".RData", sep=""))
d.test.mats <- aa # usually, aa
rm("aa")

if (task=="ANT") {
	ANTmatlab <- readMat("ANT_behav_share.mat")
	s.test.mats <- ANTmatlab[paste("ANT.", test_on, ".mats", sep="")][[1]]
} else

if (task=="ADHD") {
	s.test.mats <- mat_env$ADHD.mats.268

    missing_ROIs <- c(51,55,57,60,104,107,108,109,112,115,116,118,129,131,135,136,185,
		      189,194,196,202,214,239,240,243,246,249,250,252,256,266,268)

	s.train.mats <- s.train.mats[-missing_ROIs,-missing_ROIs,]
	d.train.mats <- d.train.mats[-missing_ROIs,-missing_ROIs,]
	s.test.mats <- s.test.mats[-missing_ROIs,-missing_ROIs,]
	d.test.mats <- d.test.mats[-missing_ROIs,-missing_ROIs,]
} else

if (task=="SST") {
	setwd("MPH/")
	MPHenv <- readMat("MethylphenidateProjectData.mat")
	if (test_on=="task") goodMPH <- which(MPHenv$incltask==1) else goodMPH <- which(MPHenv$inclrest==1)
	# goodMPH <- which(MPHenv$incltask==1 & MPHenv$inclrest==1) # testing only
	s.test.mats <- MPHenv[paste("MPH.", test_on, ".mats", sep="")][[1]][,,goodMPH]
	d.test.mats[which(is.na(d.test.mats))] <- 0 # removes NAs in dynamic matrix
}


# 3- the test behavior
if (task=="ANT") {
    test.behavs <- readMat("ANT_behav.mat")[[1]]
    test.behav <- test.behavs[,3] / test.behavs[,2] # col3/col2=RT var; col1=accuracy; col4,5,6=alert,orient,exec ctrl
} else if (task=="ADHD") test.behav <- as.vector(readMat("ADHD/ADHD_behav.mat")[[1]]) else
if (task=="SST") test.behav <- as.vector(MPHenv$behav[goodMPH,1])

n_test_sub <- length(test.behav)

# 4- train a full PLSR or GLM model

upper <- which(upper.tri(s.train.mats[,,1]))
s.train.vect <- apply(s.train.mats, 3, function(mat) as.vector(mat))[upper,]
d.train.vect <- apply(d.train.mats, 3, function(mat) as.vector(mat))[upper,]
s.test.vect <- apply(s.test.mats, 3, function(mat) as.vector(mat))[upper,]
d.test.vect <- apply(d.test.mats, 3, function(mat) as.vector(mat))[upper,]

if (opt=="c") {
	train_upper_vect <- t(rbind(s.train.vect, d.train.vect)) # n_sub by n_edge
	test_upper_vect <- t(rbind(s.test.vect, d.test.vect))
} else
if (opt=="s") {
	train_upper_vect <- t(s.train.vect) # n_sub by n_edge
	test_upper_vect <- t(s.test.vect)
} else
if (opt=="d") {
	train_upper_vect <- t(d.train.vect) # n_sub by n_edge
	test_upper_vect <- t(d.test.vect)
}

# PLSR
if (algo=="plsr") {
	m_pls <- plsr(train.behav ~ train_upper_vect, ncomp=ncomp, method="kernel")
	preds <- predict(m_pls, ncomp=ncomp, newdata=test_upper_vect)
} 

#5 - evaluate performance
corr <- cor(preds, test.behav, method=cor_method, use="complete.obs")
if (task=="ANT") corr <- corr * -1
p <- cor.test(preds, test.behav, method=cor_method)$p.value
cat("Final r from pls model:", corr, ", p =", p, "\n")

#6 - plot
# plotting?
f <- paste("combModel_ext(",task,ifelse(task=="ANT", "-RV",""), ")_", train_on, "_", test_on, measure,filt,mot,"plots.pdf", sep="")
pdf(file=f)
plot(preds ~ test.behav, xlab=paste("Observed",task), ylab="Predicted dprime", 
	main=paste("Combined Model",measure,filt,mot,taper, "(",task,ifelse(task=="ANT", "-RV",""),"):", train_on, "on", test_on, ": r =", round(corr, 2), "p =", round(p, 2)))
abline(lm(preds ~ test.behav), col="red")
dev.off()

for (i in 1:10) {
corrs <- rep(NA, 1000)
for (i in 1:1000) {corrs[i] <- cor(preds, sample(test.behav))}
cat(mean(corrs>=corr), "\n")
}


setwd("/gpfs/milgram/project/chun/hf246/Predictions")
