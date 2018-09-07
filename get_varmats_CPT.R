# saves variance matrices to a database to avoid future computation
# file is t/r_22/20/18_w

args = commandArgs(trailingOnly=TRUE)

library(R.matlab)
library(plotrix)
library(corrplot)
library(glasso)
library(psych)
library(mpmi)
library(abind)
library(spectral)
library(signal)

# can write as a loop
type <- args[1]
cat(type, "matrices", "\n")
window_size <- as.numeric(args[2])
cat("w = ", window_size, "\n")

#params
tukey <- FALSE
a=0.5
taper <- win.tukey(rep(1,window_size),a=a)

filt <- ""

if (filt=="_hpf") {
  library(signal)
  crit <- 1/window_size
  bf <- butter(2, crit, "high")
}

if (tukey) cat("Using tukey windows with a =", a, "\n")

#loop through subjects to compute matrices one by one for SST... god damn it
setwd("/gpfs/milgram/project/chun/hf246/Predictions/roimean_fullruns")
all_mats <- list.files()

n_sub <- 25

if (type=="task") curr.mats <- array(NA, c(2448, 268, n_sub)) else curr.mats <- array(NA, c(718, 268, n_sub))

for (i in 1:n_sub) {
	mtx <- as.matrix(read.table(all_mats[i*2 - (type=="rest")], skip=1, fill=TRUE))[,-1]
	#if (filt=="_hpf") mtx_f <- apply(mtx, 2, function(x) filter(bf, x))
	curr.mats[1:nrow(mtx),,i] <- mtx
}

if (filt=="_hpf") {
	if (type=="task") {
		block1 <- apply(curr.mats[c(1:180),,],c(2,3), function(x) filter(bf, x))
		block2 <- apply(curr.mats[c(213:392),,],c(2,3), function(x) filter(bf, x))
		block3 <- apply(curr.mats[c(425:604),,],c(2,3), function(x) filter(bf, x))
		block4 <- apply(curr.mats[c(637:816),,],c(2,3), function(x) filter(bf, x))
		block5 <- apply(curr.mats[c(817:996),,],c(2,3), function(x) filter(bf, x))
		block6 <- apply(curr.mats[c(1029:1208),,],c(2,3), function(x) filter(bf, x))
		block7 <- apply(curr.mats[c(1241:1420),,],c(2,3), function(x) filter(bf, x))
		block8 <- apply(curr.mats[c(1453:1632),,],c(2,3), function(x) filter(bf, x))
		block9 <- apply(curr.mats[c(1633:1812),,],c(2,3), function(x) filter(bf, x))
		block10 <- apply(curr.mats[c(1845:2024),,],c(2,3), function(x) filter(bf, x))
		block11 <- apply(curr.mats[c(2057:2236),,],c(2,3), function(x) filter(bf, x))
		block12 <- apply(curr.mats[c(2269:2448),,],c(2,3), function(x) filter(bf, x))
		run1 <- abind(block1,block2,block3,block4,along=1)
		run2 <- abind(block5,block6,block7,block8,along=1)
		run3 <- abind(block9,block10,block11,block12,along=1)
	} else {
		run2 <- run1 <- array(NA, c(359,268,n_sub))
		for (sub in 1:n_sub) {
			thisrun1 <- apply(curr.mats[c(1:359),,sub],2, function(x) filter(bf, x))
			thisrun2 <- apply(curr.mats[c(360:718),,sub],2, function(x) filter(bf, x))

			run1[1:nrow(thisrun1),,sub] <- thisrun1
			run2[1:nrow(thisrun2),,sub] <- thisrun2
		}
	}

} else {
	if (type=="task") {
		run1 <- curr.mats[c(1:180, 213:392, 425:604, 637:816),,] 
	    run2 <- curr.mats[c(817:996, 1029:1208, 1241:1420, 1453:1632),,]
	    run3 <- curr.mats[c(1633:1812, 1845:2024, 2057:2236, 2269:2448),,]
	} else {
		run1 <- curr.mats[1:359,,]
		run2 <- curr.mats[360:718,,]
	}
}

## ------------------------------------------------------------------------
n_node <- 268
mats_dur <- dim(curr.mats)[1]
n_windows <- mats_dur - window_size + 1
n_windows_run <- dim(run1)[1] - window_size + 1
## ------------------------------------------------------------------------
tmp <- matrix(data=NA, n_node, n_node)
upper <- which(upper.tri(tmp))
n_edge <- length(upper)
## ------------------------------------------------------------------------

connectomes <- array(0, c(n_node,n_node,n_sub))
for (sub in 1:n_sub) {
	
	cat("Now subject", sub, "...\n")
	run1.this.sub <- run1[,,sub]
	run2.this.sub <- run2[,,sub]
	if (type=="task") run3.this.sub <- run3[,,sub]

	run1_upper_matx <- matrix(NA, n_windows_run, n_edge)
	run2_upper_matx <- matrix(NA, n_windows_run, n_edge)
	if (type=="task") run3_upper_matx <- matrix(NA, n_windows_run, n_edge)
    
    #simultaneously across 3 runs
    for (w in 1:n_windows_run) {

		run1.now <- run1.this.sub[w:(w-1+window_size),]
		run2.now <- run2.this.sub[w:(w-1+window_size),]
		if (type=="task") run3.now <- run3.this.sub[w:(w-1+window_size),]

		if (tukey) run1.now <- apply(run1.now, 2, function(x) x*taper)
		if (tukey) run2.now <- apply(run2.now, 2, function(x) x*taper)
		if (type=="task") if (tukey) run3.now <- apply(run3.now, 2, function(x) x*taper)
	    
	    v1 <- fisherz(cor(run1.now)[upper])
	    v2 <- fisherz(cor(run2.now)[upper])
	    if (type=="task") v3 <- fisherz(cor(run3.now)[upper])
	    
		run1_upper_matx[w,] <- v1
		run2_upper_matx[w,] <- v2
		if (type=="task") run3_upper_matx[w,] <- v3
	    
    }


	upper_matx <- rbind(run1_upper_matx, run2_upper_matx)
	if (type=="task") upper_matx <- rbind(upper_matx, run3_upper_matx)
	varmats <- apply(upper_matx, 2, function(x) var(x, na.rm=TRUE))

	connectomes[,,sub][upper.tri(connectomes[,,i])] <- varmats
	connectomes[,,sub] <- connectomes[,,sub] + t(connectomes[,,sub])
		
}

setwd("/gpfs/milgram/project/chun/hf246/Predictions/paper_code")
if (tukey) tuk_label <- paste("_tuk",a,sep="") else tuk_label <- ""
cat("Saving to:", paste(ifelse(type=="task","t","r"),"25_", window_size, filt, tuk_label,".RData",sep=""), "\n")
save(connectomes, file=paste(ifelse(type=="task","t","r"),"25_", window_size, filt, tuk_label,".RData",sep=""))


