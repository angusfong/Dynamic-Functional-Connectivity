# compares anatomy of edges from sFC and dFC models
library(corrplot)
source("node2network.R")
source("node2lobe.R")
source("plotMat.R")

networks <- c("SubC", "MT", "MF", "VI", "VII", "VA", "DM", "FP")
lobes <- c("pf", "m", "ins", "par", "temp", "occ", "limb", "cer", "sub", "stem")

opt <- "d"
train_on <- "task"
test_on <- "task"

f <- paste(opt, train_on, test_on, "_edges.RData", sep="")
pdf(paste(opt, train_on, test_on, "_edges.pdf", sep=""))

load(f)
pos_edges <- sig_edges==1 & signs==1
neg_edges <- sig_edges==1 & signs==-1

# network behavior
# pos_edges
n_node <- 268
m <- matrix(0, n_node, n_node)

m[which(upper.tri(m))] <- pos_edges # can examine others
m <- m + t(m)
`diag<-`(m, 1)

d_pos_network <- node2network(m, avg=FALSE)
rownames(d_pos_network) <- networks
colnames(d_pos_network) <- networks
d_pos_network


d_pos_lobe <- node2lobe(m, avg=FALSE)
rownames(d_pos_lobe) <- lobes
colnames(d_pos_lobe) <- lobes
d_pos_lobe

plotMat(d_pos_network, "Dynamic positive edges by network", col_scale=c("white", "dark blue"))
plotMat(d_pos_lobe, "Dynamic positive edges by lobe")


# neg_edges
n_node <- 268
m <- matrix(0, n_node, n_node)

m[which(upper.tri(m))] <- neg_edges # can examine others
m <- m + t(m)
`diag<-`(m, 1)

d_neg_network <- node2network(m, avg=FALSE)
rownames(d_neg_network) <- networks
colnames(d_neg_network) <- networks
d_neg_network


d_neg_lobe <- node2lobe(m, avg=FALSE)
rownames(d_neg_lobe) <- lobes
colnames(d_neg_lobe) <- lobes
d_neg_lobe

plotMat(d_neg_network, "Dynamic negative edges by network", col_scale=c("white", "dark red"))
plotMat(d_neg_lobe, "Dynamic negative edges by lobe")

dev.off()
# static and dynamic overlap
d_pos_edges <- pos_edges
d_neg_edges <- neg_edges
f <- paste("s", train_on, test_on, "_edges.RData", sep="")

load(f)
s_pos_edges <- sig_edges==1 & signs==1
s_neg_edges <- sig_edges==1 & signs==-1

# static edges
n_node <- 268
sp <- matrix(0, n_node, n_node)
sn <- matrix(0, n_node, n_node)

sp[which(upper.tri(sp))] <- s_pos_edges # can examine others
sp <- sp + t(sp)
`diag<-`(sp, 1)

s_pos_network <- node2network(sp, avg=FALSE)
rownames(s_pos_network) <- networks
colnames(s_pos_network) <- networks
s_pos_network

sp <- node2lobe(sp, avg=FALSE)
rownames(s_pos_lobe) <- lobes
colnames(s_pos_lobe) <- lobes
s_pos_lobe

sn[which(upper.tri(sn))] <- s_neg_edges # can examine others
sn <- sn + t(sn)
`diag<-`(sn, 1)

s_neg_network <- node2network(sn, avg=FALSE)
rownames(s_neg_network) <- networks
colnames(s_neg_network) <- networks
s_neg_network


s_neg_lobe <- node2lobe(sn, avg=FALSE)
rownames(s_neg_lobe) <- lobes
colnames(s_neg_lobe) <- lobes
s_neg_lobe

plotMat(s_pos_network, "Static positive edges by network")
plotMat(s_pos_lobe, "Static positive edges by lobe")
plotMat(s_neg_network, "Static negative edges by network")
plotMat(s_neg_lobe, "Static negative edges by lobe")

# overlap
m <- matrix(0, n_node, n_node)
m[which(upper.tri(m))] <- s_neg_edges & d_neg_edges # can examine others
m <- m + t(m)
`diag<-`(m, 1)
neg_neg_network <- node2network(m, avg=FALSE)
rownames(neg_neg_network) <- networks
colnames(neg_neg_network) <- networks
plotMat(neg_neg_network, "sneg dneg")

m <- matrix(0, n_node, n_node)
m[which(upper.tri(m))] <- s_pos_edges & d_neg_edges # can examine others
m <- m + t(m)
`diag<-`(m, 1)
pos_neg_network <- node2network(m, avg=FALSE)
rownames(pos_neg_network) <- networks
colnames(pos_neg_network) <- networks
plotMat(pos_neg_network, "sneg dneg")

# alternative analysis with edge weights
m <- m_posonly <- m_negonly <- matrix(0, n_node, n_node)
# get m_pls from running combModel_ext.R up till model fitting
m[which(upper.tri(m))] <- m_posonly[which(upper.tri(m))] <- m_negonly[which(upper.tri(m))] <- m_pls$coefficients[,,3] # can examine others
m <- m + t(m)
m_posonly <- m_posonly + t(m_posonly)
m_negonly <- m_negonly + t(m_negonly)
`diag<-`(m, 1)
`diag<-`(m_posonly, 1)
`diag<-`(m_negonly, 1)
m_posonly[m <= 0] <- NA
m_negonly[m >= 0] <- NA
edges_network <- node2network(m, avg=FALSE)
edges_network_posonly <- node2network(m_posonly, avg=FALSE)
edges_network_negonly <- node2network(m_negonly, avg=FALSE)

colnames(edges_network) <- colnames(edges_network_posonly) <- colnames(edges_network_negonly) <- networks
rownames(edges_network) <- rownames(edges_network_posonly) <- rownames(edges_network_negonly) <- networks
pdf("dynamic_edge_weights.pdf")
plotMat(edges_network, "all")
plotMat(edges_network_posonly, "pos_only", col_scale=c("white", "dark blue"))
plotMat(edges_network_negonly, "neg_only", col_scale=c("dark red","white"))
dev.off()

