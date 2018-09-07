

node2network <- function(node_mtx, avg=FALSE)  {

library(R.matlab)

total_connections <- sum(node_mtx)

# either SUMS (eg for a mask) and then calculates proportion, or AVERAGES connections between nodes from canonical networks

parc <- as.vector(readMat("/gpfs/milgram/project/chun/hf246/Predictions/network_parcellation.mat")$new.label)

atlas <- list("SubC"=which(parc==1), "MT"=which(parc==2), "MF"=which(parc==3), "VI"=which(parc==4),
	      "VII"=which(parc==5), "VA"=which(parc==6), "DM"=which(parc==7), "FP"=which(parc==8))

n_node = nrow(node_mtx)
n_network = 8
     
#assume diag already set to 1     
     
network_mtx = matrix(0, n_network, n_network)
     
for (i in 1:n_network) {      
    for (j in 1:n_network) {
        thegrid <- expand.grid(atlas[[i]], atlas[[j]])
	connections <- apply(thegrid, 1, function(row) node_mtx[row[1], row[2]])
	if (avg) network_mtx[i,j] <- mean(connections) else network_mtx[i,j] <- sum(connections)
    }
}   

return(network_mtx)

}
