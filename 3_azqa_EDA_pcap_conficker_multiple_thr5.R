###########################
## MODEL WITH COVARIATES ##
###########################
setwd("D:\\Documents\\Delft\\Botnets2\\Azqa_delivered\\Re__poster_Prague_clusters_mapping/")

###########################
rm(list = ls(all = TRUE))

library(GGally)
library(network)
library(sna)
library(ggplot2)
library(igraph)
library(blockmodels)
library(Matrix)
library(parallel)
library(tseries)

# read in the data with the time window
botnet.data = read.table("D:\\Documents\\Delft\\Botnets2\\capture.botnet2.infected.1.output.csv", sep = "\t", h = T)
colnames(botnet.data) = c("Source", "Target", "TotalObs", "Protocol", "MeanLength", "SDLength", "MeanTimeDif")

output.azqa5 = read.csv("clusters-ctu91-full-5.csv")

colnames(output.azqa5)[c(4:5)] = c("Source", "Target")
write.table(output.azqa5, "output.azqa5.csv", sep = ";", quote = F, col.names = T, row.names = F)

mapping.azqa5 = read.csv("mapping2-ctu91-full-5.txt", h = F); colnames(mapping.azqa5) = "clusnum"

#########################################################
# retrieved output file from Azqa Nadeem on 12-08-2019 
setwd("D:\\Documents\\Delft\\Botnets2\\Azqa_delivered\\Fw__poster_Prague")

label5 = read.table("labels-ctu91-full-5.txt")

#########################################################
# retrieved output file from Azqa Nadeem on 12-08-2019 
bytesDist5 = read.matrix("bytesDist-ctu91-full-5.txt")
dportDist5 = read.matrix("dportDist-ctu91-full-5.txt")
gapsDist5 = read.matrix("gapsDist-ctu91-full-5.txt")
sportDist5 = read.matrix("sportDist-ctu91-full-5.txt")

bytesDist5[c(1:10),c(1:10)]
dportDist5[c(1:10),c(1:10)]
gapsDist5[c(1:10),c(1:10)]
sportDist5[c(1:10),c(1:10)]

# conduct PCA 
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/factoextra")
library(factoextra)
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

bytesDist5.pca = prcomp(bytesDist5, center = TRUE, scale. = TRUE)
fviz_eig(bytesDist5.pca)

dportDist5.pca = prcomp(dportDist5, center = TRUE, scale. = TRUE)
fviz_eig(dportDist5.pca)

gapsDist5.pca = prcomp(gapsDist5, center = TRUE, scale. = TRUE)
fviz_eig(gapsDist5.pca)

sportDist5.pca = prcomp(sportDist5, center = TRUE, scale. = TRUE)
fviz_eig(sportDist5.pca)

#############
# Eigenvalues
#eig.val <- get_eigenvalue(bytesDist5.pca)
#eig.val

# Results for Variables
res.var <- get_pca_var(bytesDist5.pca)
#res.var$coord          # Coordinates
#res.var$contrib        # Contributions to the PCs
#res.var$cos2           # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(bytesDist5.pca)
#res.ind$coord          # Coordinates
#res.ind$contrib        # Contributions to the PCs
#res.ind$cos2           # Quality of representation 

# first PCA
plot(res.ind$coord[,1])
test = dist(res.ind$coord[,c(1:2)])

# use PCA output to create covariate matrices
# bytesDist5.pca
res.ind <- get_pca_ind(bytesDist5.pca)
pca1 = as.data.frame(cbind(mapping.azqa5$clusnum, res.ind$coord[,1]))
colnames(pca1) = c("connnum", "bytesDist5.pca1")
pca.output = merge(pca1, output.azqa5, by = "connnum")

# dportDist5.pca
res.ind <- get_pca_ind(dportDist5.pca)
pca1 = as.data.frame(cbind(mapping.azqa5$clusnum, res.ind$coord[,1]))
pca2 = as.data.frame(cbind(mapping.azqa5$clusnum, res.ind$coord[,2]))
pca3 = as.data.frame(cbind(mapping.azqa5$clusnum, res.ind$coord[,3]))
colnames(pca1) = c("connnum", "dportDist5.pca1")
colnames(pca2) = c("connnum", "dportDist5.pca2")
colnames(pca3) = c("connnum", "dportDist5.pca3")
pca.output = merge(pca.output, pca1, by = "connnum")
pca.output = merge(pca.output, pca2, by = "connnum")
pca.output = merge(pca.output, pca3, by = "connnum")

# bytesDist5.pca
res.ind <- get_pca_ind(gapsDist5.pca)
pca1 = as.data.frame(cbind(mapping.azqa5$clusnum, res.ind$coord[,1]))
colnames(pca1) = c("connnum", "gapsDist5.pca1")
pca.output = merge(pca.output, pca1, by = "connnum")

# sportDist5.pca
res.ind <- get_pca_ind(sportDist5.pca)
pca1 = as.data.frame(cbind(mapping.azqa5$clusnum, res.ind$coord[,1]))
pca2 = as.data.frame(cbind(mapping.azqa5$clusnum, res.ind$coord[,2]))
pca3 = as.data.frame(cbind(mapping.azqa5$clusnum, res.ind$coord[,3]))
colnames(pca1) = c("connnum", "sportDist5.pca1")
colnames(pca2) = c("connnum", "sportDist5.pca2")
colnames(pca3) = c("connnum", "sportDist5.pca3")
pca.output = merge(pca.output, pca1, by = "connnum")
pca.output = merge(pca.output, pca2, by = "connnum")
pca.output = merge(pca.output, pca3, by = "connnum")

# remove with label -1
pca.output = pca.output[pca.output$clusnum != -1,]
output.azqa5 = output.azqa5[output.azqa5$clusnum != -1,]

# now create matrix
azqa.matrix.bern5 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa5[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.bytesDist5 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa5[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.dportDist5pca1 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa5[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.dportDist5pca2 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa5[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.dportDist5pca3 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa5[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.gapsDist5 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa5[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.sportDist5pca1 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa5[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.sportDist5pca2 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa5[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.sportDist5pca3 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa5[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.clusnum5 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa5[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))

for(i in 1:nrow(pca.output)){
  row = which(colnames(azqa.matrix.bern5) == as.character(paste(pca.output[i,]$srcip)))
  col = which(colnames(azqa.matrix.bern5) == as.character(paste(pca.output[i,]$dstip)))
  
  azqa.matrix.bytesDist5[row,col] = pca.output[i,]$bytesDist5.pca1
  azqa.matrix.dportDist5pca1[row,col] = pca.output[i,]$dportDist5.pca1
  azqa.matrix.dportDist5pca2[row,col] = pca.output[i,]$dportDist5.pca2
  azqa.matrix.dportDist5pca3[row,col] = pca.output[i,]$dportDist5.pca3
  azqa.matrix.gapsDist5[row,col] = pca.output[i,]$gapsDist5.pca1
  azqa.matrix.sportDist5pca1[row,col] = pca.output[i,]$sportDist5.pca1
  azqa.matrix.sportDist5pca2[row,col] = pca.output[i,]$sportDist5.pca2
  azqa.matrix.sportDist5pca3[row,col] = pca.output[i,]$sportDist5.pca3
  
  azqa.matrix.clusnum5[row,col] = pca.output[i,]$clusnum
}

##########################################
##########################################
##########################################
### read in the botnet data

## fit blockmodel without covariates
output.combined = BM_bernoulli("SBM",
                               adj = azqa.matrix.bern5,
                               verbosity = 2,
                               autosave = 'SBMconficker_azqa5',
                               plotting = "SBMconficker_azqa5.pdf",
                               exploration_factor = 6,
                               explore_min = 2,
                               explore_max = 10,
                               ncores=4)

output.combined$estimate()
which.max(output.combined$ICL)
Q = 4

## fit blockmodel without covariates
output.combined = BM_bernoulli_covariates_fast("SBM",
                                        adj = azqa.matrix.bern5,
                                        covariates = list(azqa.matrix.bytesDist5,
                                                          azqa.matrix.dportDist5pca1,
                                                          azqa.matrix.dportDist5pca2,
                                                          azqa.matrix.dportDist5pca3,
                                                          azqa.matrix.gapsDist5,
                                                          azqa.matrix.sportDist5pca1,
                                                          azqa.matrix.sportDist5pca2,
                                                          azqa.matrix.sportDist5pca3),
                                        verbosity = 2,
                                        autosave = 'SBMconficker_azqa5_cov_nooutliers',
                                        plotting = "SBMconficker_azqa5_cov_nooutliers.pdf",
                                        exploration_factor = 6,
                                        explore_min = 2,
                                        explore_max = 10,
                                        ncores=4)

output.combined$estimate()
which.max(output.combined$ICL)
Q = 4

# get best model
best.model = output.combined$memberships[[which.max(output.combined$ICL)]]
best.model = output.combined$memberships[[Q]]
best.model$plot()
output.combined$plot_parameters(Q)
# memberships 
SBMmemberships = as.data.frame(cbind(row.names(azqa.matrix.bytesDist5), best.model$Z))
SBMmemberships$V2 = as.numeric(paste(SBMmemberships$V2))
SBMmemberships$V3 = as.numeric(paste(SBMmemberships$V3))
SBMmemberships$V4 = as.numeric(paste(SBMmemberships$V4))
SBMmemberships$V5 = as.numeric(paste(SBMmemberships$V5))

SBMmemberships$class = 0
for(i in 1:nrow(SBMmemberships)){
  SBMmemberships[i,]$class = which(max(SBMmemberships[i,c(2:(Q+1))]) == SBMmemberships[i,c(2:(Q+1))])
}

write.table(SBMmemberships, "SBMconficker_azqa5_cov_nooutliers_postprob.csv", sep = ";", quote = F, col.names = T, row.names = F)
SBMmemberships.sel = SBMmemberships[,c(1, ncol(SBMmemberships))]
colnames(SBMmemberships.sel) = c("Id", "class")
write.table(SBMmemberships.sel, paste0("SBMconficker_azqa5_cov_nooutliers_postprob", Q, ".csv"), sep = ";", quote = F, col.names = T, row.names = F)

####
gephi_input = unique(output.azqa5[,c("srcip", "dstip")])
colnames(gephi_input) = c("Source", "Target")
write.table(gephi_input, paste0("SBMconficker_azqa5_gephi_input.csv"), sep = ";", quote = F, col.names = T, row.names = F)
