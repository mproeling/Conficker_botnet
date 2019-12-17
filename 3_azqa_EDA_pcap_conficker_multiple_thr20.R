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

output.azqa20 = read.csv("clusters-ctu91-full-20.csv")

mapping.azqa20 = read.csv("mapping2-ctu91-full-20.txt", h = F); colnames(mapping.azqa20) = "clusnum"

#########################################################
# retrieved output file from Azqa Nadeem on 12-08-2019 
setwd("D:\\Documents\\Delft\\Botnets2\\Azqa_delivered\\Fw__poster_Prague")

label20 = read.table("labels-ctu91-full-20.txt")

#########################################################
# retrieved output file from Azqa Nadeem on 12-08-2019 
bytesDist20 = read.matrix("bytesDist-ctu91-full-20.txt")
dportDist20 = read.matrix("dportDist-ctu91-full-20.txt")
gapsDist20 = read.matrix("gapsDist-ctu91-full-20.txt")
sportDist20 = read.matrix("sportDist-ctu91-full-20.txt")

bytesDist20[c(1:10),c(1:10)]
dportDist20[c(1:10),c(1:10)]
gapsDist20[c(1:10),c(1:10)]
sportDist20[c(1:10),c(1:10)]

# conduct PCA 
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/factoextra")
library(factoextra)
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

bytesDist20.pca = prcomp(bytesDist20, center = TRUE, scale. = TRUE)
fviz_eig(bytesDist20.pca)

dportDist20.pca = prcomp(dportDist20, center = TRUE, scale. = TRUE)
fviz_eig(dportDist20.pca)

gapsDist20.pca = prcomp(gapsDist20, center = TRUE, scale. = TRUE)
fviz_eig(gapsDist20.pca)

sportDist20.pca = prcomp(sportDist20, center = TRUE, scale. = TRUE)
fviz_eig(sportDist20.pca)

#############
# Eigenvalues
#eig.val <- get_eigenvalue(bytesDist20.pca)
#eig.val

# Results for Variables
res.var <- get_pca_var(bytesDist20.pca)
#res.var$coord          # Coordinates
#res.var$contrib        # Contributions to the PCs
#res.var$cos2           # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(bytesDist20.pca)
#res.ind$coord          # Coordinates
#res.ind$contrib        # Contributions to the PCs
#res.ind$cos2           # Quality of representation 

# first PCA
plot(res.ind$coord[,1])
test = dist(res.ind$coord[,c(1:2)])

# use PCA output to create covariate matrices
# bytesDist20.pca
res.ind <- get_pca_ind(bytesDist20.pca)
pca1 = as.data.frame(cbind(mapping.azqa20$clusnum, res.ind$coord[,1]))
colnames(pca1) = c("connnum", "bytesDist20.pca1")
pca.output = merge(pca1, output.azqa20, by = "connnum")

# dportDist20.pca
res.ind <- get_pca_ind(dportDist20.pca)
pca1 = as.data.frame(cbind(mapping.azqa20$clusnum, res.ind$coord[,1]))
pca2 = as.data.frame(cbind(mapping.azqa20$clusnum, res.ind$coord[,2]))
pca3 = as.data.frame(cbind(mapping.azqa20$clusnum, res.ind$coord[,3]))
colnames(pca1) = c("connnum", "dportDist20.pca1")
colnames(pca2) = c("connnum", "dportDist20.pca2")
colnames(pca3) = c("connnum", "dportDist20.pca3")
pca.output = merge(pca.output, pca1, by = "connnum")
pca.output = merge(pca.output, pca2, by = "connnum")
pca.output = merge(pca.output, pca3, by = "connnum")

# bytesDist20.pca
res.ind <- get_pca_ind(gapsDist20.pca)
pca1 = as.data.frame(cbind(mapping.azqa20$clusnum, res.ind$coord[,1]))
colnames(pca1) = c("connnum", "gapsDist20.pca1")
pca.output = merge(pca.output, pca1, by = "connnum")

# sportDist20.pca
res.ind <- get_pca_ind(sportDist20.pca)
pca1 = as.data.frame(cbind(mapping.azqa20$clusnum, res.ind$coord[,1]))
pca2 = as.data.frame(cbind(mapping.azqa20$clusnum, res.ind$coord[,2]))
pca3 = as.data.frame(cbind(mapping.azqa20$clusnum, res.ind$coord[,3]))
colnames(pca1) = c("connnum", "sportDist20.pca1")
colnames(pca2) = c("connnum", "sportDist20.pca2")
colnames(pca3) = c("connnum", "sportDist20.pca3")
pca.output = merge(pca.output, pca1, by = "connnum")
pca.output = merge(pca.output, pca2, by = "connnum")
pca.output = merge(pca.output, pca3, by = "connnum")

pca.output = pca.output[pca.output$clusnum != -1,]
output.azqa20 = output.azqa20[output.azqa20$clusnum != -1,]

# now create matrix
azqa.matrix.bern20 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa20[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.bytesDist20 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa20[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.dportDist20pca1 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa20[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.dportDist20pca2 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa20[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.dportDist20pca3 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa20[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.gapsDist20 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa20[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.sportDist20pca1 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa20[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.sportDist20pca2 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa20[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))
azqa.matrix.sportDist20pca3 = as.matrix(get.adjacency(graph.edgelist(as.matrix(unique(output.azqa20[c("srcip", "dstip")])), directed=TRUE), sparse = igraph_opt("sparsematrices")))

for(i in 1:nrow(pca.output)){
  row = which(colnames(azqa.matrix.bern20) == as.character(paste(pca.output[i,]$srcip)))
  col = which(colnames(azqa.matrix.bern20) == as.character(paste(pca.output[i,]$dstip)))
  
  azqa.matrix.bytesDist20[row,col] = pca.output[i,]$bytesDist20.pca1
  azqa.matrix.dportDist20pca1[row,col] = pca.output[i,]$dportDist20.pca1
  azqa.matrix.dportDist20pca2[row,col] = pca.output[i,]$dportDist20.pca2
  azqa.matrix.dportDist20pca3[row,col] = pca.output[i,]$dportDist20.pca3
  azqa.matrix.gapsDist20[row,col] = pca.output[i,]$gapsDist20.pca1
  azqa.matrix.sportDist20pca1[row,col] = pca.output[i,]$sportDist20.pca1
  azqa.matrix.sportDist20pca2[row,col] = pca.output[i,]$sportDist20.pca2
  azqa.matrix.sportDist20pca3[row,col] = pca.output[i,]$sportDist20.pca3
}

##########################################
##########################################
##########################################
### read in the botnet data

## fit blockmodel without covariates
output.combined = BM_bernoulli("SBM",
                               adj = azqa.matrix.bern20,
                               verbosity = 2,
                               autosave = 'SBMconficker_azqa20',
                               plotting = "SBMconficker_azqa20.pdf",
                               exploration_factor = 6,
                               explore_min = 2,
                               explore_max = 10,
                               ncores=4)

output.combined$estimate()
which.max(output.combined$ICL)
Q = 4

## fit blockmodel without covariates
output.combined = BM_bernoulli_covariates_fast("SBM",
                                        adj = azqa.matrix.bern20,
                                        covariates = list(azqa.matrix.bytesDist20,
                                                          azqa.matrix.dportDist20pca1,
                                                          azqa.matrix.dportDist20pca2,
                                                          azqa.matrix.dportDist20pca3,
                                                          azqa.matrix.gapsDist20,
                                                          azqa.matrix.sportDist20pca1,
                                                          azqa.matrix.sportDist20pca2,
                                                          azqa.matrix.sportDist20pca3),
                                        verbosity = 2,
                                        autosave = 'SBMconficker_azqa20_cov_nooutliers',
                                        plotting = "SBMconficker_azqa20_cov_nooutliers.pdf",
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
SBMmemberships = as.data.frame(cbind(row.names(azqa.matrix.bytesDist20), best.model$Z))
SBMmemberships$V2 = as.numeric(paste(SBMmemberships$V2))
SBMmemberships$V3 = as.numeric(paste(SBMmemberships$V3))
SBMmemberships$V4 = as.numeric(paste(SBMmemberships$V4))
SBMmemberships$V5 = as.numeric(paste(SBMmemberships$V5))
SBMmemberships$V6 = as.numeric(paste(SBMmemberships$V6))

SBMmemberships$class = 0
for(i in 1:nrow(SBMmemberships)){
  SBMmemberships[i,]$class = which(max(SBMmemberships[i,c(2:(Q+1))]) == SBMmemberships[i,c(2:(Q+1))])
}

write.table(SBMmemberships, "SBMconficker_azqa20_cov_nooutliers_postprob.csv", sep = ";", quote = F, col.names = T, row.names = F)
SBMmemberships.sel = SBMmemberships[,c(1, ncol(SBMmemberships))]
colnames(SBMmemberships.sel) = c("Id", "class")
write.table(SBMmemberships.sel, paste0("SBMconficker_azqa20_cov_nooutliers_postprob", Q, ".csv"), sep = ";", quote = F, col.names = T, row.names = F)


###########
SBMmemberships5 = read.table("D:\\Documents\\Delft\\Botnets2\\hybridmodel\\output5\\SBMconficker_azqa5_cov_postprob.csv", sep = ";", h = T)
SBMmemberships10 = read.table("D:\\Documents\\Delft\\Botnets2\\hybridmodel\\output10\\SBMconficker_azqa10_cov_postprob.csv", sep = ";", h = T)
SBMmemberships15 = read.table("D:\\Documents\\Delft\\Botnets2\\hybridmodel\\output15\\SBMconficker_azqa15_cov_postprob.csv", sep = ";", h = T)
SBMmemberships20 = read.table("D:\\Documents\\Delft\\Botnets2\\hybridmodel\\output20\\SBMconficker_azqa20_cov_postprob.csv", sep = ";", h = T)



