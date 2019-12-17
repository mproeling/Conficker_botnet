rm(list = ls(all = TRUE))

# start with a comparison between normal and infected nodes
normal = c("192.168.1.6","192.168.1.240","192.168.1.36","192.168.1.52","192.168.1.53","192.168.1.55","192.168.1.64","192.168.1.100","192.168.1.155","192.168.1.157")
infected = c("192.168.1.9","192.168.1.71","192.168.1.91","192.168.1.236","192.168.1.238","192.168.1.239","192.168.1.242","192.168.1.243","192.168.1.245","192.168.1.247")

SBMmemberships = read.table("D:\\Documents\\Delft\\Botnets2\\hybridmodel\\output5\\SBMconficker_azqa5_cov_nooutliers_postprob.csv", sep = ";", h = T)
SBMmemberships$label = -1
SBMmemberships[which(SBMmemberships$V1 %in% normal),]$label <- 0
SBMmemberships[which(SBMmemberships$V1 %in% infected),]$label <- 1 

table(SBMmemberships$label, SBMmemberships$class)

output = SBMmemberships[,c("V1", "class", "label")]
colnames(output) = c("Id", "class5", "label5")
write.table(output, "D:\\Documents\\Delft\\Botnets2\\hybridmodel\\output5\\SBMconficker_azqa5_gephi_input_label.csv", sep = ";", quote = F, col.names = T, row.names = F)

######
SBMmemberships = read.table("D:\\Documents\\Delft\\Botnets2\\hybridmodel\\output10\\SBMconficker_azqa10_cov_nooutliers_postprob.csv", sep = ";", h = T)
SBMmemberships$label = -1
SBMmemberships[which(SBMmemberships$V1 %in% normal),]$label <- 0
SBMmemberships[which(SBMmemberships$V1 %in% infected),]$label <- 1 

table(SBMmemberships$label, SBMmemberships$class)

output = SBMmemberships[,c("V1", "class", "label")]
colnames(output) = c("Id", "class10", "label10")
write.table(output, "D:\\Documents\\Delft\\Botnets2\\hybridmodel\\output10\\SBMconficker_azqa10_gephi_input_label.csv", sep = ";", quote = F, col.names = T, row.names = F)

######
SBMmemberships = read.table("D:\\Documents\\Delft\\Botnets2\\hybridmodel\\output15\\SBMconficker_azqa15_cov_nooutliers_postprob.csv", sep = ";", h = T)
SBMmemberships$label = -1
SBMmemberships[which(SBMmemberships$V1 %in% normal),]$label <- 0
SBMmemberships[which(SBMmemberships$V1 %in% infected),]$label <- 1 

table(SBMmemberships$label, SBMmemberships$class)

output = SBMmemberships[,c("V1", "class", "label")]
colnames(output) = c("Id", "class15", "label15")
write.table(output, "D:\\Documents\\Delft\\Botnets2\\hybridmodel\\output15\\SBMconficker_azqa15_gephi_input_label.csv", sep = ";", quote = F, col.names = T, row.names = F)

######
SBMmemberships = read.table("D:\\Documents\\Delft\\Botnets2\\hybridmodel\\output20\\SBMconficker_azqa20_cov_nooutliers_postprob.csv", sep = ";", h = T)
SBMmemberships$label = -1
SBMmemberships[which(SBMmemberships$V1 %in% normal),]$label <- 0
SBMmemberships[which(SBMmemberships$V1 %in% infected),]$label <- 1 

table(SBMmemberships$label, SBMmemberships$class)

output = SBMmemberships[,c("V1", "class", "label")]
colnames(output) = c("Id", "class20", "label20")
write.table(output, "D:\\Documents\\Delft\\Botnets2\\hybridmodel\\output20\\SBMconficker_azqa20_gephi_input_label.csv", sep = ";", quote = F, col.names = T, row.names = F)
######