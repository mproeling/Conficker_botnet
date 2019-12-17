#Build gephi input file
#Index of /publicDatasets/CTU-Malware-Capture-Botnet-91
# dataset = read.csv("C:\\Users\\Gebruiker\\Downloads\\capture.botnet2.infected.1.csv")
dataset = read.csv("D:/Documents/Delft/Botnets2/capture.botnet2.infected.1.csv")
udataset = unique(dataset[,c(3,4)])
colnames(udataset) = c("Source", "Target")
write.table(udataset, "C:\\Users\\Gebruiker\\Downloads\\capture.botnet2.infected.1_input.csv", quote = F, sep = ";", col.names = T, row.names = F)

