rm(list = ls(all=TRUE))
######################################################
######################################################
##                                                  ##
##    PCAP analyses CTU MALWARE                     ##
##    Mark Patrick Roeling, July 2019               ##                       
##                                                  ##
######################################################
######################################################

#install.packages("GGally")
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(igraph)
library(blockmodels)
library(Matrix)
library(parallel)

#setwd("D:/Documents/oxford/Botnets2/")
setwd("D:/Documents/oxford/Botnets2/")

# downloaded from https://mcfp.felk.cvut.cz/publicDatasets/CTU-Malware-Capture-Botnet-91/
# pcap opened in wireshark and exported as .csv

dataset = read.csv("capture.botnet2.infected.1.csv")

# there are 3 isolated clusters in this dataset so we remove those clusters
sourceSubset = dataset[dataset$Source == "192.168.0.118",]
DestinationSubset = dataset[dataset$Destination == "192.168.0.118",] # empty
write.table(sourceSubset, "excluded.network.192.168.0.118", sep = "\t", quote = F, col.names = T, row.names = F)
dataset.clean1 = dataset[dataset$Source != "192.168.0.118",]

sourceSubset = dataset[dataset$Source == "Broadcast",] # empty
DestinationSubset = dataset[dataset$Destination == "Broadcast",]
remove.set = as.character(unique(DestinationSubset$Source))
write.table(DestinationSubset, "excluded.network.Broadcast", sep = "\t", quote = F, col.names = T, row.names = F)
dataset.clean2 = dataset.clean1[dataset.clean1$Destination != "Broadcast" &
                                dataset.clean1$Destination != "CDP/VTP/DTP/PAgP/UDLD" &
                                dataset.clean1$Destination != "Asiarock_64:7b:1a" &
                                dataset.clean1$Destination != "Routerbo_33:45:5e" &
                                dataset.clean1$Destination != "Dell_2e:d0:6e" &
                                dataset.clean1$Destination != "Asiarock_64:7a:bc" &
                                dataset.clean1$Destination != "Asiarock_64:7b:22" &
                                dataset.clean1$Destination != "QuantaCo_0c:79:86" &
                                dataset.clean1$Destination != "Asiarock_64:7a:f2" &
                                dataset.clean1$Destination != "Asiarock_64:7a:80" &
                                dataset.clean1$Destination != "SoyoComp_a0:ad:bb" &
                                dataset.clean1$Destination != "EncoreNe_ec:1d:7f" &
                                dataset.clean1$Destination != "Asiarock_64:7b:43" &
                                dataset.clean1$Destination != "Asiarock_64:7b:41" &
                                dataset.clean1$Destination != "Asiarock_64:7a:88" &
                                dataset.clean1$Destination != "Asiarock_64:7b:47" &
                                dataset.clean1$Destination != "Asiarock_64:7b:30" &
                                dataset.clean1$Destination != "Asiarock_64:7b:0c" &
                                dataset.clean1$Destination != "Netronix_0e:ed:3a" &
                                dataset.clean1$Destination != "Netronix_0e:ef:c4" &
                                dataset.clean1$Source != "Asiarock_64:7b:1a" &
                                dataset.clean1$Source != "Routerbo_33:45:5e" &
                                dataset.clean1$Source != "Dell_2e:d0:6e" &
                                dataset.clean1$Source != "Asiarock_64:7a:bc" &
                                dataset.clean1$Source != "Asiarock_64:7b:22" &
                                dataset.clean1$Source != "QuantaCo_0c:79:86" &
                                dataset.clean1$Source != "Asiarock_64:7a:f2" &
                                dataset.clean1$Source != "Asiarock_64:7a:80" &
                                dataset.clean1$Source != "SoyoComp_a0:ad:bb" &
                                dataset.clean1$Source != "EncoreNe_ec:1d:7f" &
                                dataset.clean1$Source != "Asiarock_64:7b:43" &
                                dataset.clean1$Source != "Asiarock_64:7b:41" &
                                dataset.clean1$Source != "Asiarock_64:7a:88" &
                                dataset.clean1$Source != "Asiarock_64:7b:47" &
                                dataset.clean1$Source != "Asiarock_64:7b:30" &
                                dataset.clean1$Source != "Asiarock_64:7b:0c" &
                                dataset.clean1$Source != "Netronix_0e:ed:3a" &
                                dataset.clean1$Source != "Netronix_0e:ef:c4",]

DestinationSubset = dataset[dataset.clean2$Destination == "ff02::16" | 
                            dataset.clean2$Destination == "ff02::2" | 
                            dataset.clean2$Destination != "ff02::1:ff2e:d06e",]
write.table(DestinationSubset, "excluded.network.other", sep = "\t", quote = F, col.names = T, row.names = F)
dataset.clean3 = dataset.clean2[dataset.clean2$Destination != "ff02::16" & 
                                dataset.clean2$Destination != "ff02::2" & 
                                dataset.clean2$Destination != "ff02::1:ff2e:d06e" ,]

input = dataset.clean3[,c(3,4)]
colnames(input) = c("Source", "Target")
write.table(input, "capture.botnet2.infected.1.cleaninput.csv", sep = ";", quote = F, col.names = T, row.names = F)

# the intention of this loop is to create a file where there are no duplicate rows 
# so we obtain unique [Source to Destination] rows. But we dont want to lose information 
# about bytes (length) or protocol, so per Source-Destination pair we sum the length per protocol.
# the idea is that in the graph (visualisation) we get an edge per protocol and use length
# as weight (thick line = more traffic over that protocol)
dataset = dataset.clean3 ### OVERWRITE

uniqueIP = length(unique(dataset$Source))
for (i in 1:uniqueIP){
  sourceIP = unique(dataset$Source)[i]
  sourceSubset = dataset[dataset$Source == sourceIP,]
  # calculate unique Destination IP
  uniqueIPdest =  length(unique(sourceSubset$Destination))
  for (j in 1:uniqueIPdest){
    destinationIP = unique(sourceSubset$Destination)[j]
    destinationSubset = sourceSubset[sourceSubset$Destination == destinationIP,]
    # calculate unique Destination IP 
    uniqueIPprot = length(unique(destinationSubset$Protocol))
    for(k in 1:uniqueIPprot){
      Protocol = unique(destinationSubset$Protocol)[k]
      finalSet = destinationSubset[destinationSubset$Protocol == Protocol,]

      # calculate number of these observations 
      Tobs = nrow(finalSet)
      PacketMean = mean(finalSet$Length)
      PacketSd = sd(finalSet$Length)
      MeanTimeDiff = median(diff(finalSet$Time))
      
      output = t(as.data.frame(c(paste(sourceIP), paste(destinationIP), Tobs, paste(Protocol), PacketMean, PacketSd, MeanTimeDiff)))
      
      # write data
      if(i == 1 & j == 1 & k == 1){
        write.table(output, "capture.botnet2.infected.1.output.csv", quote = F, sep = "\t", row.names = F, col.names = T)
      } else {
        write.table(output, "capture.botnet2.infected.1.output.csv", quote = F, sep = "\t", row.names = F, col.names = F, append = T)
      }
    }
  }
  print(i)
}