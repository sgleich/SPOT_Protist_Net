### SPOT Network Analysis ###
### Figure 5: PIDA Comparison - PARASITISM ###
### By: Samantha Gleich ###
### Last Updated: 2/13/24 ###

# Load libraries
library(igraph)
library(ggplot2)
library(ggraph)
library(tidyverse)
library(reshape2)
library(stringr)
library(stringi)

# Find edges that are in surface and DCM networks
surf <- read.csv("Surf_SPOT_Feb12_2024.csv",header=TRUE,row.names=1)
dcm <- read.csv("DCM_SPOT_Feb12_2024.csv",header=TRUE,row.names=1)

namez <- rownames(surf)
namez <- str_remove_all(namez,"S_")
rownames(surf) <- namez
colnames(surf) <- namez
surf <- as.matrix(surf)
surfG <- graph_from_adjacency_matrix(surf,mode="undirected")

namez <- rownames(dcm)
namez <- str_remove_all(namez,"S_")
rownames(dcm) <- namez
colnames(dcm) <- namez
dcm <- as.matrix(dcm)
dcmG <- graph_from_adjacency_matrix(dcm,mode="undirected")

intG <- graph.intersection(surfG,dcmG)
interEdge <- data.frame(get.edgelist(intG))

# Surface graph
surfEdge <- data.frame(get.edgelist(dcmG))
tax <- read.delim("taxonomy.tsv",header=TRUE,row.names=NULL)
tax$Confidence <- NULL
colnames(tax) <- c("X1","Tax1")
surfEdge <- left_join(surfEdge,tax)
colnames(tax) <- c("X2","Tax2")
surfEdge <- left_join(surfEdge,tax)

# Load in PIDA
pida <- read.csv("PIDA_Int.csv",header=TRUE) 
pidaPar <- subset(pida,Taxonomic.interaction=="Prot - Prot" & Interaction=="par")

# Parasitism edges
surf1 <- colsplit(surfEdge$Tax1,";",c("super1","d1","k1","p1","c1","o1","f1","g1","s1"))
surf2 <- colsplit(surfEdge$Tax2,";",c("super2","d1","k2","p2","c2","o2","f2","g2","s2"))
surfEdge <- cbind(surfEdge,surf1,surf2)

# Find all unique genera in surface network
genera1 <- unique(surfEdge$g1)
genera2 <- unique(surfEdge$g2)
generaUnique <- c(genera1,genera2)
generaUnique <- unique(generaUnique)

# Find all unique genera in PIDA parasitism database
pidaGenera1 <- unique(pidaPar$Genus.org1)
pidaGenera2 <- unique(pidaPar$Genus.org2)
pidaGeneraUnique <- c(pidaGenera1,pidaGenera2)
pidaGeneraUnique <- unique(pidaGeneraUnique)

# Find genera that are in surface network and PIDA parasitism database
generaUniqueSub <- generaUnique[generaUnique %in% pidaGeneraUnique]

# The separate out Group I Syndiniales, Group II Syndiniales, and Cercozoans from surface network
synG1 <- subset(surfEdge,grepl("Group-I-",surfEdge$Tax1)|grepl("Group-I-",surfEdge$Tax2))
synG2 <- subset(surfEdge,grepl("Group-II-",surfEdge$Tax1)|grepl("Group-II-",surfEdge$Tax2))
cerc<- subset(surfEdge,grepl("Cercoz",surfEdge$Tax1)|grepl("Cercoz",surfEdge$Tax2))

# Separate PIDA relationships into Syndiniales Group I, Syndiniales Group II, and Cercozoans
pidaCerc <- subset(pidaPar,Taxonomic.level.3..org1=="Cercozoa"|Taxonomic.level.3..org2=="Cercozoa")
pidaCerc <- subset(pidaCerc,Taxonomic.level.3..org1=="Diatomea"|Taxonomic.level.3..org2=="Diatomea")
pidaGroup1 <- subset(pidaPar,Genus.org1=="Syndiniales_Group_I"|Genus.org2=="Syndiniales_Group_I"|Genus.org1=="Euduboscquella"|Genus.org2=="Euduboscquella")
pidaGroup2 <- subset(pidaPar,Genus.org1=="Syndiniales_Group_II"|Genus.org2=="Syndiniales_Group_II"|Genus.org1=="Amoebophrya"|Genus.org2=="Amoebophrya"|Genus.org1=="Syndinium"|Genus.org2=="Syndinium")

# Find genera that are parasitic hosts in PIDA and exist in surface network
synG1Out <- subset(synG1, g1 %in% pidaGroup1$Genus.org1 |g2 %in% pidaGroup1$Genus.org1)
synG2Out <- subset(synG2, g1 %in% pidaGroup2$Genus.org1 | g2 %in% pidaGroup2$Genus.org1 )
cercOut <- subset(cerc, g1 %in% pidaCerc$Genus.org1 | g2 %in% pidaCerc$Genus.org1)
cercOut <- subset(cercOut, grepl("Cryothecomonas", cercOut$f1)|grepl("Cryothecomonas", cercOut$f2))

# Combine and reorder edge table
#full <- synG2Out
#full <- rbind(synG1Out,cercOut)
full <- rbind(synG1Out,synG2Out,cercOut)

full <- full[c(1:4)]
taxz <- colsplit(full$Tax1,";",c("Supergroup","Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxz$fin <- ifelse(taxz$Class=="Dinophyceae","Dinoflagellate",NA)
taxz$fin <- ifelse(taxz$Order=="Dino-Group-I","Group I Syndiniales",taxz$fin)
taxz$fin <- ifelse(taxz$Order=="Dino-Group-II","Group II Syndiniales",taxz$fin)
taxz$fin <- ifelse(taxz$Class=="Mediophyceae"|taxz$Class=="Coscinodiscophyceae"|taxz$Class=="Bacillariophyceae","Diatom",taxz$fin)
taxz$fin <- ifelse(grepl("MAST",taxz$Family),"MAST",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Rhizaria","Rhizaria",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Chlorophyta","Chlorophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Cryptophyta","Cryptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Haptophyta","Haptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Ciliophora","Ciliate",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Metazoa","Metazoa",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Stramenopiles" & is.na(taxz$fin),"Other Stramenopiles",taxz$fin)
taxz$fin <- ifelse(taxz$Domain=="Archaeplastida" & is.na(taxz$fin),"Other Archaeplastids",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Alveolata" & is.na(taxz$fin),"Other Alveolate",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin),"Other Eukaryote",taxz$fin)

full$fin1 <- taxz$fin

taxz <- colsplit(full$Tax2,";",c("Supergroup","Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxz$fin <- ifelse(taxz$Class=="Dinophyceae","Dinoflagellate",NA)
taxz$fin <- ifelse(taxz$Order=="Dino-Group-I","Group I Syndiniales",taxz$fin)
taxz$fin <- ifelse(taxz$Order=="Dino-Group-II","Group II Syndiniales",taxz$fin)
taxz$fin <- ifelse(taxz$Class=="Mediophyceae"|taxz$Class=="Coscinodiscophyceae"|taxz$Class=="Bacillariophyceae","Diatom",taxz$fin)
taxz$fin <- ifelse(grepl("MAST",taxz$Family),"MAST",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Rhizaria","Rhizaria",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Chlorophyta","Chlorophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Cryptophyta","Cryptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Haptophyta","Haptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Ciliophora","Ciliate",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Metazoa","Metazoa",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Stramenopiles" & is.na(taxz$fin),"Other Stramenopiles",taxz$fin)
taxz$fin <- ifelse(taxz$Domain=="Archaeplastida" & is.na(taxz$fin),"Other Archaeplastids",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Alveolata" & is.na(taxz$fin),"Other Alveolate",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin),"Other Eukaryote",taxz$fin)
full$fin2 <- taxz$fin

cPlot <- full %>% select("fin1","fin2") %>% group_by(fin1,fin2) %>% tally()


taxCols <- c("#E1746D","#76C3D7","#DE5AB1","#D5E0AF","#DED3DC","#87EB58","#D4DC60","#88E5D3","#88AAE1","#DBA85C","#8B7DDA","#9A8D87","#D99CD1","#B649E3","#7EDD90")
names(taxCols) <- c("Chlorophyte","Ciliate","Cryptophyte","Diatom","Haptophyte","Dinoflagellate","MAST","Other Alveolate","Other Archaeplastida","Other Eukaryote","Other Stramenopile","Rhizaria","Group I Syndiniales","Group II Syndiniales","Unknown Eukaryote")

pdf("../DCM_Parasitism.pdf")
chordDiagram(cPlot,annotationTrack =c("grid"),annotationTrackHeight = c(0.05, 0.3),grid.col = c(taxCols),transparency = 0.6)
dev.off()


full$num1 <- 1:nrow(full)
full$num2 <- 2:(nrow(full)+1)
s <- full[1,]
s$fin1 <- "Rhizaria"
full <- rbind(full,s)

ggplot(full,aes(num1,num2,fill=fin1))+geom_point(shape=22,size=4)+theme_classic()+scale_fill_manual(name="Taxonomic Group",values=c(taxCols))
ggsave("../legend.pdf")


