### SPOT Network Analysis ###
### Figure 6: PIDA Comparison - PREDATION ###
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
surfEdge <- data.frame(get.edgelist(surfG))
tax <- read.delim("taxonomy.tsv",header=TRUE,row.names=NULL)
tax$Confidence <- NULL
colnames(tax) <- c("X1","Tax1")
surfEdge <- left_join(surfEdge,tax)
colnames(tax) <- c("X2","Tax2")
surfEdge <- left_join(surfEdge,tax)

surfEdge <-subset(surfEdge,grepl("Eukaryota",surfEdge$Tax1)|grepl("Eukaryota",surfEdge$Tax2))

# Load in PIDA
pida <- read.csv("PIDA_Int.csv",header=TRUE) 
pidaPred <- subset(pida,Taxonomic.interaction=="Prot - Prot" & Interaction=="pred")

# Parasitism edges
surf1 <- colsplit(surfEdge$Tax1,";",c("super1","d1","k1","p1","c1","o1","f1","g1","s1"))
surf2 <- colsplit(surfEdge$Tax2,";",c("super2","d2","k2","p2","c2","o2","f2","g2","s2"))
surfEdge <- cbind(surfEdge,surf1,surf2)

# Find all unique genera in surface network
genera1 <- unique(surfEdge$g1)
genera2 <- unique(surfEdge$g2)
generaUnique <- c(genera1,genera2)
generaUnique <- unique(generaUnique)

# Find all unique genera in PIDA parasitism database
pidaGenera1 <- unique(pidaPred$Genus.org1)
pidaGenera2 <- unique(pidaPred$Genus.org2)
pidaGeneraUnique <- c(pidaGenera1,pidaGenera2)
pidaGeneraUnique <- unique(pidaGeneraUnique)

# Find genera that are in surface network and PIDA parasitism database
generaUniqueSub <- generaUnique[generaUnique %in% pidaGeneraUnique]
assoc <- data.frame(a1=paste(pidaPred$Genus.org1,pidaPred$Genus.org2,sep="_"))
assoc$a2 <- paste(pidaPred$Genus.org2,pidaPred$Genus.org1,sep="_")

# Subset dataframe to include only those associations that are in PIDA
surfEdge <- subset(surfEdge,g1 %in% generaUniqueSub & g2 %in% generaUniqueSub)
surfEdge$fin <- paste(surfEdge$g1,surfEdge$g2,sep="_")
surfEdge1 <- subset(surfEdge,fin %in% assoc$a1)
surfEdge2 <- subset(surfEdge,fin %in% assoc$a2)
predFull <- rbind(surfEdge1,surfEdge2)
predFull <- predFull %>% distinct(X1,X2,.keep_all = TRUE) %>% as.data.frame()

full <- predFull[c(1:4)]
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

pdf("../Surf_Pred.pdf")
chordDiagram(cPlot,annotationTrack =c("grid"),annotationTrackHeight = c(0.05, 0.3),grid.col = c(taxCols),transparency = 0.6)
dev.off()


full$num1 <- 1:nrow(full)
full$num2 <- 2:(nrow(full)+1)
#s <- full[1,]
#s$fin1 <- "Rhizaria"
#full <- rbind(full,s)

ggplot(full,aes(num1,num2,fill=fin2))+geom_point(shape=22,size=4)+theme_classic()+scale_fill_manual(name="Taxonomic Group",values=c(taxCols))
ggsave("../legendpred.pdf")



