### SPOT Network Analysis ###
### Figure 5: PIDA Comparison - PREDATION ###
### By: Samantha Gleich ###
### Last Updated: 8/24/23 ###

# Load libraries
library(igraph)
library(ggplot2)
library(ggraph)
library(tidyverse)
library(reshape2)
library(stringr)
library(stringi)

# Find edges that are in surface and DCM networks
surf <- read.csv("Surface_SPOT_Aug2023.csv",header=TRUE,row.names=1)
dcm <- read.csv("DCM_SPOT_Aug2023.csv",header=TRUE,row.names=1)

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
tax <- read.delim("./SPOT/SPOT_2023/taxonomy_90.tsv",header=TRUE,row.names=NULL)
tax$Confidence <- NULL
colnames(tax) <- c("X1","Tax1")
surfEdge <- left_join(surfEdge,tax)
colnames(tax) <- c("X2","Tax2")
surfEdge <- left_join(surfEdge,tax)

# Load in PIDA
pida <- read.csv("./SPOT/SPOT_2023/PIDA_Int.csv",header=TRUE) 
pidaPred <- subset(pida,Taxonomic.interaction=="Prot - Prot" & Interaction=="pred")

# Parasitism edges
surf1 <- colsplit(surfEdge$Tax1,";",c("super1","k1","p1","c1","o1","f1","g1","s1"))
surf2 <- colsplit(surfEdge$Tax2,";",c("super2","k2","p2","c2","o2","f2","g2","s2"))
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
predFull <- predFull %>% distinct(X1,X2) %>% as.data.frame()
predFull <- predFull[c(1:2)]

# Prepare for plotting
outG <-graph_from_data_frame(predFull,directed=FALSE)

outDf <- data.frame(name=c(V(outG)$name))
colnames(outDf) <- "X2"
outDf <- left_join(outDf,tax)
outDf <- outDf %>% arrange(Tax2) %>% as.data.frame()

# Reorder the vertices
mat <- get.adjacency(outG)
mat <- as.matrix(mat)

mat <- mat[order(match(rownames(mat), outDf$X2)), , drop = FALSE]
mat <- mat[, c(outDf$X2)]
outG <- graph_from_adjacency_matrix(mat,mode="undirected")

outDf <- data.frame(name=c(V(outG)$name))
colnames(outDf) <- "X2"
outDf <- left_join(outDf,tax)

# Assign each ASV to taxonomic group for color coding
taxz <- colsplit(outDf$Tax2,";",c("Supergroup","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxz$fin <- ifelse(taxz$Class=="Dinophyceae","Dinoflagellate",NA)
taxz$fin <- ifelse(taxz$Order=="Dino-Group-I","Group I Syndiniales",taxz$fin)
taxz$fin <- ifelse(taxz$Order=="Dino-Group-II","Group II Syndiniales",taxz$fin)
taxz$fin <- ifelse(taxz$Class=="Bacillariophyta","Diatom",taxz$fin)
taxz$fin <- ifelse(grepl("MAST",taxz$Class),"MAST",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Rhizaria","Rhizaria",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Chlorophyta","Chlorophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Cryptophyta","Cryptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Haptophyta","Haptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Ciliophora","Ciliate",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Metazoa","Metazoa",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Stramenopiles" & is.na(taxz$fin),"Other Stramenopiles",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Archaeplastida" & is.na(taxz$fin),"Other Archaeplastids",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Alveolata" & is.na(taxz$fin),"Other Alveolate",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin),"Other Eukaryote",taxz$fin)

# Retain lowest taxonomic assignment of each ASV
outDf$namez <- taxz$fin
V(outG)$fin <- taxz$fin
taxz$Low <- paste(taxz$Genus, taxz$Species,sep="_")
taxz$Low <- ifelse(taxz$Low=="_",taxz$Family,taxz$Low)
low <- colsplit(taxz$Low,"_",c("keep","toss1","toss2"))
V(outG)$namez <- low$keep

# Which edges are in surface and DCM networks?
edgeFin <- data.frame(get.edgelist(outG))
edgeFin$fin1 <- paste(edgeFin$X1,edgeFin$X2,sep="_")
edgeFin$fin2 <- paste(edgeFin$X2,edgeFin$X1,sep="_")
interEdge$fin <- paste(interEdge$X1,interEdge$X2,sep="_")
c <- which(edgeFin$fin1 %in% interEdge$fin)
c2 <- which(edgeFin$fin2 %in% interEdge$fin)
c3 <- c(c,c2)
E(outG)$weight <- 1
E(outG)$weight[E(outG)[c3]] <- -1

# Colors 
taxCols <- c("#E1746D","#76C3D7","#DE5AB1","#D5E0AF","#DED3DC","#87EB58","#D4DC60","#88E5D3","#88AAE1","#DBA85C","#8B7DDA","#9A8D87","#D99CD1","#B649E3","#7EDD90")
names(taxCols) <- c("Chlorophyte","Ciliate","Cryptophyte","Diatom","Haptophyte","Dinoflagellate","MAST","Other Alveolate","Other Archaeplastida","Other Eukaryote","Other Stramenopile","Rhizaria","Group I Syndiniales","Group II Syndiniales","Unknown Eukaryote")

# PLOT IT UP
pdf("Pred_DCM.pdf",width=6,height=6)
ggraph(outG,layout = 'linear', circular = TRUE) + geom_edge_arc(aes(color = as.factor(weight)),alpha=0.8,width=0.7) +
  geom_node_point(shape = 21, size = 2, aes(fill = fin)) +
  theme_graph() +scale_fill_manual(name="Taxonomic Groups",values=c(taxCols))+scale_edge_color_manual(values=c("grey60","indianred"),breaks=c(1,-1))+geom_node_text(aes(label = namez, angle = node_angle(x, y)), hjust = -0.15,size=2,fontface="italic")+
  coord_fixed(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4)) 
dev.off()

# Legend
ggplot(taxz,aes(x=1:38,y=2:39,fill=fin))+geom_point(size=3,shape=21)+scale_fill_manual(name="Taxonomic Group",values=c(taxCols))+theme_classic()+theme(legend.position="bottom")
ggsave("Legend.pdf",width=6,height=4)
