### SPOT Network Analysis ###
### Figure 7: PIDA Comparison - SYMBIOSIS ###
### By: Samantha Gleich ###
### Last Updated: 2/6/24 ###

# Load libraries
library(igraph)
library(ggplot2)
library(ggraph)
library(tidyverse)
library(reshape2)
library(stringr)
library(stringi)

# Find edges that are in surface and DCM networks
surf <- read.csv("Surf_SPOT_Feb2024.csv",header=TRUE,row.names=1)
dcm <- read.csv("DCM_SPOT_Feb2024.csv",header=TRUE,row.names=1)

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
pidaSymb <- subset(pida,Taxonomic.interaction=="Prot - Prot" & Interaction=="symb")

# Symbiosis edges
surf1 <- colsplit(surfEdge$Tax1,";",c("super1","d1","k1","p1","c1","o1","f1","g1","s1"))
surf2 <- colsplit(surfEdge$Tax2,";",c("super2","d2","k2","p2","c2","o2","f2","g2","s2"))
surfEdge <- cbind(surfEdge,surf1,surf2)

subs <- subset(pidaSymb,Genus.org1 %in% surfEdge$g1 | Genus.org2 %in% surfEdge$g1 |Genus.org1 %in% surfEdge$g2 | Genus.org2 %in% surfEdge$g2)

# Subset specific host groups
poly <- subset(surfEdge,grepl("Polycystinea",surfEdge$Tax1)|grepl("Polycystinea",surfEdge$Tax2))
acanth <- subset(surfEdge,grepl("Acantharea",surfEdge$Tax1)|grepl("Acantharea",surfEdge$Tax2))
dino <- subset(surfEdge,grepl("Dinophyceae",surfEdge$Tax1)|grepl("Dinophyceae",surfEdge$Tax2))
cil <- subset(surfEdge,grepl("Ciliophor",surfEdge$Tax1)|grepl("Ciliophor",surfEdge$Tax2))
diatom <- subset(surfEdge,grepl("Mediophyceae",surfEdge$Tax1)|grepl("Mediophyceae",surfEdge$Tax2)|grepl("Coscinodiscophyceae",surfEdge$Tax1)|grepl("Coscinodiscophyceae",surfEdge$Tax2)|grepl("Bacillariophyceae",surfEdge$Tax1)|grepl("Bacillariophyceae",surfEdge$Tax2))
dict <- subset(surfEdge,grepl("Dictyochophyce",surfEdge$Tax1)|grepl("Dictyochophyce",surfEdge$Tax2))

# See if symbionts are in host groups
p <- c("Gyrodinium","Amphidinium","Scrippsiella","Brandtodinium","Gymnodinium","Gymnoxanthella","Caryotoma")
poly <- subset(poly,g1 %in% p|g2 %in% p)

a<- c("Brandtodinium","Phaeocystis","Chrysochromulina","Pelagodinium","Azadinium","Scrippsiella","Heterocapsa")
acanth <- subset(acanth,g1 %in% a|g2 %in% a)
acanth <- subset(acanth,grepl("_F",acanth$o1)|grepl("_F",acanth$o2))

ci <- c("Mychonastes","Scenedesmus","Chlorella","Micractinium","Symbiodinium","Choricystis","Meyerella","Coccomyxa","Chaetoceros","Fragilariopsis","Navicula","Synedra")
cil <- subset(cil, g1 %in% ci|g2 %in% ci)
cil <- subset(cil, grepl("Eutintin",cil$g1)|grepl("Eutintin",cil$g2))

diatom1 <- subset(diatom,grepl("Eutinn",diatom$g1)|grepl("Eutinn",diatom$g2))
diatom2 <- subset(diatom,grepl("MAST-3",diatom$g1)|grepl("MAST-3",diatom$g2))
diatom3 <- subset(diatom,grepl("MAST-4",diatom$g1)|grepl("MAST-4",diatom$g2))
diatom2 <- subset(diatom2,grepl("Leptocylindrus",diatom2$g1)|grepl("Leptocylindrus",diatom2$g2))
diatom3 <- subset(diatom3,grepl("Leptocylindrus",diatom3$g1)|grepl("Leptocylindrus",diatom3$g2))

# Combine all
full <- rbind(poly,acanth,cil,diatom1,diatom2,diatom3)

# Prepare for plotting
outG <-graph_from_data_frame(full,directed=FALSE)
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
taxz <- colsplit(outDf$Tax2,";",c("Supergroup","Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
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
unique(taxz$fin)

# Retain lowest taxonomic assignment of each ASV
V(outG)$fin <- taxz$fin
taxz$Low <- paste(taxz$Genus, taxz$Species,sep="_")
taxz$Low <- ifelse(taxz$Low=="_",taxz$Order,taxz$Low)
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

d <- data.frame(Name=c(V(outG)$namez))
namez <- read.csv("PidaSymbLetters.csv",header=TRUE)
d <- left_join(d,namez)
V(outG)$namez <- d$Letter


# PLOT IT UP
pdf("../DCM_Other.pdf",width=6,height=4)
ggraph(outG, layout = 'linear', circular = TRUE) + geom_edge_arc(aes(color = as.factor(weight)),alpha=0.95,width=0.7) +
  geom_node_point(shape = 21, size = 6, aes(fill = fin)) +
  theme_graph() +scale_fill_manual(name="Taxonomic Groups",values=c(taxCols))+scale_edge_color_manual(values=c("grey60","red"),breaks=c(1,-1))+geom_node_text(aes(label = namez),size=4,fontface="bold")
dev.off()

# Legend
ggplot(taxz,aes(x=1:11,y=2:12,fill=fin))+geom_point(size=6,shape=21)+scale_fill_manual(name="Taxonomic Group",values=c(taxCols))+theme_classic()+theme(legend.position="bottom")
ggsave("../Legend.pdf")
