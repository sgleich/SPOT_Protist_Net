### SPOT Network Analysis ###
### Add SCC to Glasso-predicted edges ###
### By: Samantha Gleich ###
### Last Updated: 4/3/23 ###

# Load libraries
library(reshape2)
library(psych)
library(tidyverse)

# Add SCC as weights to glasso networks
s5 <- read.csv("s5_SPOT_APRIL2023_80p.csv",header=TRUE,row.names=1)
sDCM <- read.csv("sDCM_SPOT_APRIL2023_80p.csv",header=TRUE,row.names=1)

# Counts
counts <- read.delim("feature-table.tsv",header=TRUE,row.names=1)
counts <- as.data.frame(t(counts))

# Manifest
man <- read.delim("manifest.txt",header=FALSE,sep=",")
namez <- colsplit(man$V2,"-",c("a","b","c","d","e","f","g"))
namez$depth <- ifelse(grepl("5m",namez$g),"5m","DCM")
namez$fin <- paste(namez$c,namez$d,namez$e,namez$f,namez$depth,sep="_")
fin <- data.frame(namez=c(as.character(man$V1)),new=c(namez$fin))

fin$new <- paste("SPOT",fin$new,sep="_")
fin <- fin %>% distinct(new,.keep_all = TRUE)

# Join
counts$namez <- rownames(counts)
dfNew <- left_join(counts,fin)
rownames(dfNew) <- dfNew$new
dfNew$namez <- NULL
dfNew$new <- NULL

colnames(s5) <- str_remove_all(colnames(s5),"S_")
rownames(s5) <- str_remove_all(rownames(s5),"S_")
colnames(sDCM) <- str_remove_all(colnames(sDCM),"S_")
rownames(sDCM) <- str_remove_all(rownames(sDCM),"S_")

dfNew <- subset(dfNew,rownames(dfNew)!="SPOT_115_2_16_12_5m"& rownames(dfNew)!="SPOT_115_2_16_12_DCM")

dfCLR <- as.data.frame(clr(dfNew))
dfCLR <- as.data.frame(t(dfCLR))

# Change to 5 or DCM
dfCLR <- subset(dfCLR,rownames(dfCLR) %in% rownames(sDCM))
dfCLR <- as.data.frame(t(dfCLR))

dfCLR <- dfCLR[colnames(sDCM)]

dfCLR5 <- subset(dfCLR,grepl("5m",rownames(dfCLR)))
dfCLRDCM <- subset(dfCLR,grepl("DCM",rownames(dfCLR)))

dfCor5m <- corr.test(dfCLR5, method="spearman",adjust="bonferroni")
corOut5m <- as.data.frame(dfCor5m$r)

dfCorDcm <- corr.test(dfCLRDCM, method="spearman",adjust="bonferroni")
corOutDcm <- as.data.frame(dfCorDcm$r)

s5New <- s5*corOut5m
sDCMNew <- sDCM*corOutDcm

write.csv(s5New,"Glasso_5m_SCC_APRIL2023.csv")
write.csv(sDCMNew,"Glasso_DCM_SCC_APRIL2023.csv")

s5New <- as.matrix(s5New)
g5 <- graph_from_adjacency_matrix(s5New,mode="undirected",weighted=TRUE)

sDCMNew <- as.matrix(sDCMNew)
gd <- graph_from_adjacency_matrix(sDCMNew,mode="undirected",weighted=TRUE)

s5Edge <- data.frame(get.edgelist(g5))
s5Edge$weight <- E(g5)$weight
s5Edge$weight <-ifelse(s5Edge$weight>0,1,-1)

sDCMEdge <- data.frame(get.edgelist(gd))
sDCMEdge$weight <- E(gd)$weight
sDCMEdge$weight <-ifelse(sDCMEdge$weight>0,1,-1)

s5Edge$fin <- paste(s5Edge$X1,s5Edge$X2,s5Edge$weight,sep="_")
sDCMEdge$fin1 <- paste(sDCMEdge$X1,sDCMEdge$X2,sDCMEdge$weight,sep="_")
sDCMEdge$fin2 <- paste(sDCMEdge$X2,sDCMEdge$X1,sDCMEdge$weight,sep="_")

sInt <- subset(s5Edge,fin %in% sDCMEdge$fin1)
sInt2 <- subset(s5Edge,fin %in% sDCMEdge$fin2)
sInt <- rbind(sInt,sInt2)
sInt$fin <- NULL
write.csv(sInt,"InterEdge_APRIL2023.csv")
