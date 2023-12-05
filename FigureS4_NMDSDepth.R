### SPOT Network Analysis ###
### Figure S2: Surface and DCM NMDS plot ###
### By: Samantha Gleich ###
### Last Updated: 12/5/23 ###

# Libraries
library(reshape2)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpubr)
library(patchwork)

# Set SPOT working directory
setwd("/Users/samanthagleich/Desktop/SPOT/SPOT_2023")

# Read in ASVs and set up date info
counts <- read.delim("feature-table.tsv",header=TRUE,row.names=1)
manifest <- read.delim("manifest.txt",header=FALSE,sep=",")
counts <- as.data.frame(t(counts))

# We need to rename our samples in the counts dataframe so that they contain meaningful information
namez <- colsplit(manifest$V2,"-",c("a","b","c","d","e","f","g"))
namez$depth <- ifelse(grepl("5m",namez$g),"5m","DCM")
# Making a "fin" column that contains the date and depth of each sample
namez$fin <- paste(namez$c,namez$d,namez$e,namez$f,namez$depth,sep="_")
# Make a new dataframe with old sample name info and new sample name info
fin <- data.frame(namez=c(manifest$V1),new=c(namez$fin))
fin$new <- paste("SPOT",fin$new,sep="_")
fin <- fin %>% distinct(new,.keep_all = TRUE)

# Let's add the new (meaningful) sample names to our counts dataframe 
counts$namez <- rownames(counts)
dfNew <- left_join(counts,fin)
rownames(dfNew) <- dfNew$new
dfNew$namez <- NULL
dfNew$new <- NULL

dfNew <- subset(dfNew,rownames(dfNew)!="SPOT_115_2_16_12_5m"& rownames(dfNew)!="SPOT_115_2_16_12_DCM")

tax <- read.delim("taxonomy_90.tsv",header=TRUE)
tax <- subset(tax,grepl("Eukaryot",tax$Taxon))

dfNew <- as.data.frame(t(dfNew))
dfNew <- subset(dfNew,rownames(dfNew) %in% tax$Feature.ID)
dfNew <- as.data.frame(t(dfNew))

dfNorm <- decostand(dfNew,method="total")
dfBray <- vegdist(dfNorm,method="bray")

nmdsOut <- metaMDS(dfBray,distance="bray",k=2)
nmdsOut$stress #0.213 kind of high :/

nmdsDf <- data.frame(nmdsOut$points)
colz <- colsplit(rownames(nmdsDf),"_",c("spot","num","month","day","year","depth"))
nmdsDf$Depth <- colz$depth

ggplot(nmdsDf,aes(x=MDS1,y=MDS2,fill=Depth))+geom_point(shape=21,size=2)+geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+theme_classic()+scale_fill_manual(values=c("grey30","grey80"),labels=c("Surface (5 m)","DCM"))+xlab("NMDS1")+ylab("NMDS2")+stat_ellipse(aes(color=Depth),level=0.95)+scale_color_manual(values=c("grey30","grey80"),labels=c("Surface (5 m)","DCM"))
ggsave("../../FigureS2_NEW.pdf")

anosimOut <- vegan::anosim(dfBray,nmdsDf$Depth)
