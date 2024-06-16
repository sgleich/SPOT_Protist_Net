# Libraries
library(reshape2)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpmisc)
library(patchwork)
library(factoextra)

setwd("~/Desktop/export_dir_feb2024")

counts <- read.delim("feature-table.tsv",header=TRUE,skip=1)
manifest <- read.delim("manifest.txt",header=FALSE)
manifest_cols <- colsplit(manifest$V1,",",c("name","sample","direction"))
manifest_cols <- manifest_cols %>% distinct(name,.keep_all = TRUE) %>% as.data.frame()
rownames(counts) <- counts$X.OTU.ID
counts$X.OTU.ID <- NULL
d <- data.frame(name=colnames(counts))
d <- left_join(d,manifest_cols)
sample_cols <- colsplit(d$sample,"_",c("name","other"))
colnames(counts) <- sample_cols$name

env <- read.csv("SPOT_Env_Feb2024_Final.csv",header=TRUE)
env_keep <- c("Cruise","Depth","SLA","MEI","SST","Beam","Oxy","Sal","Temp","O2Wink","NH4","SiO3","PO4","NO2.NO3","DayOfYear","DayLength","DayDiff","CSDepth")

env_5 <- subset(env,select=c(env_keep))
env_5 <- subset(env_5,Depth=="5m")
env_5$CSDepth <- NULL
env_5$Depth <- NULL
rownames(env_5) <- paste("SPOT",env_5$Cruise,sep="_")
env_5$Cruise <- NULL

env_dcm <- subset(env,select=c(env_keep))
env_dcm <- subset(env_dcm,Depth=="DCM")
env_dcm$SiO3 <- NULL
env_dcm$Depth <- NULL
rownames(env_dcm) <- paste("SPOT",env_dcm$Cruise,sep="_")
env_dcm$Cruise <- NULL


set.seed(100)
envTable5 <- missForest::missForest(env_5)
envTable5 <- data.frame(envTable5$ximp)
envTable5 <- decostand(envTable5, method = "standardize")


set.seed(100)
envTabledcm <- missForest::missForest(env_dcm)
envTabledcm <- data.frame(envTabledcm$ximp)
envTabledcm <- decostand(envTabledcm, method = "standardize")


pca5 <- prcomp(envTable5)
pcaDCM <- prcomp(envTabledcm)
pca5Df <- data.frame(pca5$x)
rownames(pca5Df) <- rownames(envTable5)
pcaDCMDf <- data.frame(pcaDCM$x)
rownames(pcaDCMDf) <- rownames(envTabledcm)
pca5Dist <- vegdist(pca5Df$PC1,method="euclidean")
pca5Dist <- as.matrix(pca5Dist)
rownames(pca5Dist) <- rownames(pca5Df)
colnames(pca5Dist) <- rownames(pca5Df)
pcaDCMDist <- vegdist(pcaDCMDf$PC1,method="euclidean")
pcaDCMDist <- as.matrix(pcaDCMDist)
rownames(pcaDCMDist) <- rownames(pcaDCMDf)
colnames(pcaDCMDist) <- rownames(pcaDCMDf)

counts <- as.data.frame(t(counts))
dfNew <- subset(counts,rownames(counts) !="SPOTRe-DNA-115-02-16-12-5m" & rownames(counts)!="SPOTRe-DNA-115-02-16-12-DCM")

tax <- read.delim("taxonomy.tsv",header=TRUE)
tax <- data.frame(name=c(tax$Feature.ID),tax=c(tax$Taxon))
dfNew <- as.data.frame(t(dfNew))
dfNew$name <- rownames(dfNew)
dfNew <- left_join(dfNew,tax)
dfNew <- subset(dfNew,grepl("Eukaryota",dfNew$tax))
rownames(dfNew) <- dfNew$name
dfNew$name <- NULL
dfNew$tax <- NULL

# Subset 5m and DCM samples
asvSurf <- subset(dfNew,select=c(grepl("5m",colnames(dfNew))))
asvDCM <- subset(dfNew,select=c(grepl("DCM",colnames(dfNew))))

# Make sequence of dates for day of time-series parameter
dayz <- seq(as.Date("2003/9/18"), as.Date("2019/1/1"), "days")
dayz<- as.data.frame(dayz)
dayz$num <- 1:nrow(dayz)

# Surface bray-curtis dissimilarity matrix
asvSurf <- as.data.frame(t(asvSurf))
surfDate <- asvSurf$date
asvSurf$date <- NULL
asvSurf <- mutate_all(asvSurf, function(x) as.numeric(as.character(x)))
asvSurf <- decostand(asvSurf,method="total")
braySurf <- vegdist(asvSurf,method="bray")
braySurf <- 1-braySurf
range(braySurf)
braySurf <- as.matrix(braySurf)

# Set up data for plot of bray-curtis similarity vs. lag
braySurfDf <- as.data.frame(braySurf)
meltBray <- melt(braySurf)


# Bray-cutris dissimilarity matrix 
asvDCM <- as.data.frame(t(asvDCM))
dcmDate <- asvDCM$date
asvDCM$date <- NULL
asvDCM <- mutate_all(asvDCM, function(x) as.numeric(as.character(x)))
asvDCM <- decostand(asvDCM,method="total")
brayDCM <- vegdist(asvDCM,method="bray")
brayDCM <- 1-brayDCM
range(brayDCM)
brayDCM <- as.matrix(brayDCM)

# Set up data for plot of bray-curtis similarity vs. lag
brayDCMDf <- as.data.frame(brayDCM)
meltBrayDCM <- melt(brayDCM)

meltpca5Dist <- melt(pca5Dist)
meltpcadcmDist <- melt(pcaDCMDist)

cols <- colsplit(meltBray$Var1,"-",c("SPOT","DCM","Cruise","Other"))
cols2 <- colsplit(meltBray$Var2,"-",c("SPOT","DCM","Cruise","Other"))
meltBray$Var1 <- paste("SPOT",cols$Cruise,sep="_")
meltBray$Var2 <- paste("SPOT",cols2$Cruise,sep="_")
meltBray <- subset(meltBray,value >0)
meltBrayB <- meltBray[c(1:2)]
meltBrayB <- data.frame(t(apply(meltBrayB,1,sort)))
meltBrayB <- meltBrayB[!duplicated(meltBrayB),]
colnames(meltBrayB) <- c("Var1","Var2")
meltBray <- left_join(meltBrayB,meltBray)

meltpca5Dist <- subset(meltpca5Dist,value >0)
meltpca5DistB <- meltpca5Dist[c(1:2)]
meltpca5DistB <- data.frame(t(apply(meltpca5DistB,1,sort)))
meltpca5DistB <- meltpca5DistB[!duplicated(meltpca5DistB),]
colnames(meltpca5DistB) <- c("Var1","Var2")
meltpca5Dist <- left_join(meltpca5DistB,meltpca5Dist)
colnames(meltpca5Dist)[3] <- "PCA"

meltBray <- left_join(meltBray,meltpca5Dist)


p5 <- ggplot(meltBray,aes(x=PCA,y=value))+geom_point(size=2,shape=21,fill="grey90")+theme_classic()+xlab("Euclidean PC1 Distance")+ylab("Bray-Curtis Similarity")+stat_poly_line(color="red")+stat_poly_eq(use_label(c("eq", "R2")),hjust = -0.5)+ggtitle("Surface")

var5 <- get_pca_var(pca5)


cols <- colsplit(meltBrayDCM$Var1,"-",c("SPOT","DCM","Cruise","Other"))
cols2 <- colsplit(meltBrayDCM$Var2,"-",c("SPOT","DCM","Cruise","Other"))
meltBrayDCM$Var1 <- paste("SPOT",cols$Cruise,sep="_")
meltBrayDCM$Var2 <- paste("SPOT",cols2$Cruise,sep="_")
meltBrayDCM <- subset(meltBrayDCM,value >0)
meltBrayDCMB <- meltBrayDCM[c(1:2)]
meltBrayDCMB <- data.frame(t(apply(meltBrayDCMB,1,sort)))
meltBrayDCMB <- meltBrayDCMB[!duplicated(meltBrayDCMB),]
colnames(meltBrayDCMB) <- c("Var1","Var2")
meltBrayDCM <- left_join(meltBrayDCMB,meltBrayDCM)

meltpcaDCMDist <- subset(meltpcaDCMDist,value >0)
meltpcaDCMDistB <- meltpcaDCMDist[c(1:2)]
meltpcaDCMDistB <- data.frame(t(apply(meltpcaDCMDistB,1,sort)))
meltpcaDCMDistB <- meltpcaDCMDistB[!duplicated(meltpcaDCMDistB),]
colnames(meltpcaDCMDistB) <- c("Var1","Var2")
meltpcaDCMDist <- left_join(meltpcaDCMDistB,meltpcaDCMDist)
colnames(meltpcaDCMDist)[3] <- "PCA"

meltBrayDCM <- left_join(meltBrayDCM,meltpcaDCMDist)


pdcm <- ggplot(meltBrayDCM,aes(x=PCA,y=value))+geom_point(size=2,shape=21,fill="grey90")+theme_classic()+xlab("Euclidean PC1 Distance")+ylab("Bray-Curtis Similarity")+stat_poly_line(color="red")+stat_poly_eq(use_label(c("eq", "R2")),hjust = -0.5)+ggtitle("DCM")

vardcm <- get_pca_var(pcaDCM)

p5+pdcm
ggsave("../NewSPOT.pdf",width=12,height=5)
