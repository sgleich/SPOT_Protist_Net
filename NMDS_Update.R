### SPOT Network Analysis ###
### Taxa barplots and NMDS ###
### By: Samantha Gleich ###
### Last Updated: 1/3/23 ###

# Load libraries
library(tidyverse)
library(reshape2)
library(ggplot2)
library(randomcoloR)
library(vegan)

# Colors to use for all taxonomic analyses
taxCols <- c("#E1746D","#76C3D7","#DE5AB1","#D5E0AF","#DED3DC","#87EB58","#D4DC60","#88E5D3","#88AAE1","#DBA85C","#8B7DDA","#9A8D87","#D99CD1","#B649E3","#7EDD90")
names(taxCols) <- c("Chlorophyte","Ciliate","Cryptophyte","Diatom","Dinoflagellate","Fungi","Haptophyte","MAST","Other Alveolate","Other Archaeplastida","Other Eukaryote","Other Stramenopile","Rhizaria","Syndiniales","Unknown Eukaryote")

# Load in counts data and manifest file
counts <- read.delim(file.choose(),header=TRUE,row.names=1)
manifest <- read.delim(file.choose(),header=FALSE,sep=",")

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
counts <- as.data.frame(t(counts))
counts$namez <- rownames(counts)
dfNew <- left_join(counts,fin)
rownames(dfNew) <- dfNew$new
dfNew$namez <- NULL
dfNew$new <- NULL

# Add taxonomy information
tax <- read.delim(file.choose(),header=TRUE)
tax <- tax[c(1:2)]
dfNew <- as.data.frame(t(dfNew))
dfNew$Feature.ID <- rownames(dfNew)
dfTax <- left_join(dfNew,tax)
dfTax <- subset(dfTax,grepl("Eukaryot",dfTax$Taxon))
dfTax$Feature.ID <- NULL
dfTax$Taxon <- NULL
dfTax <- as.data.frame(t(dfTax))


dfTax <- subset(dfTax,grepl("DCM",rownames(dfTax)))
dfTax <- subset(dfTax,rownames(dfTax)!="SPOT_115_2_16_12_5m"& rownames(dfTax)!="SPOT_115_2_16_12_DCM")
# Make NMDS plot with bray-curtis dissimilarity
normDf <- decostand(dfTax,method="total")
distDf <- vegdist(normDf, method = "bray")
distMat <- as.matrix(distDf, labels = T)
runNMDS <-metaMDS(distMat,distance = "bray",k = 2,maxit = 999, trymax = 500,wascores = TRUE)
runNMDS$stress # Record stress = 0.21
mdsDf <- data.frame(runNMDS$points)
colz <- colsplit(rownames(mdsDf),"_",c("spot","Cruise","month","day","year","Depth"))
env <- read.csv(file.choose(),header=TRUE)
colz <- left_join(colz,env)


mdsDf$Month <- colz$Month
#colre <- randomcoloR::distinctColorPalette(12)
mdsDf[is.na(mdsDf)] <- 2

fun_color_range <- colorRampPalette(c("darkgreen","white","darkgoldenrod3")) 
my_colors <- fun_color_range(12)    

nmds1 <- ggplot(mdsDf,aes(MDS1,MDS2,fill=factor(Month)))+geom_point(size=2,color="black",shape=21)+scale_fill_manual(values=c(my_colors),labels=c("January","February","March","April","May","June","July","August","September","October","November","December"))+guides(fill = guide_legend(override.aes=list(shape=21)))+theme_classic()+xlab("NMDS1")+ylab("NMDS2")+geom_vline(xintercept = 0,linetype="dotted")+geom_hline(yintercept = 0,linetype="dotted")+ggtitle("")
nmds

mytable <- colz[c(11,13:19,21:27,28:34)]
mytable$DayLengthYest <- NULL
mytable <- missForest::missForest(mytable)
mytable <- mytable$ximp
fit <- envfit(runNMDS, mytable)

out <- p.adjust(fit$vectors$pvals,method="fdr")
out <- out[out < 0.01]
out <- data.frame(out)

fitDf <- data.frame(fit$vectors$arrows)
fitDf$r <- fit$vectors$r
fitDf <- subset(fitDf,rownames(fitDf) %in% rownames(out))
fitDf$originX <- mean(mdsDf$MDS1)
fitDf$originY <- mean(mdsDf$MDS2)
fitDf <- fitDf %>% arrange(desc(abs(r)))

outP <- ggplot(mdsDf,aes(MDS1,MDS2,fill=factor(Month)))+geom_point(size=2,color="black",shape=21)+scale_fill_manual(name="Month",values=c(my_colors),labels=c("January","February","March","April","May","June","July","August","September","October","November","December"))+guides(fill = guide_legend(override.aes=list(shape=21)))+theme_classic()+xlab("NMDS1")+ylab("NMDS2")+geom_vline(xintercept = 0,linetype="dotted")+geom_hline(yintercept = 0,linetype="dotted")+ggtitle("Surface (5 m)")+geom_segment(aes(x = fitDf[1,4], y = fitDf[1,5], xend = fitDf[1,1] , yend = fitDf[1,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[2,4], y = fitDf[2,5], xend = fitDf[2,1] , yend = fitDf[2,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[3,4], y = fitDf[3,5], xend = fitDf[3,1] , yend = fitDf[3,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[4,4], y = fitDf[4,5], xend = fitDf[4,1] , yend = fitDf[4,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[5,4], y = fitDf[5,5], xend = fitDf[5,1] , yend = fitDf[5,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)
outP
outP2

ggsave("../trySurf.pdf")

tableSurf <- fitDf
