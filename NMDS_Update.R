### SPOT Network Analysis ###
### Taxa barplots and NMDS ###
### By: Samantha Gleich ###
### Last Updated: 2/7/23 ###

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
# normDf <- clr(dfTax)
# normDf <- as.data.frame(normDf)
env <- read.csv(file.choose(),header=TRUE)
colz <- colsplit(rownames(normDf),"_",c("spot","Cruise","month","day","year","Depth"))
colz <- left_join(colz,env)
mytable <- colz[c(10:32)]
# mytable$CSDepth <- NULL
mytable <- missForest::missForest(mytable)
mytable <- mytable$ximp

decorana (normDf) # Axis length > 3 means CCA vs RDA

ccaOut <- cca(normDf~.,mytable)
ccaDf <- data.frame(ccaOut$CCA$u)
try <- ccaOut$CCA$eig
try

ccaDf$Month <- colz$Month
#colre <- randomcoloR::distinctColorPalette(12)
ccaDf[is.na(ccaDf)] <- 2

finalmodel<- ordistep(ccaOut, scope=formula(ccaOut))
vif.cca(finalmodel)
anovOut <- anova.cca(finalmodel, by="terms")
anovOutSub <- subset(anovOut,`Pr(>F)` < 0.01)
keep <- rownames(anovOutSub)

fitDf <- data.frame(ccaOut$CCA$biplot)
fitDf <- fitDf[c(1:2)]
fitDf$x <- mean(ccaDf$CCA1)
fitDf$y <- mean(ccaDf$CCA2)

t <- envfit(ccaOut,mytable)
t <- data.frame(t$vectors$arrows)
fitDf$CCA1 <- t$CCA1
fitDf$CCA2 <- t$CCA2

fitDf <- subset(fitDf,rownames(fitDf) %in% keep)

outP <- ggplot(ccaDf,aes(CCA1,CCA2,fill=factor(Month),shape=as.factor(Month)))+geom_point(size=2,color="black")+scale_fill_manual(name="Month",values=c("dodgerblue","dodgerblue","darkolivegreen4","darkolivegreen4","darkolivegreen4","goldenrod1","goldenrod1","goldenrod1","firebrick2","firebrick2","firebrick2","dodgerblue"),labels=c("January","February","March","April","May","June","July","August","September","October","November","December"))+scale_shape_manual(name="Month",values=c(21,22,21,22,24,21,22,24,21,22,24,24),labels=c("January","February","March","April","May","June","July","August","September","October","November","December"))+theme_classic()+xlab("CCA1 (41.5%)")+ylab("CCA2 (34.5%)")+geom_vline(xintercept = 0,linetype="dotted")+geom_hline(yintercept = 0,linetype="dotted")+ggtitle("DCM")+geom_segment(aes(x = fitDf[1,3], y = fitDf[1,4], xend = fitDf[1,1] , yend = fitDf[1,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[2,3], y = fitDf[2,4], xend = fitDf[2,1] , yend = fitDf[2,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[3,3], y = fitDf[3,4], xend = fitDf[3,1] , yend = fitDf[3,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[4,3], y = fitDf[4,4], xend = fitDf[4,1] , yend = fitDf[4,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[5,3], y = fitDf[5,4], xend = fitDf[5,1] , yend = fitDf[5,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[6,3], y = fitDf[6,4], xend = fitDf[6,1] , yend = fitDf[6,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)#+geom_segment(aes(x = fitDf[7,3], y = fitDf[7,4], xend = fitDf[7,1] , yend = fitDf[7,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)#+geom_segment(aes(x = fitDf[8,3], y = fitDf[8,4], xend = fitDf[8,1] , yend = fitDf[8,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)
outP
outP2+outP+plot_layout(ncol=2,guides = "collect")
ggsave("../../DCM.pdf",width=8,height=6)
