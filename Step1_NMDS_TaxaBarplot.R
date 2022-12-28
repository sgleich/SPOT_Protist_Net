### SPOT Network Analysis ###
### Taxa barplots and NMDS ###
### By: Samantha Gleich ###
### Last Updated: 11/11/22 ###

# Load libraries
library(tidyverse)
library(reshape2)
library(ggplot2)
library(randomcoloR)
library(vegan)

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

# Make NMDS plot with bray-curtis dissimilarity
normDf <- decostand(dfTax,method="total")
distDf <- vegdist(normDf, method = "bray")
distMat <- as.matrix(distDf, labels = T)
runNMDS <-metaMDS(distMat,distance = "bray",k = 2,maxit = 999, trymax = 500,wascores = TRUE)
runNMDS$stress # Record stress = 0.21
mdsDf <- data.frame(runNMDS$points)
colz <- colsplit(rownames(mdsDf),"_",c("spot","num","month","day","year","depth"))
mdsDf$Depth <- colz$depth

nmds <- ggplot(mdsDf,aes(MDS1,MDS2,fill=Depth))+geom_point(size=2,color="black",shape=21)+scale_fill_manual(values=c("grey70","grey20"),name="Depth",labels=c("5 m","DCM"))+guides(fill = guide_legend(override.aes=list(shape=21)))+theme_classic()+xlab("NMDS1")+ylab("NMDS2")+geom_vline(xintercept = 0,linetype="dotted")+geom_hline(yintercept = 0,linetype="dotted")+ggtitle("")+stat_ellipse(level=0.95,aes(color=Depth))+scale_color_manual(values=c("grey70","grey20"))
nmds
# ggsave("NMDS_Nov2022.pdf") # Save figure

# Let's add the taxonomy information associated with each ASV to our dataframe
# tax <- read.delim(file.choose(),header=TRUE)
# tax <- tax[c(1:2)]
# dfNew$Feature.ID <- rownames(dfNew)
dfTax <- left_join(dfNew,tax)

# Now let's choose specific taxonomic groups to visualize in our plots
dfTax <- subset(dfTax,grepl("Eukaryot",dfTax$Taxon))
taxz <- colsplit(dfTax$Taxon,";",c("Supergroup","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxz$fin <- ifelse(taxz$Class=="Syndiniales","Syndiniales",NA)
taxz$fin <- ifelse(taxz$Class=="Dinophyceae","Dinoflagellate",taxz$fin)
taxz$fin <- ifelse(taxz$Class=="Bacillariophyta","Diatom",taxz$fin)
taxz$fin <- ifelse(grepl("MAST",taxz$Class),"MAST",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Rhizaria","Rhizaria",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Chlorophyta","Chlorophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Haptophyta","Haptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Metazoa","Metazoa",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Stramenopiles" & is.na(taxz$fin),"Other Stramenopiles",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Archaeplastida" & is.na(taxz$fin),"Other Archaeplastids",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin),"Other Eukaryote",taxz$fin)
# unique(taxz$fin)

# Add new taxonomy groups to main dataframe
dfTax$fin <- taxz$fin
dfTax$Taxon <- NULL

# Melt dataframe
dfMelt <- melt(dfTax,id.vars=c("fin","Feature.ID"))

# Split variable column up so that we can get month and year info
colz <- colsplit(dfMelt$variable,"_",c("SPOT","Num","Month","Day","Year","Depth"))
dfMelt$Month <- colz$Month
dfMelt$Year <- colz$Year
dfMelt$Day <- colz$Day
dfMelt$Depth <- colz$Depth

# Summarize dataframe prior to plotting
dfSum <- dfMelt %>% group_by(fin,Month,Year,Day,Depth) %>% summarize(s=sum(value)) %>% as.data.frame()
dfAvg <- dfSum %>% group_by(fin,Month,Depth) %>% summarize(m=mean(s))%>%as.data.frame()

# Let's label our months as opposed to having #s
dfAvg$Month <- factor(dfAvg$Month,levels=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c("January","February","March","April","May","June","July","August","September","October","November","December"))

dfAvg$Depth <- ifelse(dfAvg$Depth=="5m","5 m",dfAvg$Depth)

colrs <- randomcoloR::distinctColorPalette(length(unique(dfAvg$fin)))
 
dfAvg %>% ggplot(aes(x=Month,y=m))+geom_area(aes(fill = fin, group = fin),position="fill")+theme_classic()+xlab("Month")+ylab("Relative Abundance")+facet_wrap(~Depth)+theme(axis.text.x = element_text(angle=45, hjust=1))+scale_fill_manual(name="Taxonomic Group",values=c(colrs))
# ggsave("Taxa_Barplot_Nov2022.pdf")
