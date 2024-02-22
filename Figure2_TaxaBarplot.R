### SPOT Network Analysis ###
### Figure 2: Taxa barplots ###
### By: Samantha Gleich ###
### Last Updated: 2/2/24 ###

# Load libraries
library(tidyverse)
library(reshape2)
library(ggplot2)
library(randomcoloR)
library(vegan)

# Set working directory
setwd("~/Desktop/export_dir_feb2024")

# Load in counts data and manifest file
counts <- read.delim("feature-table.tsv",header=TRUE,row.names=1,skip=1)
manifest <- read.delim("manifest.txt",header=FALSE,sep=",")

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
counts$sum <- rowSums(counts)
counts <- counts/counts$sum
counts$sum <- NULL
rowSums(counts)
counts$namez <- rownames(counts)
dfNew <- left_join(counts,fin)
rownames(dfNew) <- dfNew$new
dfNew$namez <- NULL
dfNew$new <- NULL

# Add taxonomy information
tax <- read.delim("taxonomy.tsv",header=TRUE)
tax <- tax[c(1:2)]
dfNew <- as.data.frame(t(dfNew))
dfNew$Feature.ID <- rownames(dfNew)

# Let's add the taxonomy information associated with each ASV to our dataframe
# tax <- read.delim(file.choose(),header=TRUE)
# tax <- tax[c(1:2)]
# dfNew$Feature.ID <- rownames(dfNew)
dfTax <- left_join(dfNew,tax)

# Now let's choose specific taxonomic groups to visualize in our plots
dfTax <- subset(dfTax,grepl("Eukaryot",dfTax$Taxon))
taxz <- colsplit(dfTax$Taxon,";",c("Supergroup","Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxz$fin <- ifelse(taxz$Class=="Dinophyceae","Dinoflagellates",NA)
taxz$fin <- ifelse(taxz$Order=="Dino-Group-II","Group II Syndiniales",taxz$fin)
taxz$fin <- ifelse(taxz$Order=="Dino-Group-I","Group I Syndiniales",taxz$fin)
taxz$fin <- ifelse(taxz$Class=="Bacillariophyceae","Diatoms",taxz$fin)
taxz$fin <- ifelse(grepl("MAST",taxz$Family),"MAST groups",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Rhizaria","Rhizaria",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Chlorophyta","Chlorophytes",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Cryptophyta","Cryptophytes",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Haptophyta","Haptophytes",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Ciliophora","Ciliates",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Metazoa","Metazoans",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Stramenopiles" & is.na(taxz$fin),"Other Stramenopiles",taxz$fin)
taxz$fin <- ifelse(taxz$Domain=="Archaeplastida" & is.na(taxz$fin),"Other Archaeplastids",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Alveolata" & is.na(taxz$fin),"Other Alveolates",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin),"Other Eukaryotes",taxz$fin)
unique(taxz$fin)

# Add new taxonomy groups to main dataframe
dfTax$fin <- taxz$fin
dfTax$Taxon <- NULL

# Melt dataframe
dfMelt <- melt(dfTax,id.vars=c("fin","Feature.ID"))

# Split variable column up so that we can get month and year info
colz <- colsplit(dfMelt$variable,"_",c("SPOT","Num","Month","Day","Year","Depth"))
dfMelt$Cruise <- colz$Num
dfMelt$Depth <- colz$Depth

env <- read.csv("../SPOT/SPOT_2023/SPOT_Env_NewJan11.csv",header=TRUE)
dfMelt <- left_join(dfMelt,env)
dfMelt <- dfMelt[c(1,4,6,8:9)]
dfMelt <- subset(dfMelt,!is.na(Month))

# Summarize dataframe prior to plotting
dfSum <- dfMelt %>% group_by(fin,Month,Date,Depth) %>% summarize(s=sum(value)) %>% as.data.frame()

surfNum <- dfSum %>% filter(Depth=="5m" & fin=="Chlorophytes") %>% group_by(Month) %>% tally() %>% mutate(Depth="5m") %>% as.data.frame()
dcmNum <- dfSum %>% filter(Depth=="DCM" & fin=="Chlorophytes") %>% group_by(Month) %>% tally() %>% mutate(Depth="DCM") %>% as.data.frame()
numDf <- rbind(surfNum,dcmNum)

dfAvg <- dfSum %>% group_by(fin,Month,Depth) %>% summarize(m=mean(s))%>%as.data.frame()
dfAvg <- left_join(dfAvg,numDf)

# Let's label our months as opposed to having #s
dfAvg$Month <- factor(dfAvg$Month,levels=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c("January","February","March","April","May","June","July","August","September","October","November","December"))

dfAvg$Depth <- ifelse(dfAvg$Depth=="5m","Surface",dfAvg$Depth)
dfAvg$Depth <- factor(dfAvg$Depth,levels=c("Surface","DCM"))
dfAvg$Month2 <- paste(dfAvg$Month," (",dfAvg$n,")",sep="")

# colrs <- randomcoloR::distinctColorPalette(length(unique(dfAvg$fin)))
taxCols <- c("#E1746D","#76C3D7","#DE5AB1","#D5E0AF","#DED3DC","#87EB58","#D4DC60","#88E5D3","#88AAE1","#DBA85C","#8B7DDA","#9A8D87","#D99CD1","#B649E3","#7EDD90","#4FC4D0")
names(taxCols) <- c("Chlorophytes","Ciliates","Cryptophytes","Diatoms","Haptophytes","Dinoflagellates","MAST groups","Other Alveolates","Other Archaeplastids","Other Eukaryotes","Other Stramenopiles","Rhizaria","Group I Syndiniales","Group II Syndiniales","Unknown Eukaryotes","Metazoans")


dfAvg$Month2 <- factor(dfAvg$Month2,levels=c("January (11)","February (10)", "March (10)","April (10)","May (11)","June (10)","July (11)","July (9)","August (9)","September (10)","October (10)","November (10)","December (10)" ))

d <- dfAvg %>% filter(Depth=="DCM") %>% ggplot(aes(x=Month2,y=m))+geom_area(aes(fill = fin, group = fin),position="fill")+theme_classic()+xlab("Month")+ylab("Relative Abundance")+facet_wrap(~Depth)+theme(axis.text.x = element_text(angle=45, hjust=1))+scale_fill_manual(name="Taxonomic Group",values=c(taxCols))

s/d+plot_annotation(tag_levels="a")+plot_layout(guides="collect")
ggsave("../Figure2_Feb2024.pdf",width=6,height=8)
