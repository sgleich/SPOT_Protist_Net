### SPOT Network Analysis ###
### Figure 3: Bray-curtis similarity as a function of months between samples ###
### By: Samantha Gleich ###
### Last Updated: 2/23/24 ###

# Libraries
library(lubridate)
library(reshape2)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ape)
library(patchwork)

# Set SPOT working directory
setwd("~/Desktop/export_dir_2024")

# Read in ASVs and set up date info
counts <- read.delim("feature-table.tsv",header=TRUE,row.names=1,skip=1)
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

# Let's add the new (meaningful) sample names to our counts dataframe counts <- as.data.frame(t(counts))
counts$namez <- rownames(counts)
dfNew <- left_join(counts,fin)
rownames(dfNew) <- dfNew$new
dfNew$namez <- NULL
dfNew$new <- NULL

dfNew <- subset(dfNew,rownames(dfNew)!="SPOT_115_2_16_12_5m"& rownames(dfNew)!="SPOT_115_2_16_12_DCM")

# Make a "date" column in ASV counts dataframe
datez <- colsplit(rownames(dfNew),"_",c("spot","Cruise","m","d","y","Depth"))
env <- read.csv("SPOT_Env_Feb2024_Final.csv",header=TRUE)
datez <- left_join(datez,env)
datez2 <- colsplit(datez$Date,"/",c("month","day","year"))
dfNew$month <- datez2$month
dfNew$day <- datez2$day
dfNew$year <- datez2$year
dfNew$date <- as.Date(with(dfNew, paste(year, month, day,sep="-")), "%y-%m-%d")
dfNew$row <- paste("SPOT",datez$Cruise,datez$Depth,as.character(dfNew$date),sep="_")

# Remove month, day, year columns since we have a date column now
dfNew$month <- NULL
dfNew$day <- NULL
dfNew$year <- NULL
rownames(dfNew) <- dfNew$row
dfNew$row <- NULL
dfNew <- as.data.frame(t(dfNew))

tax <- read.delim("taxonomy.tsv",header=TRUE)
tax <- data.frame(name=c(tax$Feature.ID),tax=c(tax$Taxon))
dfNew$name <- rownames(dfNew)
v <- c("date","date")
tax <- rbind(tax,v)
dfNew <- left_join(dfNew,tax)
dfNew <- subset(dfNew,grepl("Eukaryota",dfNew$tax) |grepl("date",dfNew$tax) )
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

# Surface temporal distance matrix
surfDate <- as.data.frame(surfDate)
colnames(surfDate) <- "dayz"
surfDate$dayz <- as.Date(surfDate$dayz)
surfDate <- left_join(surfDate,dayz)
surfDate$names <- rownames(asvSurf)
surfDate$dayz <- NULL
rownames(surfDate) <- surfDate$names
surfDate$names <- NULL
distMat <- vegdist(surfDate,method="euclidean")
distMat <- as.matrix(distMat)

# Surface mantel test
braySurfDis <- vegdist(asvSurf,method="bray")
braySurfDis <- as.matrix(braySurfDis)
surfM <- mantel.test(distMat, braySurfDis, method = "spearman", permutations = 1000, na.rm = TRUE)

# Set up data for plot of bray-curtis similarity vs. lag
braySurfDf <- as.data.frame(braySurf)
meltBray <- melt(braySurf)
distMatDf <- as.data.frame(distMat)
meltDist <- melt(distMat)
colnames(meltDist)<- c("Var1","Var2","dist")

# Month lag
total <- left_join(meltBray,meltDist)

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

# DCM temporal distance matrix
dcmDate <- as.data.frame(dcmDate)
colnames(dcmDate) <- "dayz"
dcmDate$dayz <- as.Date(dcmDate$dayz)
dcmDate <- left_join(dcmDate,dayz)
dcmDate$names <- rownames(asvDCM)
dcmDate$dayz <- NULL
rownames(dcmDate) <- dcmDate$names
dcmDate$names <- NULL
distMat2 <- vegdist(dcmDate,method="euclidean")
distMat2 <- as.matrix(distMat2)

# DCM mantel test
brayDCMDis <- vegdist(asvDCM,method="bray")
brayDCMDis <- as.matrix(brayDCMDis)
dcmM <- mantel.test(distMat2, brayDCMDis, method = "spearman", permutations = 1000, na.rm = TRUE)

# Set up data for plot of bray-curtis similarity vs. lag
brayDCMDf <- as.data.frame(brayDCM)
meltBray2 <- melt(brayDCM)
distMatDf2 <- as.data.frame(distMat2)
meltDist2 <- melt(distMat2)
colnames(meltDist2)<- c("Var1","Var2","dist")

# Month Lag
total2 <- left_join(meltBray2,meltDist2)

# Combine surface and DCM
total$depth <- "Surface"
total2$depth <- "DCM"
totalAll <- rbind(total,total2)

totalAll$month <- totalAll$dist/30.4
totalAll$month <- round(totalAll$month)

totalAll2 <- totalAll[c(1:2)]
newDf <- data.frame(t(apply(totalAll2,1,sort)))
newDf <- newDf[!duplicated(newDf),]
colnames(newDf) <- c("Var1","Var2")
totalAll <- left_join(newDf,totalAll)

# Summarize data
totalNum <- totalAll%>%filter(month!=0)%>%group_by(month,depth)%>%tally()
total.month <- totalAll%>%filter(month!=0)%>%group_by(month,depth)%>%summarize(m=mean(value),s=sd(value))
total.month <- left_join(total.month,totalNum)
total.month$se <- total.month$s/sqrt(total.month$n)

s <- seq(from=0,to=177,by=6)

# Plot Bray-curtis similarity vs month lag
dcmP <- total.month%>%filter(depth=="DCM") %>% ggplot(aes(x=month,y=m))+geom_point(shape=21,fill="grey",size=2)+geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2,position=position_dodge(.1))+theme_classic()+scale_x_continuous(breaks=c(s))+xlab("Number of Months Between Samples")+ylab("Mean Bray-Curis Similarity (+/- SE)")+ggtitle("DCM")+theme(axis.text.x = element_text(angle = 45,hjust=1))+ylim(0.08,0.4)


surfP+dcmP+plot_layout(ncol=1)+plot_layout(guides = "collect")+plot_annotation(tag_levels="a")

#ggsave("Figure3.pdf",width=6,height=8)
