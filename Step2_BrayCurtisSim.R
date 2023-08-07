### SPOT Network Analysis ###
### Figure 1: Bray-curtis similarity as a function of month and year lags ###
### By: Samantha Gleich ###
### Last Updated: 8/7/23 ###

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
setwd("~/Desktop/SPOT/export_dir")

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

# Let's add the new (meaningful) sample names to our counts dataframe counts <- as.data.frame(t(counts))
counts$namez <- rownames(counts)
dfNew <- left_join(counts,fin)
rownames(dfNew) <- dfNew$new
dfNew$namez <- NULL
dfNew$new <- NULL

dfNew <- subset(dfNew,rownames(dfNew)!="SPOT_115_2_16_12_5m"& rownames(dfNew)!="SPOT_115_2_16_12_DCM")

# Make a "date" column in ASV counts dataframe
datez <- colsplit(rownames(dfNew),"_",c("spot","Cruise","m","d","y","Depth"))
env <- read.csv("../SPOT_Env_NewJan11.csv",header=TRUE)
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
total$month <- ifelse(total$dist>=15 & total$dist<=45,1,0)
total$month <- ifelse(total$dist>=46 & total$dist<=76,2,total$month)
total$month <- ifelse(total$dist>=77 & total$dist<=107,3,total$month)
total$month <- ifelse(total$dist>=108 & total$dist<=138,4,total$month)
total$month <- ifelse(total$dist>=139 & total$dist<=169,5,total$month)
total$month <- ifelse(total$dist>=170 & total$dist<=200,6,total$month)
total$month <- ifelse(total$dist>=201 & total$dist<=231,7,total$month)
total$month <- ifelse(total$dist>=232 & total$dist<=262,8,total$month)
total$month <- ifelse(total$dist>=263 & total$dist<=293,9,total$month)
total$month <- ifelse(total$dist>=294 & total$dist<=324,10,total$month)
total$month <- ifelse(total$dist>=325 & total$dist<=355,11,total$month)

# Year Lag
total$year <- ifelse(total$dist>=350 & total$dist<=380,1,0)
total$year <- ifelse(total$dist>=715 & total$dist<=745,2,total$year)
total$year <- ifelse(total$dist>=1080 & total$dist<=1110,3,total$year)
total$year <- ifelse(total$dist>=1445 & total$dist<=1475,4,total$year)
total$year <- ifelse(total$dist>=1810 & total$dist<=1840,5,total$year)
total$year <- ifelse(total$dist>=2175 & total$dist<=2205,6,total$year)
total$year <- ifelse(total$dist>=2540 & total$dist<=2570,7,total$year)
total$year <- ifelse(total$dist>=2905 & total$dist<=2935,8,total$year)
total$year <- ifelse(total$dist>=3270 & total$dist<=3300,9,total$year)
total$year <- ifelse(total$dist>=3635 & total$dist<=3665,10,total$year)
total$year <- ifelse(total$dist>=4000 & total$dist<=4030,11,total$year)
total$year <- ifelse(total$dist>=4365 & total$dist<=4395,12,total$year)
total$year <- ifelse(total$dist>=4730 & total$dist<=4760,13,total$year)
total$year <- ifelse(total$dist>=5095 & total$dist<=5125,14,total$year)
total$year <- ifelse(total$dist>=5460 & total$dist<=5490,15,total$year)

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
total2$month <- ifelse(total2$dist>=15 & total2$dist<=45,1,0)
total2$month <- ifelse(total2$dist>=46 & total2$dist<=76,2,total2$month)
total2$month <- ifelse(total2$dist>=77 & total2$dist<=107,3,total2$month)
total2$month <- ifelse(total2$dist>=108 & total2$dist<=138,4,total2$month)
total2$month <- ifelse(total2$dist>=139 & total2$dist<=169,5,total2$month)
total2$month <- ifelse(total2$dist>=170 & total2$dist<=200,6,total2$month)
total2$month <- ifelse(total2$dist>=201 & total2$dist<=231,7,total2$month)
total2$month <- ifelse(total2$dist>=232 & total2$dist<=262,8,total2$month)
total2$month <- ifelse(total2$dist>=263 & total2$dist<=293,9,total2$month)
total2$month <- ifelse(total2$dist>=294 & total2$dist<=324,10,total2$month)
total2$month <- ifelse(total2$dist>=325 & total2$dist<=355,11,total2$month)

# Year Lag
total2$year <- ifelse(total2$dist>=350 & total2$dist<=380,1,0)
total2$year <- ifelse(total2$dist>=715 & total2$dist<=745,2,total2$year)
total2$year <- ifelse(total2$dist>=1080 & total2$dist<=1110,3,total2$year)
total2$year <- ifelse(total2$dist>=1445 & total2$dist<=1475,4,total2$year)
total2$year <- ifelse(total2$dist>=1810 & total2$dist<=1840,5,total2$year)
total2$year <- ifelse(total2$dist>=2175 & total2$dist<=2205,6,total2$year)
total2$year <- ifelse(total2$dist>=2540 & total2$dist<=2570,7,total2$year)
total2$year <- ifelse(total2$dist>=2905 & total2$dist<=2935,8,total2$year)
total2$year <- ifelse(total2$dist>=3270 & total2$dist<=3300,9,total2$year)
total2$year <- ifelse(total2$dist>=3635 & total2$dist<=3665,10,total2$year)
total2$year <- ifelse(total2$dist>=4000 & total2$dist<=4030,11,total2$year)
total2$year <- ifelse(total2$dist>=4365 & total2$dist<=4395,12,total2$year)
total2$year <- ifelse(total2$dist>=4730 & total2$dist<=4760,13,total2$year)
total2$year <- ifelse(total2$dist>=5095 & total2$dist<=5125,14,total2$year)
total2$year <- ifelse(total2$dist>=5460 & total2$dist<=5490,15,total2$year)

# Combine surface and DCM
total$depth <- "Surface"
total2$depth <- "DCM"
totalAll <- rbind(total,total2)

# Summarize data
totalNum <- totalAll%>%filter(month!=0)%>%group_by(month,depth)%>%tally()
total.month <- totalAll%>%filter(month!=0)%>%group_by(month,depth)%>%summarize(m=mean(value),s=sd(value))
total.month <- left_join(total.month,totalNum)
total.month$se <- total.month$s/sqrt(total.month$n)

# Plot Bray-curtis similarity vs month lag
monthPlot <- total.month%>%ggplot(aes(x=month,y=m,color=depth,shape=depth,linetype=depth))+geom_point()+geom_line()+geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2,position=position_dodge(.3))+theme_classic()+scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12))+xlab("Month Lag")+ylab("Mean Bray-Curis Similarity (+/- SE)")+ggtitle("")+scale_y_continuous(limits=c(0.20,0.40),breaks=c(0.2,0.25,0.3,0.35,0.4))+scale_color_manual(name="Depth",breaks=c("Surface","DCM"),values=c("black","grey60"))+scale_shape_manual(name="Depth",breaks=c("Surface","DCM"),values=c(19,15))+scale_linetype_manual(name="Depth",breaks=c("Surface","DCM"),values=c("solid","dashed"))

# Plot Bray-curtis similarity vs. year lag
total.year <- totalAll%>%filter(year!=0 & year !=15)%>%group_by(year,depth)%>%summarize(m=mean(value),s=sd(value))
totalNum <- totalAll%>%filter(year!=0 & year != 15)%>%group_by(year,depth)%>%tally()
total.year <- left_join(total.year,totalNum)
total.year$se <- total.year$s/sqrt(total.year$n)

yearPlot<- total.year%>%ggplot(aes(x=year,y=m,color=depth,shape=depth,linetype=depth))+geom_point()+geom_line()+geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.2,position=position_dodge(.3))+theme_classic()+scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14))+xlab("Year Lag")+ylab("Mean Bray-Curis Similarity (+/- SE)")+ggtitle("")+scale_y_continuous(limits=c(0.2,0.4),breaks=c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5))+scale_color_manual(name="Depth",breaks=c("Surface","DCM"),values=c("black","grey60"))+scale_shape_manual(name="Depth",breaks=c("Surface","DCM"),values=c(19,15))+scale_linetype_manual(name="Depth",breaks=c("Surface","DCM"),values=c("solid","dashed"))

monthPlot+yearPlot+plot_layout(ncol=1)+plot_layout(guides = "collect")+plot_annotation(tag_levels="a")

ggsave("Figure1_Aug2023.pdf",width=6,height=8)
