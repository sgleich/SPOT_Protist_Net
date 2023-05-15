### SPOT Network Analysis ###
### Figure 5: Bray-curtis similarity vs. temeprature distance ###
### By: Samantha Gleich ###
### Last Updated: 5/12/23 ###

# Libraries
library(lubridate)
library(reshape2)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ape)
library(ggpmisc)
library(patchwork)

# Set SPOT working directory
setwd("~/Desktop/SPOT/SPOT_2023")

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
env <- read.csv("SPOT_Env_NewJan11.csv",header=TRUE)
datez <- left_join(datez,env)
datez2 <- colsplit(datez$Date,"/",c("month","day","year"))
dfNew$month <- datez2$month
dfNew$day <- datez2$day
dfNew$year <- datez2$year
dfNew$date <- as.Date(with(dfNew, paste(year, month, day,sep="-")), "%y-%m-%d")

# Remove month, day, year columns since we have a date column now
dfNew$month <- NULL
dfNew$day <- NULL
dfNew$year <- NULL
dfNew <- as.data.frame(t(dfNew))

rownames(datez) <- paste(datez$spot,datez$Cruise,datez$m,datez$d,datez$y,datez$Depth,sep="_")
datez <- subset(datez,rownames(datez)!="SPOT_115_2_16_12_5m"& rownames(datez)!="SPOT_115_2_16_12_DCM")

keepEnvSurf <- c("DayOfYear","O2Wink","NH4","SiO3","PO4","NO2.NO3","CTDTMP","CTDBEAM","CTDFLUOR","CTDOXY","CTDSAL","DayLength","DayDiff","SLA","cyanos","sar11","cyanos_PA","sar11_PA","MEI","SST","Chla","PP")

keepEnvDCM <- c("DayOfYear","CSDepth","O2Wink","NH4","SiO3","PO4","NO2.NO3","CTDTMP","CTDBEAM","CTDFLUOR","CTDOXY","CTDSAL","DayLength","DayDiff","SLA","cyanos","sar11","cyanos_PA","sar11_PA","MEI","SST","Chla","PP")

envSurf <- subset(datez,select=c(keepEnvSurf))
envDCM <- subset(datez,select=c(keepEnvDCM))

envSurf <- subset(envSurf,grepl("5m",rownames(envSurf)))
envDCM <- subset(envDCM,grepl("DCM",rownames(envDCM)))

set.seed(100)
envTableSurf <- missForest::missForest(envSurf)
envTableSurf <- data.frame(envTableSurf$ximp)

set.seed(100)
envTableDCM <- missForest::missForest(envDCM)
envTableDCM <- data.frame(envTableDCM$ximp)

# Subset 5m and DCM samples
asvSurf <- subset(dfNew,select=c(grepl("5m",colnames(dfNew))))
asvDCM <- subset(dfNew,select=c(grepl("DCM",colnames(dfNew))))

# Make sequence of dates for day of time-series parameter
dayz <- seq(as.Date("2003/9/8"), as.Date("2019/1/1"), "days")
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

# Set up data for plot of bray-curtis similarity vs. lag
braySurfDf <- as.data.frame(braySurf)
meltBray <- melt(braySurf)
distMatDf <- as.data.frame(distMat)
meltDist <- melt(distMat)
colnames(meltDist)<- c("Var1","Var2","dist")

pickEnv <- function(df,envParam){
  total <- left_join(meltBray,meltDist)
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
  newEnv <- subset(df,select=c(paste(envParam)))
  newEnv <- vegdist(newEnv,method="euclidean")
  newEnv <- as.matrix(newEnv)
  newEnv <- as.data.frame(newEnv)
  newEnv$var <- rownames(newEnv)
  newEnv <- melt(newEnv,id.vars="var")
  colnames(newEnv)<- c("Var1","Var2","param")
  total <- left_join(total,newEnv)
  return(total)} 

ctdtmp <- pickEnv(envTableSurf,"SST")
ctdtmp2 <- ctdtmp %>% filter(year!=0) %>% group_by(year) %>% summarize(mBray=mean(value),sdBray=sd(value),mParam=mean(param),sdParam=sd(param)) %>% as.data.frame()
ctdtmpNum <- ctdtmp %>% filter(year!=0) %>% group_by(year) %>% tally()
ctdtmp2 <- left_join(ctdtmp2,ctdtmpNum)
ctdtmp2$sdBray <- ctdtmp2$sdBray/ctdtmp2$n
ctdtmp2$sdParam <- ctdtmp2$sdParam/ctdtmp2$n

ggplot(ctdtmp2,aes(x=mParam,y=mBray,label=as.character(year)))+geom_point()+theme_classic()+geom_errorbar(aes(ymin=mBray-sdBray, ymax=mBray+sdBray), width=.01)+geom_errorbar(aes(xmin=mParam-sdParam, xmax=mParam+sdParam), width=.001)+geom_text(hjust=-0.95, vjust=1.2)+ylim(0.25,0.35)+xlab("Euclidean Temperature Distance (°C)")+ylab("Bray-curtis Similarity")+geom_smooth(method="gam",linewidth=0.5,se=FALSE,color="grey40")
ggsave("Figure5_May2023.pdf")

gamTry <- mgcv::gam(ctdtmp2$mBray~te(ctdtmp2$mParam,k=3),data=ctdtmp2,method="REML")
summary(gamTry)
plot(gamTry,all.terms = TRUE)
