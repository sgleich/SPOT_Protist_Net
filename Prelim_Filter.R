### SPOT Network Analysis ###
## Pre-processing (Filtering + GAM-transforming) ###
library(tidyverse)
library(compositions)
library(NetGAM)
library(mgcv)
library(reshape2)

df <- read.csv(file.choose(),header=TRUE,row.names=1)
df2 <- df[1:244,]
df2 <- mutate_all(df2, function(x) as.numeric(as.character(x)))

df2clr <- clr(df2)

df5 <- subset(df2,grepl("5m", rownames(df2)))
dfDCM<- subset(df2,grepl("DCM", rownames(df2)))

df5 <- as.data.frame(t(df5))
dfDCM <- as.data.frame(t(dfDCM))

df5b <- subset(df5,rowSums(df5==0) <= 25) # 25 is ~20% not zero
df5b$count <- rowSums(df5b==0)
df5b$count <- NULL

dfDCMb <- subset(dfDCM,rowSums(dfDCM==0) <= 25) # 25 is ~20% not zero
dfDCMb$count <- rowSums(dfDCMb==0)
df5b$count <- NULL

namez5 <- rownames(df5b)
namezDCM <- rownames(dfDCMb)
namez <- c(namez5,namezDCM)
namez <- unique(namez)

df2clr <- as.data.frame(df2clr)

dfclrFilt <- subset(df2clr,select=c(namez))
colz <- colsplit(rownames(dfclrFilt),"\\.",c("spot","dna","num","month","day","year","other"))

vec <- rep(1:12, length=192)
vec <- as.data.frame(vec)
vec$year <- rep(3:18, each=12)
vec$DayOTS <- 1:192
colnames(vec)<- c("month","year","day")

colz <- data.frame(colz$month,colz$year)
colnames(colz)<- c("month","year")

colz <- left_join(colz,vec)

dfclrFilt$month <- colz$month
dfclrFilt$day <- colz$day

df5clr <- subset(dfclrFilt,grepl("5m", rownames(dfclrFilt)))
dfDCMclr <- subset(dfclrFilt,grepl("DCM", rownames(dfclrFilt)))

# Remove dups 
df5clr <- subset(df5clr,rownames(df5clr)!="SPOTRe.DNA.115.02.16.12.5m_S31_L001_R1_trimmed.fastq")
dfDCMclr <- subset(dfDCMclr,rownames(dfDCMclr)!="SPOTRe.DNA.115.02.16.12.DCM_S32_L001_R1_trimmed.fastq")

df5Month <- df5clr$month
df5Day <- df5clr$day
df5clr$month <- NULL
df5clr$day <- NULL

netGAM5 <- netGAM.df(df5clr,MOY=df5Month,MCount=df5Day,clrt=FALSE)

dfDCMMonth <- dfDCMclr$month
dfDCMDay <- dfDCMclr$day
dfDCMclr$month <- NULL
dfDCMclr$day <- NULL

netGAMDCM <- netGAM.df(dfDCMclr,MOY=dfDCMMonth,MCount=dfDCMDay,clrt=FALSE)

netGAM5<- as.data.frame(t(netGAM5))
netGAMDCM <- as.data.frame(t(netGAMDCM))

write.csv(netGAM5,"SPOT_5m_Filtered_GAM_LM.csv")
write.csv(netGAMDCM,"SPOT_DCM_Filtered_GAM_LM.csv")

### Run eLSA ###
