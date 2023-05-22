### SPOT LASSO ANALYSIS ###
### Script by: Samantha Gleich ###
### Last modified: May 22, 2023 ###

# Load libraries
library(glmnet)
library(missForest)
library(tidyverse)
library(reshape2)
library(caret)
library(yardstick)
library(compositions)
library(Metrics)
library(vegan)

# Set SPOT working directory
setwd("~/Desktop/SPOT/SPOT_2023")

# Set seed for reproducibility
set.seed(100)

# Load in ASV table and manifest
counts <- read.delim("feature-table.tsv",header=TRUE,row.names=1)
manifest <- read.delim("manifest.txt",header=FALSE,sep=",")

# We need to rename our samples in the counts dataframe so that they contain meaningful information
namez <- colsplit(manifest$V2,"-",c("a","b","c","d","e","f","g"))
namez$depth <- ifelse(grepl("5m",namez$g),"5m","DCM")
# Making a "fin" column that contains the date and depth of each sample
namez$fin <- paste(namez$c,namez$d,namez$e,namez$f,namez$depth,sep="_")

# Make a new dataframe with old sample name info and new sample name info
fin <- data.frame(namez=c(as.character(manifest$V1)),new=c(namez$fin))
fin$new <- paste("SPOT",fin$new,sep="_")
fin <- fin %>% distinct(new,.keep_all = TRUE)

# Let's add the new (meaningful) sample names to our counts dataframe 
counts <- as.data.frame(t(counts))
counts$namez <- rownames(counts)
dfNew <- left_join(counts,fin)
rownames(dfNew) <- dfNew$new
dfNew$namez <- NULL
dfNew$new <- NULL
dfNew<- subset(dfNew,rownames(dfNew)!="SPOT_115_2_16_12_5m"& rownames(dfNew)!="SPOT_115_2_16_12_DCM")

# CLR transform data (because they're relative abundances)
dfCLR <- as.data.frame(clr(dfNew))

# Let's only keep those ASVs that we included in our networks
netOut5 <- read.csv("NetGAM_Out5_May16.csv",header=TRUE,row.names=1) 
netOutDCM <- read.csv("NetGAM_OutDCM_May16.csv",header=TRUE,row.names=1)
namez <- c(colnames(netOut5),colnames(netOutDCM))
namez <- unique(namez)
namez <- str_remove(namez,"S_")
dfCLR <- subset(dfCLR,select=c(namez))
namez <- colnames(dfCLR)

# Load environmental data 
env <- read.csv("SPOT_Env_NewJan11.csv",header=TRUE)

# Match up environmental data to order of ASV table samples (not necessarily ordered by time)
colz <- colsplit(rownames(dfCLR),"_",c("spot","Cruise","Month","Day","Year","Depth"))
dfCLR$Depth <- colz$Depth
dfCLR$Cruise <- colz$Cruise
dfCLR <- dfCLR %>% arrange(Cruise)
all <- left_join(dfCLR,env)
dfCLR$Depth <- NULL
dfCLR$Cruise <- NULL

# Now we can select just the environmental variables we care about/ variable that do not display multicollinearity
set.seed(100)
keepEnv <- c("DayOfYear","O2Wink","NH4","SiO3","PO4","NO2.NO3","CTDTMP",    "CTDBEAM","CTDFLUOR","CTDOXY","CTDSAL","DayLength","DayDiff","SLA","cyanos","sar11","cyanos_PA","sar11_PA","MEI","SST","Chla","PP")
envNew <- subset(all,select=c(keepEnv))
envNew <- missForest(envNew)
envNew <- envNew$ximp

# Rename CLR transformed dataframe
total <- dfCLR
total <- cbind(total,envNew)
envNamez <- colnames(envNew)
`%ni%` <- Negate(`%in%`)

# Set up training and test dataset for ASV-ASV predictions
finalTest <- total[181:242,] # 25% of data used for test
finalTrain <- total[1:180,] # 75% of data used for train

out <- NULL
coefs <- data.frame(a=c(1:ncol(total)))
for (i in 1:ncol(dfCLR)){
  y_train <- finalTrain[, i] # Target ASV
  X_train <- finalTrain[, -i] # Everything but target ASV
  y_test <- finalTest[,i] # Target ASV test
  X_test <- finalTest[,-i] # Everything but target ASV test
  
  X_train <- as.matrix(X_train)
  f <- rep(1:15,each=12)
  
  l <- seq(0,1,by=0.01)
  
  cvOut <- cv.glmnet(X_train,y_train,nfolds=15,foldid=f,keep=TRUE,type.measure = "mse",lambda=c(l),standardize=TRUE)
  
  coefOut <- coef(cvOut,s=cvOut$lambda.min)
  coefOut <- as.matrix(coefOut)
  coefOut <- as.data.frame(coefOut)
  coefs <- cbind(coefs,coefOut)
  coefOut <- subset(coefOut,s1!=0 & rownames(coefOut)!="(Intercept)")
  
  envNum <- sum(rownames(coefOut) %in% envNamez)
  asvNum <- sum(rownames(coefOut) %ni% envNamez)
  
  X_test <- as.matrix(X_test)
  outPred <- predict(cvOut$glmnet.fit,newx=X_test)
  mseOut <- mse(y_test,outPred)
  
  vec <- c(colnames(total)[i],envNum,asvNum,mseOut)
  out <- rbind(out,vec)}

coefs$a <- NULL
colnames(coefs) <- colnames(dfCLR)


out <- as.data.frame(out)
rownames(out)<- out$V1
out$V1 <- NULL

out$V2 <- as.numeric(out$V2)
out$V3 <- as.numeric(out$V3)
out$V4 <- as.numeric(out$V4)
colnames(out) <- c("Env","Asv","MSE")

out$percEnv <- out$Env/ncol(envNew)
out$percAsv <- out$Asv/(ncol(dfCLR)-1)
mean(out$percEnv)
mean(out$percAsv)
median(out$percEnv)
median(out$percAsv)

percComp <- data.frame(num=c(rownames(out)),env=c(out$percEnv),asv=c(out$percAsv))
percComp <- melt(percComp,id.vars=c("num"))

# Boxplots
p2 <- ggplot(percComp,aes(x=variable,y=value))+geom_boxplot()+theme_bw()+xlab("Predictor Group")+ylab("Percent of Total Predictors per Model")+scale_x_discrete(labels=c("Environmental/Biological Parameters\n(22 Total)","ASVs\n(1310 Total)"))
p1+p2

# ggsave("../../TRY.pdf",width=12,height=5)


fin <- read.delim("taxonomy_90.tsv",header=TRUE)
out$Feature.ID <- rownames(out)
out <- left_join(out,fin)

taxz <- colsplit(out$Taxon,";",c("Supergroup","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxz$fin <- ifelse(taxz$Class=="Syndiniales","Syndiniales",NA)
taxz$fin <- ifelse(taxz$Class=="Dinophyceae","Dinoflagellate",taxz$fin)
taxz$fin <- ifelse(taxz$Class=="Bacillariophyta","Diatom",taxz$fin)
taxz$fin <- ifelse(grepl("MAST",taxz$Class),"MAST",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Rhizaria","Rhizaria",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Chlorophyta","Chlorophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Cryptophyta","Cryptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Haptophyta","Haptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Ciliophora","Ciliate",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Stramenopiles" & is.na(taxz$fin),"Other Stramenopiles",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Archaeplastida" & is.na(taxz$fin),"Other Archaeplastids",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Alveolata" & is.na(taxz$fin),"Other Alveolate",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin),"Other Eukaryote",taxz$fin)
unique(taxz$fin)

out$fin <- taxz$fin

# Histograms
hist1 <- ggplot(out, aes(x=Env)) + geom_histogram(aes(y=..density..), binwidth=1,colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666",adjust = 2)+theme_classic()+xlab("Number of Environmental/Biological Predictors Per Model")+ylab("Frequency")

hist2 <- ggplot(out, aes(x=Asv)) + geom_histogram(aes(y=..density..), binwidth=20,colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666",adjust=2)+theme_classic()+xlab("Number of ASV Predictors Per Model")+ylab("Frequency")
hist1 +hist2
# ggsave("../../NEWESTTOP.pdf",width=10,height=4)

# Look at important predictors
envCoef <- subset(coefs,rownames(coefs) %in% envNamez)
envCoef$var <- rownames(envCoef)
envCoefM <- melt(envCoef,id.vars="var") 
sumz <- envCoefM %>% group_by(var) %>% tally(value!=0) %>% arrange(desc(n))

ggplot(sumz,aes(x=reorder(var,-n),y=n))+geom_bar(stat="identity",color="black",fill="grey")+theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+xlab("Biological/Environmental Predictor")+ylab("# Models Containing Biological/Environmental Predictor")+scale_x_discrete(labels=c("Day length","Chlorophyll a (satellite)","Oxygen","Silica","Cyanobacteria","Ammonium","Oxygen (Winkler)","Nitrate + Nitrite","Day of year","SSH (satellite)","Chlorophyll a","MEI","Particle-associated SAR11","SST (satellite)","Temperature","Primary production (satellite)","Beam attenuation","Salinity","Particle-associated cyanobacteria","Rate of change of day length","Phosphate","SAR11"))
# ggsave("../../NEW2.pdf",width=10,height=7)

dayLength <- subset(envCoefM,var=="DayLength")
dayLength <- subset(dayLength,value!=0)
colnames(dayLength)[2] <- "Feature.ID"
dayLength <- left_join(dayLength,fin)
