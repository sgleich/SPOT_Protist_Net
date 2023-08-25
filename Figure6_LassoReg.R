### SPOT LASSO ANALYSIS ###
### Script by: Samantha Gleich ###
### Figure 6: Lasso Reg ###
### Last modified: 8/25/23 ###

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

# Tax color palette
taxCols <- c("#E1746D","#76C3D7","#DE5AB1","#D5E0AF","#DED3DC","#87EB58","#D4DC60","#88E5D3","#88AAE1","#DBA85C","#8B7DDA","#9A8D87","#D99CD1","#B649E3","#7EDD90","#4FC4D0")
names(taxCols) <- c("Chlorophyte","Ciliate","Cryptophyte","Diatom","Haptophyte","Dinoflagellate","MAST","Other Alveolates","Other Archaeplastids","Other Eukaryote","Other Stramenopiles","Rhizaria","Group I Syndiniales","Group II Syndiniales","Unknown Eukaryote","Metazoa")

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
netOut5 <- read.csv("../../Surface_SPOT_Aug2023.csv",header=TRUE,row.names=1) 
netOutDCM <- read.csv("../../DCM_SPOT_Aug2023.csv",header=TRUE,row.names=1)
namez <- c(colnames(netOut5),colnames(netOutDCM))
namez <- unique(namez)
namez <- str_remove(namez,"S_")
dfCLR <- subset(dfCLR,select=c(namez))
namez <- colnames(dfCLR)

# Load environmental data 
env <- read.csv("env_impute_jun2.csv",header=TRUE)
env2 <- read.csv("SPOT_Env_NewJan11.csv",header=TRUE)
env$Cruise <- env2$Cruise
env$Depth <- env2$Depth

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
colnames(envNew) <- c("Day of year","Oxygen (Winkler)","Ammonium","Silica","Phosphate","Nitrate + Nitrite","Temperature (CTD)","Beam attenuation (CTD)","Chlorophyll a fluorescence (CTD)","Oxygen (CTD)","Salinity (CTD)","Day length","Rate of change of day length","Sea surface height","Cyanobacteria","SAR 11","Particle-associated cyanobacteria","Particle-associated SAR11","MEI","Sea surface temperature (satellite)","Chlorophyll a fluorescence (satellite)","Primary production (satellite)")

# Rename CLR transformed dataframe
total <- dfCLR
total <- cbind(total,envNew)
envNamez <- colnames(envNew)
`%ni%` <- Negate(`%in%`)

# Set up training and test dataset for ASV-ASV predictions
total <- decostand(total,method='standardize')
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
  set.seed(100)
  cvOut <- cv.glmnet(X_train,y_train,nfolds=20,foldid=f,keep=TRUE,type.measure = "mse",lambda=c(l),standardize=FALSE)
  
  coefOut <- coef(cvOut,s=cvOut$lambda.min)
  coefOut <- as.matrix(coefOut)
  coefOut <- as.data.frame(coefOut)
  coefs <- cbind(coefs,coefOut)
  coefOut <- subset(coefOut,s1!=0 & rownames(coefOut)!="(Intercept)")
  
  envNum <- sum(rownames(coefOut) %in% envNamez)
  asvNum <- sum(rownames(coefOut) %ni% envNamez)
  
  X_test <- as.matrix(X_test)
  if(length(unique(y_test))>1){
    finalMod <- glmnet(X_test,y_test,type.measure = "mse",lambda=cvOut$lambda.min,standardize=FALSE)
    rSq <- finalMod$dev.ratio
  # outPredTest <- predict(cvOut$glmnet.fit,newx=X_test)
  # cvOut$glmnet.fit$dev.ratio
  
  # outPredTrain <- predict(cvOut$glmnet.fit,newx=X_train)
  
  # mseOutTrain <- Metrics::rmse(y_train,outPredTrain)
  # mseOutTest <- Metrics::rmse(y_test,outPredTest)
    vec <- c(colnames(total)[i],envNum,asvNum,rSq)
    out <- rbind(out,vec)
    print(i)}}

# Reorganize output - only keep models with R^2 > 0.7
coefs$a <- NULL
colnames(coefs) <- colnames(dfCLR)
coefs <- coefs[2:nrow(coefs),]
out <- as.data.frame(out)
rownames(out)<- out$V1
out$V1 <- NULL
out$V2 <- as.numeric(out$V2)
out$V3 <- as.numeric(out$V3)
out$V4 <- as.numeric(out$V4)
colnames(out) <- c("Env","Asv","Rsq")
out2 <- out
out <- subset(out,Rsq>0.7)
# out <- subset(out,Asv!=0)

out$percEnv <- out$Env/ncol(envNew)
out$percAsv <- out$Asv/(ncol(dfCLR)-1)
mean(out$percEnv)
mean(out$percAsv)
median(out$percEnv)
median(out$percAsv)

percComp <- data.frame(num=c(rownames(out)),env=c(out$percEnv),asv=c(out$percAsv))
percComp <- melt(percComp,id.vars=c("num"))

# Boxplots
fun_mean <- function(x){
  return(data.frame(y=mean(x),label=scales::percent(round(mean(x,na.rm=T),4),0.01)))}

p <- ggplot(percComp,aes(x=variable,y=value))+geom_boxplot(outlier.shape=21)+theme_bw(base_size=14)+xlab("Predictor Group")+ylab("Percent of Total Predictors per Model")+scale_x_discrete(labels=c("Environmental/Biological Parameters\n(22 Total)","ASVs\n(1309 Total)"))+stat_summary(fun.data = fun_mean, geom="text", hjust=1.5,vjust=-2,color="red")+ scale_y_continuous(labels = scales::percent)
p

# Histogram
out$total <- out$Env+out$Asv
hist <- ggplot(out, aes(x=total)) + geom_histogram(aes(y=..density..), binwidth=20,colour="black", fill="white") +
  geom_density(alpha=.2, fill="grey30",adjust=2)+theme_classic(base_size=14)+xlab("Number of Predictors Per Model")+ylab("Frequency")
hist


# Look at important environmental predictors
fin <- read.delim("taxonomy_90.tsv",header=TRUE)
coefs <- subset(coefs,select=c(namez))
envCoef <- subset(coefs,rownames(coefs) %in% envNamez)
asvCoef <- subset(coefs,rownames(coefs) %in% namez)
envCoef$var <- rownames(envCoef)
asvCoef$var <- rownames(asvCoef)
envCoefM <- melt(envCoef,id.vars="var") 
asvCoefM <- melt(asvCoef,id.vars="var") 
colnames(envCoefM)[2] <- "Feature.ID"
envCoefM <- left_join(envCoefM,fin)
taxz <- colsplit(envCoefM$Taxon,";",c("Supergroup","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxz$fin <- ifelse(taxz$Class=="Dinophyceae","Dinoflagellate",NA)
taxz$fin <- ifelse(taxz$Order=="Dino-Group-I","Group I Syndiniales",taxz$fin)
taxz$fin <- ifelse(taxz$Order=="Dino-Group-II","Group II Syndiniales",taxz$fin)
taxz$fin <- ifelse(taxz$Class=="Bacillariophyta","Diatom",taxz$fin)
taxz$fin <- ifelse(grepl("MAST",taxz$Class),"MAST",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Rhizaria","Rhizaria",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Chlorophyta","Chlorophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Cryptophyta","Cryptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Haptophyta","Haptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Ciliophora","Ciliate",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Metazoa","Metazoa",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Stramenopiles" & is.na(taxz$fin),"Other Stramenopiles",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Archaeplastida" & is.na(taxz$fin),"Other Archaeplastids",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Alveolata" & is.na(taxz$fin),"Other Alveolate",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin),"Other Eukaryote",taxz$fin)
envCoefM$fin <- taxz$fin


sumz <- envCoefM %>% group_by(var,fin) %>% tally(value!=0) %>% arrange(desc(n))

env <- ggplot(sumz,aes(x=reorder(var,-n),y=n,fill=fin))+geom_bar(stat="identity",color="black")+theme_classic(base_size=14)+theme(axis.text.x = element_text(angle = 45,hjust=1))+xlab("Environmental/Biological Predictor")+ylab("Number of Models Containing\nEnvironmental/Biological Predictor")+scale_fill_manual(name=c("Taxonomic Group"),values=c(taxCols))
# ggsave("../../NEW2.pdf",width=10,height=7)

(p+hist)/env+plot_annotation(tag_levels="a")
ggsave("../../Figure6_Aug2023.pdf",width=16,height=12)
