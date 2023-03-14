### SPOT LASSO ANALYSIS ###
### Script by: Samantha Gleich ###
### Last modified: March 14, 2023 ###

# Load libraries
library(glmnet)
library(missForest)
library(tidyverse)
library(reshape2)
library(caret)
library(yardstick)
library(compositions)
library(Metrics)

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
netOut <- read.csv("Glasso_5m_SPOT_Dec28.csv",header=TRUE,row.names=1) 
namez <- colnames(netOut)
namez <- str_remove(namez,"S_")
dfCLR <- subset(dfCLR,select=c(namez))
namez <- colnames(dfCLR)
# s <- c(paste("V",1:389,sep="_"))
# colnames(dfCLR) <- s

# Load environmental data 
env <- read.csv("../SPOT_Env_NewJan11.csv",header=TRUE)

# Match up environmental data to order of ASV table samples (not necessarily ordered by time)
colz <- colsplit(rownames(dfCLR),"_",c("spot","Cruise","Month","Day","Year","Depth"))
dfCLR$Depth <- colz$Depth
dfCLR$Cruise <- colz$Cruise
all <- left_join(dfCLR,env)
dfCLR$Depth <- NULL
dfCLR$Cruise <- NULL

# Now we can select just the environmental variables we care about
set.seed(100)
envNew <- all[c(395,397:417)]
envNew <- missForest(envNew)
envNew <- envNew$ximp

# Rename CLR transformed dataframe
total <- dfCLR
#total <- mutate_all(total, function(x) as.numeric(as.character(x)))
total <- cbind(total,envNew)
envNamez <- colnames(envNew)
`%ni%` <- Negate(`%in%`)

# Set up training and test dataset for ASV-ASV predictions
finalTest <- total[181:242,] # 25% of data used for test
finalTrain <- total[1:180,] # 75% of data used for train

out <- NULL
coefs <- data.frame(a=c(1:411))
for (i in 1:389){
  y_train <- finalTrain[, i] # Target ASV
  X_train <- finalTrain[, -i] # Everything but target ASV
  y_test <- finalTest[,i] # Target ASV test
  X_test <- finalTest[,-i] # Everything but target ASV test
  
  X_train <- as.matrix(X_train)
  f <- rep(1:15,each=12)
  
  l <- seq(0,1,by=0.01)
  
  cvOut <- cv.glmnet(X_train,y_train,nfolds=15,foldid=f,standardize=FALSE,keep=TRUE,type.measure = "mse",lambda=c(l))
  
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
out <- subset(out, rowSums(out)!=0)
colnames(out) <- c("Env","Asv","MSE")

out$win <- ifelse(out$Env>out$Asv,"env","asv")
out$win <- ifelse(out$Env==out$Asv,"tie",out$win)
sum(out$win=="env")
sum(out$win=="asv")
sum(out$win=="tie")

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
hist1 <- ggplot(out,aes(x=Env))+geom_histogram(color="black",fill="grey")+theme_classic()+xlab("# Significant Environmental Predictors Per ASV")+ylab("Frequency")+ggtitle("Environmental Predictors")+scale_x_continuous(breaks = seq(0, 10, by = 1))

hist1b <- ggplot(out, aes(x=Env)) + geom_histogram(aes(y=..density..), binwidth=1,colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")+theme_classic()+xlab("# Significant Non-ASV Predictors Per ASV")+ylab("Frequency")

hist2 <- ggplot(out,aes(x=Asv))+geom_histogram(color="black",fill="grey")+theme_classic()+xlab("# Significant ASV Predictors Per ASV")+ylab("Frequency")+ggtitle("ASV Predictors")+scale_x_continuous(breaks = seq(0, 130, by = 20))

hist2b <- ggplot(out, aes(x=Asv)) + geom_histogram(aes(y=..density..), binwidth=6,colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")+theme_classic()+xlab("# Significant ASV Predictors Per ASV")+ylab("Frequency")

hist1b +hist2b
ggsave("../../glmnetOut_MARCH2023.pdf",width=10,height=4)

# Plot taxon-specific trends
out$rat <- out$Asv/out$Env
countz <- out %>% group_by(fin) %>%tally()
out <- left_join(out,countz)
asvP <- ggplot(out,aes(x=as.factor(fin),y=Asv))+geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust=1))+xlab("Taxonomic Group")+ylab("# Significant ASV Predictors")
envP <- ggplot(out,aes(x=as.factor(fin),y=Env))+geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust=1))+xlab("Taxonomic Group")+ylab("# Significant Environmental Predictors")
mseP <- ggplot(out,aes(x=as.factor(fin),y=MSE))+geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust=1))+xlab("Taxonomic Group")+ylab("Model MSE")
mseP

envP+asvP+mseP
ggsave("../../glmnetOut2.pdf",width=12,height=4)

## Look at important predictors
envCoef <- subset(coefs,rownames(coefs) %in% envNamez)
envCoef$var <- rownames(envCoef)
envCoefM <- melt(envCoef,id.vars="var") 
envCoefM <- subset(envCoefM, value!=0)

sumz <- envCoefM %>% group_by(var) %>% tally(value!=0) %>% arrange(desc(n))

ggplot(sumz,aes(x=reorder(var,-n),y=n))+geom_bar(stat="identity",color="black",fill="grey")+theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+xlab("Environmental Parameter")+ylab("# ASV Models Containing Environmental Predictor")+scale_x_discrete(labels=c("Primary production (Satellite)","Day of year","Day length","Nitrate + Nitrite","Chlorophyll a (Satellite)","Temperature","SST (Satellite)","Oxygen","Chlorophyll a fluorescence","MEI","Silica","Rate of change of day length","Oxygen (Winkler)","Salinity","Phosphate","Beam attenuation","Relative abundance cyanobacteria","Relative abundance particle-associated cyanobacteria","Ammonium","Relative abundance SAR11","Relative abundance particle-associated SAR11","SSH"))

ggsave("../../glmnetOut3.pdf",width=10,height=7)

