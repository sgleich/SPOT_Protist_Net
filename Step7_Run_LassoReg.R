### SPOT LASSO ANALYSIS ###
### Script by: Samantha Gleich ###
### Last modified: April 4, 2023 ###

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
netOut <- read.csv(file.choose(),header=TRUE,row.names=1) 
namez <- c(netOut$V1,netOut$V2)
namez <- unique(namez)
namez <- str_remove(namez,"S_")
dfCLR <- subset(dfCLR,select=c(namez))
namez <- colnames(dfCLR)

# Load environmental data 
env <- read.csv("../SPOT_Env_NewJan11.csv",header=TRUE)

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
envNew <- all[c(1030,1032:1052)]
envNew <- missForest(envNew)
envNew <- envNew$ximp
envNew$DayDiff <- NULL
envNew$cyanos_PA <- NULL
envNew$CTDFLUOR <- NULL

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
out <- subset(out, rowSums(out)!=0)
colnames(out) <- c("Env","Asv","MSE")

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

# Look at relative number of environmental vs. ASV predictors
out$ratEnv <- out$Env/19
out$ratAsv <- out$Asv/1024
out$mod <- 1:nrow(out)


newDf <- data.frame(envRat =c(out$ratEnv),asvRat=c(out$ratAsv),mod=c(out$mod))
newDf <- newDf %>% pivot_longer(cols = c("envRat","asvRat"))

rat <- ggplot(newDf,aes(x=name,y=value))+geom_boxplot()+ylim(0,0.4)+theme_classic()+ylab("# of Predictors in Final Model/\nTotal # of Input Predictors")+scale_x_discrete(labels=c("ASVs","Environmental/\nBiological Variables"))+xlab("Predictor Group")


newDf2 <- data.frame(env =c(out$Env),asv=c(out$Asv),mod=c(out$mod))
newDf2 <- newDf2 %>% pivot_longer(cols = c("env","asv"))

reg <- ggplot(newDf2,aes(x=name,y=value))+geom_boxplot()+theme_classic()+ylab("# of Predictors in Final Model")+scale_x_discrete(labels=c("ASVs","Environmental/\nBiological Variables"))+xlab("Predictor Group")+ylim(0,200)

reg+rat
ggsave("../../NEWLassoReg.pdf")

# Look at important predictors
envCoef <- subset(coefs,rownames(coefs) %in% envNamez)
envCoef$var <- rownames(envCoef)
envCoefM <- melt(envCoef,id.vars="var") 
sumz <- envCoefM %>% group_by(var) %>% tally(value!=0) %>% arrange(desc(n))

ggplot(sumz,aes(x=reorder(var,-n),y=n))+geom_bar(stat="identity",color="black",fill="grey")+theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+xlab("Biological/Environmental Predictor")+ylab("# ASV Models Containing Biological/Environmental Predictor")+scale_x_discrete(labels=c("Day length","Oxygen","Silica","Chlorophyll a (satellite)","Cyanobacteria","Day of year","Oxygen (Winkler)","MEI","Ammonium","Nitrate","Primary production (satellite)","SSH (satellite)","SST (satellite)","Beam attenuation","SAR11","Salinity","Temperature","Particle-associated SAR11","Phosphate"))
ggsave("../../NEW2.pdf",width=10,height=7)

out2 <- out %>% arrange(desc(total))
out2 <- subset(out2,Asv!=1023)
out2 <- subset(out2,Asv!=0)
top10 <- out2 %>% top_frac(0.1,Asv)
bottom10 <- out2 %>% top_frac(-0.1,Asv)
median(top10$Env)
median(bottom10$Env)
wilcox.test(top10$Env,bottom10$Env)
full <- rbind(top10,bottom10)


hist <- out %>% filter (total < 1000 & total > 0) %>% ggplot(aes(x=total)) + geom_histogram(aes(y=..density..), binwidth=10,colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666",adjust=2)+theme_classic()+xlab("Number of Predictors Per Model")+ylab("Frequency")
hist
ggsave("Histogram.pdf")
