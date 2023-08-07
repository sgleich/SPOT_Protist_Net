### SPOT Network Analysis ###
### Figure 2: RDA ###
### By: Samantha Gleich ###
### Last Updated: 8/7/23 ###

# Load libraries
library(tidyverse)
library(reshape2)
library(ggplot2)
library(randomcoloR)
library(vegan)
library(compositions)
library(missForest)

# Set SPOT working directory
setwd("~/Desktop/SPOT/export_dir")

# Load in counts data and manifest file
counts <- read.delim("feature-table.tsv",header=TRUE,row.names=1)
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
counts$namez <- rownames(counts)
dfNew <- left_join(counts,fin)
rownames(dfNew) <- dfNew$new
dfNew$namez <- NULL
dfNew$new <- NULL

# Add taxonomy information
tax <- read.delim("taxonomy_90.tsv",header=TRUE)
tax <- tax[c(1:2)]
dfNew <- as.data.frame(t(dfNew))
dfNew$Feature.ID <- rownames(dfNew)
dfTax <- left_join(dfNew,tax)
dfTax <- subset(dfTax,grepl("Eukaryot",dfTax$Taxon))
dfTax$Feature.ID <- NULL
dfTax$Taxon <- NULL
dfTax <- as.data.frame(t(dfTax))

runRDA <- function(df,depth){
  df <- subset(df,grepl(paste(depth),rownames(df)))
  df <- subset(df,rownames(df)!="SPOT_115_2_16_12_5m"& rownames(df)!="SPOT_115_2_16_12_DCM")
  
  normDf <- as.data.frame(clr(df))
  env <- read.csv("../SPOT_Env_NewJan11.csv",header=TRUE)
  colz <- colsplit(rownames(normDf),"_",c("spot","Cruise","month","day","year","Depth"))
  colz <- left_join(colz,env)
  
  if (depth=="5m"){
    keepEnv <- c("DayOfYear","O2Wink","NH4","SiO3","PO4","NO2.NO3","CTDTMP","CTDBEAM","CTDFLUOR","CTDOXY","CTDSAL","DayLength","DayDiff","SLA","cyanos","sar11","cyanos_PA","sar11_PA","MEI","SST","Chla","PP")}
  
  if (depth=="DCM"){
    keepEnv <- c("DayOfYear","CSDepth","O2Wink","NH4","SiO3","PO4","NO2.NO3","CTDTMP","CTDBEAM","CTDFLUOR","CTDOXY","CTDSAL","DayLength","DayDiff","SLA","cyanos","sar11","cyanos_PA","sar11_PA","MEI","SST","Chla","PP")}
  
  envTable <- subset(colz,select=c(keepEnv))
  set.seed(100)
  envTable <- missForest::missForest(envTable)
  envTable <- data.frame(envTable$ximp)
  
  envTable <- decostand(envTable, method = "standardize")
  
  tryAgain=TRUE
  while(tryAgain==TRUE){
    rdaOut<- rda(normDf ~ ., data = envTable)
    rdaDf <- data.frame(rdaOut$CCA$u)
    rdaDf$Month <- colz$Month
    
    finalmodel<- ordistep(rdaOut, scope=formula(rdaOut))
    t <- vif.cca(finalmodel)
    t2 <- which(t>=10)
    tryAgain <- ifelse(length(t2)>=1,TRUE,FALSE)
    
    if(length(t2)>=1){
      t <- which.max(t)
      envTable <- envTable[,!names(envTable) %in% c(names(t))]}}
  
  
  rda.0 <- rda (normDf ~ 1, data = envTable) 
  rda.all<- rda(normDf ~ ., data = envTable)
  
  adjR2.tbrda <- RsquareAdj(rda.all)$adj.r.squared
  
  forSel <- ordiR2step(rda.0, scope = formula (rda.all), R2scope = adjR2.tbrda, direction = 'forward', permutations = 999)
  
  forSel$anova$`Pr(>F)` <- p.adjust (forSel$anova$`Pr(>F)`, method = 'holm', n = ncol (envTable))
  
  import <- data.frame(forSel$anova)
  
  sub <- rownames(forSel$anova)
  sub <- str_remove_all(sub,"\\+ ")
  len <- length(sub)
  sub <- sub[1:len-1]
  envTable2 <- subset(envTable,select=c(sub))
  
  finalModforSel <- rda (normDf ~ ., data = envTable2) 
  summar <- summary(finalModforSel)
  summar <- data.frame(summar$cont$importance)
  
  rdaDf <- data.frame(finalModforSel$CCA$u)
  
  fitDf <- data.frame(finalModforSel$CCA$biplot)
  fitDf <- fitDf[c(1:2)]
  fitDf$x <- mean(rdaDf$RDA1)
  fitDf$y <- mean(rdaDf$RDA2)
  
  return(list(rdaDf,fitDf,import,summar))}

outSurf <- runRDA(dfTax,"5m")
mainDf <- outSurf[[1]]
fitDf <- outSurf[[2]]
impDf <- outSurf[[3]]
statDf <- outSurf[[4]]

outDCM <- runRDA(dfTax,"DCM")
mainDf <- outDCM[[1]]
fitDf <- outDCM[[2]]
impDf <- outDCM[[3]]
statDf <- outDCM[[4]]

colz <- colsplit(rownames(mainDf),"_",c("spot","Cruise","month","day","year","Depth"))
mainDf$Month <- colz$month
statDf <- subset(statDf,select=c(grepl("RDA",colnames(statDf))))
rowSums(statDf)

outP <- ggplot(mainDf,aes(RDA1,RDA2,fill=factor(Month),shape=as.factor(Month)))+geom_point(size=2,color="black")+scale_fill_manual(name="Month",values=c("dodgerblue","dodgerblue","darkolivegreen4","darkolivegreen4","darkolivegreen4","goldenrod1","goldenrod1","goldenrod1","firebrick2","firebrick2","firebrick2","dodgerblue"),labels=c("January","February","March","April","May","June","July","August","September","October","November","December"))+scale_shape_manual(name="Month",values=c(21,22,21,22,24,21,22,24,21,22,24,24),labels=c("January","February","March","April","May","June","July","August","September","October","November","December"))+theme_classic()+xlab("RDA1 (3.35%)")+ylab("RDA2 (1.70%)")+geom_vline(xintercept = 0,linetype="dotted")+geom_hline(yintercept = 0,linetype="dotted")+ggtitle("DCM")+geom_segment(aes(x = fitDf[1,3], y = fitDf[1,4], xend = fitDf[1,1] , yend = fitDf[1,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[2,3], y = fitDf[2,4], xend = fitDf[2,1] , yend = fitDf[2,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[3,3], y = fitDf[3,4], xend = fitDf[3,1] , yend = fitDf[3,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[4,3], y = fitDf[4,4], xend = fitDf[4,1] , yend = fitDf[4,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[5,3], y = fitDf[5,4], xend = fitDf[5,1] , yend = fitDf[5,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)#+geom_segment(aes(x = fitDf[6,3], y = fitDf[6,4], xend = fitDf[6,1] , yend = fitDf[6,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)+geom_segment(aes(x = fitDf[7,3], y = fitDf[7,4], xend = fitDf[7,1] , yend = fitDf[7,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)#+geom_segment(aes(x = fitDf[8,3], y = fitDf[8,4], xend = fitDf[8,1] , yend = fitDf[8,2]),arrow = arrow(length=unit(3, "mm")),linewidth=0.5)
outP
ggsave("../../DCMRDA_Aug2023.pdf",width=6,height=4)

