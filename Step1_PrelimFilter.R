### SPOT Network Analysis ###
### Pre-processing (Filtering + GAM-transforming) ###
### By: Samantha Gleich ###
### Last Updated: 12/14/22 ###

# Load libraries
library(tidyverse)
library(compositions)
library(NetGAM)
library(mgcv)
library(reshape2)

# Load in the counts dataframe from qiime2 and the manifest used to run the qiime2 pipeline
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

# Now let's apply a clr transformation to our dataframe
dfCLR <- clr(dfNew)

# Let's separate the non-clr transformed samples by depth 
df5 <- subset(dfNew,grepl("5m", rownames(dfNew)))
dfDCM<- subset(dfNew,grepl("DCM", rownames(dfNew)))
df5 <- as.data.frame(t(df5))
dfDCM <- as.data.frame(t(dfDCM))

# This is our ASV filtration step. If we include all ~30,000 ASVs in the network analysis it would take too long to run and it would be difficult to interpret the output. Instead, we will filter our surface (5 m) and DCM datasets so that only ASVs that are non-zero in > 20% of samples are included in the network analysis
df5Filt <- subset(df5,rowSums(df5==0) <= 62) # ~60% of samples have to be non-zero
# Optional: Check to see how many 0s are in the dataframe for each ASV
# df5Filt$count <- rowSums(df5Filt==0) 
# df5Filt$count <- NULL

dfDCMFilt <- subset(dfDCM,rowSums(dfDCM==0) <= 61) # ~60 % of samples have to be non-zero
# Optional: Check to see how many 0s are in the dataframe for each ASV
# dfDCMFilt$count <- rowSums(dfDCMFilt==0)
# dfDCMFilt$count <- NULL

# For this analysis, we will include ASVs in both networks that passed our filtering threshold in either network. 
namez5 <- rownames(df5Filt)
namezDCM <- rownames(dfDCMFilt)
namez <- c(namez5,namezDCM)
namez <- unique(namez) # These are all of the ASVs that will be included in both networks

# Now let's grab the ASVs that will be included in our networks from the CLR-transformed dataframe.
dfCLR <- as.data.frame(dfCLR)
dfCLRFilt <- subset(dfCLR,select=c(namez))

# Now we need to step up some vectors for our NetGAM time-series transformation (vectors are MOY and DayofTS -- see NetGAM documentation)
vec <- rep(1:12, length=192)
vec <- as.data.frame(vec)
vec$year <- rep(3:18, each=12)
vec$DayOTS <- 1:192
colnames(vec)<- c("month","year","day")
vec <- vec[9:nrow(vec),]
vec$day <- 1:nrow(vec)

env <- read.csv("../SPOT_NEW.csv",header=TRUE)
env$year <- stri_sub(env$Date, -2,-1)
env$year <- as.numeric(as.character(env$year))
colz <- colsplit(env$Date,"/",c("m","d","y"))
env <- data.frame(month=c(env$Month),year=c(env$year),Cruise=c(env$Cruise),Depth=c(env$Depth))
env <- left_join(env,vec)

dfCLRFilt <- subset(dfCLRFilt,rownames(dfCLRFilt)!="SPOT_115_2_16_12_5m"& rownames(dfCLRFilt)!="SPOT_115_2_16_12_DCM")

# Add month of year and day of time-series information to our CLR-transformed dataframe
colz <- colsplit(rownames(dfCLRFilt),"_",c("SPOT","Cruise","Month","Day","Year","Depth"))
dfCLRFilt$Cruise <- colz$Cruise
dfCLRFilt$Depth <- colz$Depth
env$year <- NULL
namez <- rownames(dfCLRFilt)
dfCLRFilt <- left_join(dfCLRFilt,env)
dfCLRFilt$Depth <- NULL
dfCLRFilt$Cruise <- NULL
rownames(dfCLRFilt) <- namez

df5CLR <- subset(dfCLRFilt,grepl("5m", rownames(dfCLRFilt)))
dfDCMCLR<- subset(dfCLRFilt,grepl("DCM", rownames(dfCLRFilt)))

# NetGAM expects month of year and day of time-series vectors (not columns). Set up vectors for 5m samples. 
df5Month <- df5CLR$month
df5Day <- df5CLR$day
df5CLR$month <- NULL
df5CLR$day <- NULL

# NetGAM also doesn't like column names that start with numbers. So, we'll add an "S" to all column names.
namez <- colnames(df5CLR)
namez <- paste("S",namez,sep="_")
colnames(df5CLR) <- namez

# Run NetGAM for 5m samples
netGAM5 <- netGAM.df(df5CLR,MOY=df5Month,MCount=df5Day,clrt=FALSE)

# NetGAM expects month of year and day of time-series vectors (not columns). Set up vectors for DCM samples.
dfDCMMonth <- dfDCMCLR$month
dfDCMDay <- dfDCMCLR$day
dfDCMCLR$month <- NULL
dfDCMCLR$day <- NULL

# NetGAM also doesn't like column names that start with numbers. So, we'll add an "S" to all column names.
namez <- colnames(dfDCMCLR)
namez <- paste("S",namez,sep="_")
colnames(dfDCMCLR) <- namez

# Run NetGAM for DCM samples
netGAMDCM <- netGAM.df(dfDCMCLR,MOY=dfDCMMonth,MCount=dfDCMDay,clrt=FALSE)

# Make transformed dataframes into matricies
netGAM5 <- as.matrix(netGAM5)
netGAMDCM <- as.matrix(netGAMDCM)

# Run graphical lasso network
lams  <- pulsar::getLamPath(pulsar::getMaxCov(netGAMDCM), .01, len=30)
hugeargs <- list(lambda=lams, verbose=FALSE)
out5 <- pulsar::batch.pulsar(netGAMDCM, fun=huge::huge, fargs=hugeargs,rep.num=50, criterion = "stars")
opt <- out5$stars
n <- opt$opt.index
# Get output adjacency matrix from graphical lasso model
fit <- pulsar::refit(out5)
fit2 <- fit$refit
fit.fin <- fit2$stars
fit.fin <- as.matrix(fit.fin)
fit.fin <- as.data.frame(fit.fin)
colnames(fit.fin) <- colnames(netGAM5)
rownames(fit.fin)<- colnames(netGAM5)
fit.fin <- as.matrix(fit.fin)

dim(fit.fin)

g <- graph_from_adjacency_matrix(fit.fin,mode="undirected")
plot(g,vertex.label=NA,vertex.size=6)
length(V(g))

write.csv(fit.fin,"Glasso_DCM_SPOT.csv")


### Run eLSA ###
# Command to run eLSA on server - networks for surface and DCM were run separately. 
# lsa_compute SPOT_5m_Filtered_GAM_Nov2022.txt SPOT_5m_Filtered_GAM_Nov2022_out.txt -p perm -r 1 -s 122 -d 1 -n none 
