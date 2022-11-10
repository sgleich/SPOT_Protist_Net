### SPOT Network Analysis ###
### Pre-processing (Filtering + GAM-transforming) ###
### By: Samantha Gleich ###
### Last Updated: 11/10/22 ###

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
fin <- data.frame(namez=c(namez$V1),new=c(namez$fin))
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
df5Filt <- subset(df5,rowSums(df5==0) <= 25) # 25 == ~20% not zero
# Optional: Check to see how many 0s are in the dataframe for each ASV
# df5Filt$count <- rowSums(df5Filt==0) 
# df5Filt$count <- NULL

dfDCMFilt <- subset(dfDCM,rowSums(dfDCM==0) <= 25) # 25 == ~20% not zero
# Optional: Check to see how many 0s are in the dataframe for each ASV
dfDCMFilt$count <- rowSums(dfDCMFilt==0)
dfDCMFilt$count <- NULL

# For this analysis, we will include ASVs in both networks that passed our filtering threshold in either network. 
namez5 <- rownames(df5Filt)
namezDCM <- rownames(dfDCMFilt)
namez <- c(namez5,namezDCM)
namez <- unique(namez) # These are all of the ASVs that will be included in both networks

# Now let's grab the ASVs that will be included in our networks from the CLR-transformed dataframe.
dfCLR <- as.data.frame(dfCLR)
dfCLRFilt <- subset(dfCLR,select=c(namez))

# Now we need to step up some vectors for our NetGAM time-series transformation (vectos are MOY and DayofTS -- see NetGAM documentation)
colz <- colsplit(rownames(dfclrFilt),"_",c("spot","num","month","day","year","depth"))
vec <- rep(1:12, length=192)
vec <- as.data.frame(vec)
vec$year <- rep(3:18, each=12)
vec$DayOTS <- 1:192
colnames(vec)<- c("month","year","day")
colz <- data.frame(colz$month,colz$year)
colnames(colz)<- c("month","year")
colz <- left_join(colz,vec)

# Add month of year and day of time-series information to our CLR-transformed dataframe
dfCLRFilt$month <- colz$month
dfCLRFilt$day <- colz$day
df5CLR <- subset(dfCLRFilt,grepl("5m", rownames(dfCLRFilt)))
dfDCMCLR<- subset(dfCLRFilt,grepl("DCM", rownames(dfCLRFilt)))

# Remove duplicate samples.
df5CLR <- subset(df5CLR,rownames(df5CLR)!="SPOT_115_2_16_12_5m")
dfDCMCLR <- subset(dfDCMCLR,rownames(dfDCMCLR)!="SPOT_115_2_16_12_DCM")

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

# Save dataframes that will be used for eLSA network runs
netGAM5<- as.data.frame(t(netGAM5))
netGAMDCM <- as.data.frame(t(netGAMDCM))
write.csv(netGAM5,"SPOT_5m_Filtered_GAM_Nov2022.csv")
write.csv(netGAMDCM,"SPOT_DCM_Filtered_GAM_Nov2022.csv")

### Run eLSA ###
# Command to run eLSA on server - networks for surface and DCM were run separately. 
# lsa_compute SPOT_5m_Filtered_GAM_Nov2022.txt SPOT_5m_Filtered_GAM_Nov2022_out.txt -p perm -r 1 -s 122 -d 1 -n none 
