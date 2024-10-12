# Load in PIDA
pida <- read.csv("PIDA_Int.csv",header=TRUE) 
pidaPred <- subset(pida,Taxonomic.interaction=="Prot - Prot" & Interaction=="par")

# Load in the counts dataframe from qiime2 and the manifest used to run the qiime2 pipeline
counts <- read.delim("feature-table.tsv",header=TRUE,row.names=1,skip=1)
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
dfNew <- subset(dfNew,rownames(dfNew)!="SPOT_115_2_16_12_5m"& rownames(dfNew)!="SPOT_115_2_16_12_DCM")
df5 <- subset(dfNew,grepl("5m", rownames(dfNew)))
dfDCM<- subset(dfNew,grepl("DCM", rownames(dfNew)))
df5 <- as.data.frame(t(df5))
dfDCM <- as.data.frame(t(dfDCM))

# This is our ASV filtration step. If we include all ~30,000 ASVs in the network analysis it would take too long to run and it would be difficult to interpret the output. Instead, we will filter our surface (5 m) and DCM datasets so that only ASVs that are non-zero in > 20% of samples are included in the network analysis
df5Filt <- subset(df5,rowSums(df5==0) <= 98) # ~20% of samples have to be non-zero
# Optional: Check to see how many 0s are in the dataframe for each ASV
# df5Filt$count <- rowSums(df5Filt==0) 
# df5Filt$count <- NULL

dfDCMFilt <- subset(dfDCM,rowSums(dfDCM==0) <= 96) # ~20 % of samples have to be non-zero
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

dfCLRFilt <- subset(dfCLRFilt,rownames(dfCLRFilt)!="SPOT_115_2_16_12_5m"& rownames(dfCLRFilt)!="SPOT_115_2_16_12_DCM")

# Get ASVs for network
df5CLR <- subset(dfCLRFilt,grepl("5m", rownames(dfCLRFilt)))
dfDCMCLR<- subset(dfCLRFilt,grepl("DCM", rownames(dfCLRFilt)))

namez <- colnames(df5CLR)
namez <- data.frame(Feature.ID=namez)
tax <- read.delim("taxonomy.tsv",header=TRUE)
namez <- left_join(namez,tax)

# Compare
cols <- colsplit(namez$Taxon,";",c("d","k","p","c","o","f","g","g2","sp"))
namez$G <- cols$g2

# Symb
sub <- subset(pidaPred,Genus.org1 %in% namez$G | grepl("Acantharea",pidaPred$Taxonomic.level.3..org1) | grepl("MAST",pidaPred$Taxonomic.level.3..org1) | grepl("Polycyst",pidaPred$Taxonomic.level.3..org1))
sub <- subset(sub,Genus.org2 %in% namez$G | grepl("Acantharea",sub$Taxonomic.level.3..org2) | grepl("MAST",sub$Taxonomic.level.3..org2) | grepl("Polycyst",sub$Taxonomic.level.3..org2))

# Pred
sub <- subset(pidaPred, Genus.org1 %in% namez$G)
sub <- subset(sub,Genus.org2 %in% namez$G)

# Para
sub <- subset(pidaPred, Genus.org1 %in% namez$G | grepl("Syndiniales",pidaPred$Genus.org1))
sub <- subset(sub,Genus.org2 %in% namez$G | grepl("Syndiniales",sub$Genus.org2))
