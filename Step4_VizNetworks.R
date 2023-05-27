### SPOT Network Analysis ###
### Visualize Networks + Calculate Network Statistics ###
### By: Samantha Gleich ###
### Last Updated: 5/27/2023 ###

taxCols <- c("#E1746D","#76C3D7","#DE5AB1","#D5E0AF","#DED3DC","#87EB58","#D4DC60","#88E5D3","#88AAE1","#DBA85C","#8B7DDA","#9A8D87","#D99CD1","#B649E3","#7EDD90")
names(taxCols) <- c("Chlorophyte","Ciliate","Cryptophyte","Diatom","Dinoflagellate","Fungi","Haptophyte","MAST","Other Alveolate","Other Archaeplastida","Other Eukaryote","Other Stramenopiles","Rhizaria","Syndiniales","Unknown Eukaryote")

# Set working directory
setwd("~/Desktop/SPOT/SPOT_2023")

s5 <- read.csv("Glasso_5m_SCC_MAY2023.csv",header=TRUE,row.names=1)
sDCM <- read.csv("Glasso_DCM_SCC_MAY2023.csv",header=TRUE,row.names=1)
colnames(s5) <- str_remove_all(colnames(s5),"^X")
colnames(sDCM) <- str_remove_all(colnames(sDCM),"^X")

# Let's make the edge weights positive (1) or negative (-1)
s5 <- as.matrix(s5)
sDCM <- as.matrix(sDCM)

s5[s5<0] <- -1
s5[s5>0] <- 1
sDCM[sDCM<0] <- -1
sDCM[sDCM>0] <- 1

# Make igraph object
g5 <- graph_from_adjacency_matrix(s5,mode="undirected",weighted=TRUE)
gdcm <- graph_from_adjacency_matrix(sDCM,mode="undirected",weighted=TRUE)

# Get intersection graph i.e., edges that are in both networks
int <- graph.intersection(g5,gdcm)

# Let's only include edges in the intersection network with the same SCC direction (i.e., positive or negative)
interEdge <- as.data.frame(get.edgelist(int))
g5Edge <- as.data.frame(get.edgelist(g5))
g5Edge$weight5 <- E(g5)$weight
gdcmEdge <- as.data.frame(get.edgelist(gdcm))
gdcmEdge$weightDCM <- E(gdcm)$weight

interEdge$fin1 <- paste(interEdge$V1,interEdge$V2,sep="_")
g5Edge$fin1 <- paste(g5Edge$V1,g5Edge$V2,sep="_")
g5Edge$fin2 <- paste(g5Edge$V2,g5Edge$V1,sep="_")
gdcmEdge$fin1 <- paste(gdcmEdge$V1,gdcmEdge$V2,sep="_")
gdcmEdge$fin2 <- paste(gdcmEdge$V2,gdcmEdge$V1,sep="_")

interEdge <- left_join(interEdge,g5Edge)
interEdge <- left_join(interEdge,gdcmEdge)

interEdge <- subset(interEdge,weight5==weightDCM)

tax <- read.delim("taxonomy_90.tsv",header=TRUE,row.names=NULL)
tax$Confidence <- NULL
colnames(tax) <- c("V1","Tax1")
interEdge <- left_join(interEdge,tax)
colnames(tax) <- c("V2","Tax2")
interEdge <- left_join(interEdge,tax)

# Make plot of intersection edges
interEdge <- subset(interEdge,select=c("V1","V2","weight5"))
gInterOnly <- graph_from_data_frame(interEdge,directed=FALSE)
gInterOnly$Weight <- interEdge$weight5

# pdf("Intersect_Network.pdf")
# plot(gInterOnly,vertex.size=4,vertex.label=NA,edge.width=0.5,layout=layout_with_kk(gInterOnly),vertex.color="grey",vertex.color="grey20",edge.color="olivedrab4",main="Associations found at surface and DCM (merged network)")
# dev.off()

# Network statistics
length(V(gInterOnly))
length(E(gInterOnly))
V(gInterOnly)$deg <- degree(gInterOnly)
mean(V(gInterOnly)$deg)
mean_distance(gInterOnly)
transitivity(gInterOnly,type="global")

# Find the degree of each node - see taxonomic affiliation
dfDeg <- data.frame(name=c(V(gInterOnly)$name),deg=c(V(gInterOnly)$deg))
dfDeg <- dfDeg %>% arrange(desc(deg)) %>% as.data.frame()
dfDeg$name <- str_replace_all(dfDeg$name,"S_","")
t <- read.delim("taxonomy_90.tsv",header=TRUE)
t$Confidence <- NULL
colnames(t)[1] <- "name"
dfDeg <- left_join(dfDeg,t)

# Add taxonomic information to visualize degree by taxonomic group
taxz <- colsplit(dfDeg$Taxon,";",c("Supergroup","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
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
dfDeg$tax <- taxz$fin

dfDeg %>% filter(!grepl("Other",dfDeg$tax))%>%ggplot(aes(x=tax,y=deg))+geom_boxplot()+theme_classic() # Maybe not that interesting

# Save intersection edgelist
# write.csv(interEdge,"Surf_DCM_edges_MAY2023.csv")

# Plot node size by degree
V(gInterOnly)$tax <- dfDeg$tax
taxColsDf <- data.frame(taxCols)
taxColsDf$tax <- rownames(taxColsDf)
dfDeg <- left_join(dfDeg,taxColsDf)
V(gInterOnly)$col <- dfDeg$taxCols

# pdf("../../Degree2.pdf")
plot(gInterOnly,vertex.size=V(gInterOnly)$deg*1.1,vertex.label=NA,edge.width=0.5,layout=layout_with_kk(gInterOnly),vertex.color=V(gInterOnly)$col,edge.color="grey30",main="")
# dev.off()

# Plot node size by average abundance
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
dfNew <- as.data.frame(clr(dfNew))

dfNew <- as.data.frame(t(dfNew))
dfNew <- subset(dfNew,rownames(dfNew) %in% V(gInterOnly)$name)
dfNew$sumz <- rowMeans(dfNew)
c <- data.frame(name=c(rownames(dfNew)),mean=c(dfNew$sumz))
dfDeg <- left_join(dfDeg,c)
m <- min(dfDeg$mean)
V(gInterOnly)$ab <- (dfDeg$mean)+abs(m)+0.01

# pdf("../../Abundance.pdf")
plot(gInterOnly,vertex.size=V(gInterOnly)$ab*3,vertex.label=NA,edge.width=0.5,layout=layout_with_kk(gInterOnly),vertex.color=V(gInterOnly)$col,edge.color="grey30",main="")
# dev.off()

ggplot(dfDeg,aes(x=deg,y=mean,fill=tax))+geom_point(shape=21,size=2)+scale_fill_manual(name=c("Taxonomic Group"),values=c(taxCols))+theme_classic(base_size = 14)+xlab("Number of Associations")+ylab("CLR-transformed Abundance")+theme(legend.position='top')+scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14))
# ggsave("../../plot1.pdf",width=8,height=4)
