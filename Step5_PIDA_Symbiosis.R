### SPOT Network Analysis ###
### PIDA Comparison - MUTUALISM###
### By: Samantha Gleich ###
### Last Updated: 4/5/23 ###

# Colors to use for all taxonomic analyses
taxCols <- c("#E1746D","#76C3D7","#DE5AB1","#D5E0AF","#DED3DC","#87EB58","#D4DC60","#88E5D3","#88AAE1","#DBA85C","#8B7DDA","#9A8D87","#D99CD1","#B649E3","#7EDD90")
names(taxCols) <- c("Chlorophyte","Ciliate","Cryptophyte","Diatom","Dinoflagellate","Fungi","Haptophyte","MAST","Other Alveolate","Other Archaeplastida","Other Eukaryote","Other Stramenopile","Rhizaria","Syndiniales","Unknown Eukaryote")

interEdge <- read.csv(file.choose(),header=TRUE)

tax <- read.delim("taxonomy_90.tsv",header=TRUE,row.names=NULL)
tax$Confidence <- NULL
colnames(tax) <- c("V1","Tax1")
interEdge <- left_join(interEdge,tax)
colnames(tax) <- c("V2","Tax2")
interEdge <- left_join(interEdge,tax)

pida <- read.csv(file.choose(),header=TRUE) 
pida <- subset(pida,Taxonomic.interaction=="Prot - Prot")
pidaPar <- subset(pida,Interaction=="par")
pidaSymb <- subset(pida,Interaction=="symb")
pidaPred <- subset(pida,Interaction=="pred")

# Let's start by making our parasitism subnetwork
inter1 <- colsplit(interEdge$Tax1,";",c("super1","k1","p1","c1","o1","f1","g1","s1"))
inter2 <- colsplit(interEdge$Tax2,";",c("super2","k2","p2","c2","o2","f2","g2","s2"))
interEdge <- cbind(interEdge,inter1,inter2)

subs <- subset(pidaSymb,Genus.org1 %in% interEdge$g1 | Genus.org2 %in% interEdge$g1 |Genus.org1 %in% interEdge$g2 | Genus.org2 %in% interEdge$g2)

acanth <- c("Acanthometron","Lithoptera","Arthracanthida","Phyllostaurus","Lonchostaurus","Amphilonche","Xiphacantha","Lonchostaurus","Stauracantha","Amphibelone","Acanthostaurus") # Clade F
other <- c("Askenasia","Leptocylindrus","Chaetoceros","Dictyocha")
#######
subsAcanth <- subset(subs,grepl("Acanthar",Taxonomic.level.3..org1)| grepl("Acanthar",Taxonomic.level.3..org2))
subsAcanth <- unique(subsAcanth$Genus.org2)
acanthEdges <- subset(interEdge,grepl("Acantharea",interEdge$Tax2) | grepl("Acantharea",interEdge$Tax1))
acanthEdges1 <- subset(acanthEdges,grepl("Chrysochrom",acanthEdges$Tax2) | grepl("Chrysochrom",acanthEdges$Tax1))
acanthEdges2 <- subset(acanthEdges,grepl("Phaeo",acanthEdges$Tax2) | grepl("Phaeo",acanthEdges$Tax1))
subsOther <- subset(subs,Genus.org1 %in% other| Genus.org2 %in% other)
otherEdges <- subset(interEdge,grepl("Lepto",interEdge$Tax1) & grepl("MAST-4",interEdge$Tax2))
outSymb<- rbind(otherEdges,acanthEdges1,acanthEdges2)
######
acan <- subset(interEdge,grepl("Acantharea_F",interEdge$Tax1)|grepl("Acantharea_F",interEdge$Tax2))
poly <- subset(interEdge,grepl("Polycystine",interEdge$Tax1)|grepl("Polycystine",interEdge$Tax2))
acanDino <- subset(acan, grepl("Dinophyceae",acan$Tax1)|grepl("Dinophyceae",acan$Tax2))
polyDino <- subset(poly,grepl("Dinophyceae",poly$Tax1)|grepl("Dinophyceae",poly$Tax2))
polyPrym <- subset(poly,grepl("Prymnesiophy",poly$Tax1)|grepl("Prymnesiophy",poly$Tax2))

outSymb <- rbind(acanDino,polyDino,polyPrym)
outSymb$X <- NULL
outSymb <- outSymb[c(1:3)]

outSymb$Tax1 <- NULL
outSymb$Tax2 <- NULL

# Plot 
outG <-graph_from_data_frame(outSymb,directed=FALSE)
E(outG)$weight <- outSymb$weight
outDf <- data.frame(name=c(V(outG)$name))
colnames(outDf) <- "V2"
outDf <- left_join(outDf,tax)

taxz <- colsplit(outDf$Tax2,";",c("Supergroup","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxz$fin <- ifelse(taxz$Class=="Syndiniales","Syndiniales",NA)
taxz$fin <- ifelse(taxz$Class=="Dinophyceae","Dinoflagellate",taxz$fin)
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

taxz <- colsplit(outDf$Tax2,";",c("Supergroup","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxz$fin <- ifelse(taxz$Class=="Syndiniales","Syndiniales",NA)
taxz$fin <- ifelse(taxz$Class=="Dinophyceae","Dinoflagellate",taxz$fin)
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

V(outG)$fin <- taxz$fin
taxz$Low <- paste(taxz$Genus, taxz$Species,sep="_")
taxz$Low <- ifelse(taxz$Low=="_",taxz$Family,taxz$Low)
taxz$Low <- ifelse(taxz$Low=="",taxz$Order,taxz$Low)
taxz$Low <- ifelse(taxz$Low=="",taxz$Class,taxz$Low)
low <- colsplit(taxz$Low,"_",c("keep","keep2","toss"))
low$keep <- ifelse(low$keep=="Acantharea","Acantharea Clade F",low$keep)
low$keep <- ifelse(low$keep=="Cannobotryidae","Cannobotryidae (Polycystine)",low$keep)
low$keep <- ifelse(low$keep=="Lophophaenidae","Lophophaenidae (Polycystine)",low$keep)
V(outG)$namez <- low$keep

pdf("../../try2.pdf",width=20,height=11)
ggraph(outG, layout = 'linear', circular = TRUE) + geom_edge_arc(color="dodgerblue",alpha=0.6,width=1) + 
  geom_node_text(aes(label = namez, angle = node_angle(x, y)), hjust = -0.1,size=3) +
  geom_node_point(shape = 21, size = 4, aes(fill = fin)) +
  theme_graph() +scale_edge_color_gradientn(colours = c('indianred', 'white', 'dodgerblue'),values=c(-1,0,1))+
  coord_fixed(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4)) +scale_fill_manual(name="Taxonomic Groups",values=c(taxCols))
dev.off()


# LEGEND
low$fin <- c("Haptophyte","Haptophyte","Rhizaria","Rhizaria")
low$x <- 1:4
low$y <- 2:5
ggplot(low,aes(x=x,y=y,fill=fin))+geom_point(shape=21,size=3)+scale_fill_manual(name="Taxonomic Groups",values=c(taxCols))+theme_classic()
ggsave("../../MutualismLegend.pdf")
