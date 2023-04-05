### SPOT Network Analysis ###
### PIDA Comparison - PARASITISM ###
### By: Samantha Gleich ###
### Last Updated: 4/5/23 ###

# Colors to use for all taxonomic analyses
taxCols <- c("#E1746D","#76C3D7","#DE5AB1","#D5E0AF","#DED3DC","#87EB58","#D4DC60","#88E5D3","#88AAE1","#DBA85C","#8B7DDA","#9A8D87","#D99CD1","#B649E3","#7EDD90")
names(taxCols) <- c("Chlorophyte","Ciliate","Cryptophyte","Diatom","Dinoflagellate","Fungi","Haptophyte","MAST","Other Alveolate","Other Archaeplastida","Other Eukaryote","Other Stramenopile","Rhizaria","Syndiniales","Unknown Eukaryote")

interEdge <- read.csv(file.choose(),header=TRUE)

tax <- read.delim("taxonomy_90.tsv",header=TRUE,row.names=NULL)
tax$Confidence <- NULL
colnames(tax) <- c("V1","Tax1")
interEdge$V1 <- str_remove_all(interEdge$V1,"S_")
interEdge$V2 <- str_remove_all(interEdge$V2,"S_")
interEdge <- left_join(interEdge,tax)
colnames(tax) <- c("V2","Tax2")
interEdge <- left_join(interEdge,tax)
# write.csv(interEdge,"../PIDA_Prelim.csv")

pida <- read.csv(file.choose(),header=TRUE) 
pida <- subset(pida,Taxonomic.interaction=="Prot - Prot")
pidaPar <- subset(pida,Interaction=="par")
pidaSymb <- subset(pida,Interaction=="symb")
pidaPred <- subset(pida,Interaction=="pred")

# Let's start by making our parasitism subnetwork
inter1 <- colsplit(interEdge$Tax1,";",c("super1","k1","p1","c1","o1","f1","g1","s1"))
inter2 <- colsplit(interEdge$Tax2,";",c("super2","k2","p2","c2","o2","f2","g2","s2"))
interEdge <- cbind(interEdge,inter1,inter2)

genera <- c(unique(interEdge$g1,unique(interEdge$g2)))
genera <- unique(genera)
genera <- colsplit(genera,"_",c("keep","other"))
genera <- genera$keep

pidaGenera <- c(unique(pidaPar$Genus.org1,unique(pidaPar$Genus.org2)))
pidaGenera <- unique(pidaGenera)
keeperz <- NULL
for (i in 1:length(genera)){
  g <- genera[i]
  if (g %in% pidaGenera){
    print(g)
    keeperz <- c(keeperz,g)
  }
}


synG1 <- subset(interEdge,grepl("Group-I-",interEdge$Tax1)|grepl("Group-I-",interEdge$Tax2))
synG2 <- subset(interEdge,grepl("Group-II-",interEdge$Tax1)|grepl("Group-II-",interEdge$Tax2))
cerc<- subset(interEdge,grepl("Cercoz",interEdge$Tax1)|grepl("Cercoz",interEdge$Tax2))

# Group I Syndiniales (Infect ciliates)
Strob <- subset(synG1,grepl("Strombidium",synG1$Tax1)|grepl("Strombidium",synG1$Tax2))
Pelago <- subset(synG1,grepl("Pelagostrobilidium",synG1$Tax1)|grepl("Pelagostrobilidium",synG1$Tax2))
Eut <- subset(synG1,grepl("Eutintinnus",synG1$Tax1)|grepl("Eutintinnus",synG1$Tax2))

# Group II Syndiniales (Infect Dinoflagellate)
Proro <- subset(synG2,grepl("Prorocentrum",synG2$Tax1)|grepl("Prorocentrum",synG2$Tax2))
Gyro <- subset(synG2,grepl("Gyrodinium",synG2$Tax1)|grepl("Gyrodinium",synG2$Tax2))
Gony <- subset(synG2,grepl("Gonyaulax",synG2$Tax1)|grepl("Gonyaulax",synG2$Tax2))
Trip <- subset(synG2,grepl("Tripos",synG2$Tax1)|grepl("Tripos",synG2$Tax2))
Hetero<- subset(synG2,grepl("Heterocapsa",synG2$Tax1)|grepl("Heterocapsa",synG2$Tax2))
Gymno<- subset(synG2,grepl("Gymnodinium",synG2$Tax1)|grepl("Gymnodinium",synG2$Tax2))
Karlo<- subset(synG2,grepl("Karlodinium",synG2$Tax1)|grepl("Karlodinium",synG2$Tax2))
Lepid<- subset(synG2,grepl("Lepidodinium",synG2$Tax1)|grepl("Lepidodinium",synG2$Tax2))
Alex<- subset(synG2,grepl("Alexandrium",synG2$Tax1)|grepl("Alexandrium",synG2$Tax2))
Proto<- subset(synG2,grepl("Protoperidinium",synG2$Tax1)|grepl("Protoperidinium",synG2$Tax2))

# Cercozoa no matches to either diatom genus
Cerc<- subset(interEdge,grepl("Cryothe",interEdge$Tax1)|grepl("Cryotheco",interEdge$Tax2))
Cerc<- subset(Cerc,grepl("Bacillari",Cerc$Tax1)|grepl("Bacillari",Cerc$Tax2))

full <- rbind(Strob,Pelago,Eut,Proro,Gyro,Gony,Trip,Hetero,Gymno,Karlo,Lepid,Alex,Proto,Cerc)
full <- full[c(2:4)]


# Plot 
outG <-graph_from_data_frame(full,directed=FALSE)
E(outG)$weight <- full$weight
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

V(outG)$fin <- taxz$fin
taxz$Low <- paste(taxz$Genus, taxz$Species,sep="_")
taxz$Low <- ifelse(taxz$Low=="_",taxz$Family,taxz$Low)
low <- colsplit(taxz$Low,"_",c("keep","keep2","toss"))
V(outG)$namez <- low$keep

outG2 <- as.data.frame(get.edgelist(outG))
outG2$edge <- E(outG)$weight
outG3 <- data.frame(V1=c(V(outG)$name)) 
outG3$fin1 <- V(outG)$fin
# outG3$namez <- V(outG)$namez
outG2 <- left_join(outG2,outG3)
outG2 <- outG2 %>% arrange(desc(fin1)) %>% as.data.frame()
colnames(outG3) <- c("V2","fin2")
outG2 <- left_join(outG2,outG3)
outG2 <- outG2 %>% arrange(desc(fin2)) %>% as.data.frame()

outG2$swap <- ifelse(outG2$fin1=="Syndiniales","yes","no")
outG2$V1b <- outG2$V1
outG2$V2b <- outG2$V2
outG2$V1 <- ifelse(outG2$swap=="yes",outG2$V2b,outG2$V1)
outG2$V2 <- ifelse(outG2$swap=="yes",outG2$V1b,outG2$V2)
outG2 <- outG2[c(1:6)]
outG2$fin1b <- outG2$fin1
outG2$fin2b <- outG2$fin2
outG2$fin1 <- ifelse(outG2$swap=="yes",outG2$fin2b,outG2$fin1)
outG2$fin2 <- ifelse(outG2$swap=="yes",outG2$fin1b,outG2$fin2)
outG2 <- outG2[c(1:6)]

finalG <- graph_from_data_frame(outG2,directed=FALSE)
E(finalG)$weight <- outG2$edge

outDf <- data.frame(name=c(V(finalG)$name))
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

V(finalG)$fin <- taxz$fin
taxz$Low <- paste(taxz$Genus, taxz$Species,sep="_")
taxz$Low <- ifelse(taxz$Low=="_",taxz$Family,taxz$Low)
low <- colsplit(taxz$Low,"_",c("keep","keep2","toss"))
V(finalG)$namez <- low$keep

# PLOT
pdf("../../PARASITISM_APRIL5.pdf",width=20,height=15)
ggraph(finalG, layout = 'linear', circular = TRUE) + geom_edge_arc(aes(color = weight),alpha=0.6,width=1) + 
  geom_node_text(aes(label = namez, angle = node_angle(x, y)), hjust = -0.15,size=3) +
  geom_node_point(shape = 21, size = 6, aes(fill = fin)) +
  theme_graph() +
  scale_edge_color_gradientn(colours = c('indianred', 'white', 'dodgerblue'),values=c(-1,0,1)) +
  coord_fixed(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4)) +scale_fill_manual(name="Taxonomic Groups",values=c(taxCols))
dev.off()

# LEGEND
outG3$x <- 1:21
outG3$y <- 2:22
ggplot(outG3,aes(x=x,y=y,fill=fin2))+geom_point(shape=21,size=3)+scale_fill_manual(name="Taxonomic Groups",values=c(taxCols))+theme_classic()
ggsave("../../ParasitismLegend.pdf")
