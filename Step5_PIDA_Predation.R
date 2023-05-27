### SPOT Network Analysis ###
### PIDA Comparison - PREDATION ###
### By: Samantha Gleich ###
### Last Updated: 5/27/2023 ###

# Colors to use for all taxonomic analyses
taxCols <- c("#E1746D","#76C3D7","#DE5AB1","#D5E0AF","#DED3DC","#87EB58","#D4DC60","#88E5D3","#88AAE1","#DBA85C","#8B7DDA","#9A8D87","#D99CD1","#B649E3","#7EDD90")
names(taxCols) <- c("Chlorophyte","Ciliate","Cryptophyte","Diatom","Dinoflagellate","Fungi","Haptophyte","MAST","Other Alveolate","Other Archaeplastida","Other Eukaryote","Other Stramenopile","Rhizaria","Syndiniales","Unknown Eukaryote")

interEdge <- read.csv("Surf_DCM_edges_MAY2023.csv",header=TRUE)

tax <- read.delim("taxonomy_90.tsv",header=TRUE,row.names=NULL)
tax$Confidence <- NULL
colnames(tax) <- c("V1","Tax1")
interEdge <- left_join(interEdge,tax)
colnames(tax) <- c("V2","Tax2")
interEdge <- left_join(interEdge,tax)
# write.csv(interEdge,"../PIDA_Prelim.csv")

pida <- read.csv("PIDA_Int.csv",header=TRUE) 
pida <- subset(pida,Taxonomic.interaction=="Prot - Prot")
pidaPar <- subset(pida,Interaction=="par")
pidaSymb <- subset(pida,Interaction=="symb")
pidaPred <- subset(pida,Interaction=="pred")

# Let's start by making our parasitism subnetwork
inter1 <- colsplit(interEdge$Tax1,";",c("super1","k1","p1","c1","o1","f1","g1","s1"))
inter2 <- colsplit(interEdge$Tax2,";",c("super2","k2","p2","c2","o2","f2","g2","s2"))
interEdge <- cbind(interEdge,inter1,inter2)

subs1 <- subset(pidaPred,Genus.org1 %in% interEdge$g1 & Genus.org2 %in% interEdge$g2)
subs2 <- subset(pidaPred,Genus.org2 %in% interEdge$g1 & Genus.org1 %in% interEdge$g2)
subsFin <- data.frame(Var1=c(subs1$Genus.org1,subs2$Genus.org1),Var2=c(subs1$Genus.org2,subs2$Genus.org2))
subsFin$Var3 <- paste(subsFin$Var1,subsFin$Var2,sep="_")
subsFin$Var4 <- paste(subsFin$Var2,subsFin$Var1,sep="_")

interEdge$Var <- paste(interEdge$g1,interEdge$g2,sep="_")
outPred1 <- subset(interEdge,Var %in% subsFin$Var3)
outPred2 <- subset(interEdge,Var %in% subsFin$Var4)
outPred <- rbind(outPred1,outPred2)
outPred$X <- NULL

# Plot 
outG <-graph_from_data_frame(outPred,directed=FALSE)
E(outG)$weight <- outPred$weight
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
taxz$Low <- taxz$Species
taxz$Low <- ifelse(taxz$Species=="",taxz$Genus,taxz$Low)
taxz$Low <- str_replace_all(taxz$Low,"_sp.;","")
taxz$Low <- str_replace_all(taxz$Low,"_"," ")
taxz$Low <- str_replace_all(taxz$Low,";","")
V(finalG)$namez <- taxz$Low



# PLOT
pdf("../../PREDATION_APRIL5.pdf",width=20,height=15)
ggraph(finalG, layout = 'linear', circular = TRUE) + geom_edge_arc(aes(color = weight),alpha=0.6,width=1) + 
  geom_node_text(aes(label = namez, angle = node_angle(x, y)), hjust = -0.15,size=6,fontface="italic") +
  geom_node_point(shape = 21, size = 6, aes(fill = fin)) +
  theme_graph() +
  scale_edge_color_gradientn(colours = c('indianred', 'white', 'dodgerblue'),values=c(-1,0,1)) +
  coord_fixed(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4)) +scale_fill_manual(name="Taxonomic Groups",values=c(taxCols))
dev.off()

# LEGEND
outG2$x <- 1:11
outG2$y <- 2:12
ggplot(outG2,aes(x=x,y=y,fill=fin1))+geom_point(shape=21,size=3)+scale_fill_manual(name="Taxonomic Groups",values=c(taxCols))+theme_classic()
ggsave("../../PredationLegend.pdf")
