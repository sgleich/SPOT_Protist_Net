### SPOT Network Analysis ###
### PIDA Comparison ###
### By: Samantha Gleich ###
### Last Updated: 1/2/23 ###

# Colors to use for all taxonomic analyses
taxCols <- c("#E1746D","#76C3D7","#DE5AB1","#D5E0AF","#DED3DC","#87EB58","#D4DC60","#88E5D3","#88AAE1","#DBA85C","#8B7DDA","#9A8D87","#D99CD1","#B649E3","#7EDD90")
names(taxCols) <- c("Chlorophyte","Ciliate","Cryptophyte","Diatom","Dinoflagellate","Fungi","Haptophyte","MAST","Other Alveolate","Other Archaeplastida","Other Eukaryote","Other Stramenopile","Rhizaria","Syndiniales","Unknown Eukaryote")

interEdge <- read.csv("Surf_DCM_Edges.csv",header=TRUE,row.names=1)
interEdge$V1 <- str_remove_all(interEdge$V1,"S_")
interEdge$V2 <- str_remove_all(interEdge$V2,"S_")
tax <- read.delim("taxonomy_90.tsv",header=TRUE,row.names=NULL)
tax$Confidence <- NULL
colnames(tax) <- c("V1","Tax1")
interEdge <- left_join(interEdge,tax)
colnames(tax) <- c("V2","Tax2")
interEdge <- left_join(interEdge,tax)
write.csv(interEdge,"../PIDA_Prelim.csv")

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
keeperz <- keeperz[-(4)]
keeperz <- keeperz[-(7)]

dinoSyn <- keeperz[c(3:6,8:10)]
cilSyn <- keeperz[c(1:2)]
diaCer <- keeperz[c(7,11)]

interEdge <- interEdge[c(1:4)]
synG1 <- subset(interEdge,grepl("Group-I-",interEdge$Tax1)|grepl("Group-I-",interEdge$Tax2))
synG2 <- subset(interEdge,grepl("Group-II-",interEdge$Tax1)|grepl("Group-II-",interEdge$Tax2))
cerc<- subset(interEdge,grepl("Cercoz",interEdge$Tax1)|grepl("Cercoz",interEdge$Tax2))

# Group I Syndiniales (Infect ciliates)
Strob <- subset(synG1,grepl("Strombidium",synG1$Tax1)|grepl("Strombidium",synG1$Tax2))
Pelago <- subset(synG1,grepl("Pelagostrobilidium",synG1$Tax1)|grepl("Pelagostrobilidium",synG1$Tax2))

# Group II Syndiniales (Infect Dinoflagellate)
Proro <- subset(synG2,grepl("Prorocentrum",synG2$Tax1)|grepl("Prorocentrum",synG2$Tax2))
Gyro <- subset(synG2,grepl("Gyrodinium",synG2$Tax1)|grepl("Gyrodinium",synG2$Tax2))
Gony <- subset(synG2,grepl("Gonyaulax",synG2$Tax1)|grepl("Gonyaulax",synG2$Tax2))
Trip <- subset(synG2,grepl("Tripos",synG2$Tax1)|grepl("Tripos",synG2$Tax2))
Hetero<- subset(synG2,grepl("Heterocapsa",synG2$Tax1)|grepl("Heterocapsa",synG2$Tax2))
Gymno<- subset(synG2,grepl("Gymnodinium",synG2$Tax1)|grepl("Gymnodinium",synG2$Tax2))
Karlo<- subset(synG2,grepl("Karlodinium",synG2$Tax1)|grepl("Karlodinium",synG2$Tax2))

# Cercozoa no matches to either diatom genus

full <- rbind(Strob,Pelago,Proro,Gyro,Gony,Trip,Hetero,Gymno,Karlo)
full <- full[c(1:2)]

outPara <- as.matrix(full)

# Plot 
outG <-graph_from_data_frame(outPara,directed=FALSE)
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

outDf$Final2 <- taxz$fin
taxz$Low <- paste(taxz$Genus,taxz$Species,sep="_")
taxz$Low <- ifelse(taxz$Class=="Syndiniales",taxz$Order,taxz$Low)
outDf$Low <- taxz$Low

outDfNames <- data.frame(outDfNames=c(unique(outDf$Low)),num=c(1))
for (i in 1:nrow(outDf)){
  m <- outDf[i,4]
  num <- subset(outDfNames,outDfNames==m)
  outDf$Low[i] <- paste(outDf$Low[i],num$num,sep="_")
  row <- which(outDfNames==m)
  outDfNames[row,2] <- outDfNames[row,2]+1
}

outSum <- outDf %>% group_by(Final2) %>% tally()

chlor <- subset(outDf,Final2=="Chlorophyte")
cil <- subset(outDf,Final2=="Ciliate")
cry <- subset(outDf,Final2=="Cryptophyte")
dia <- subset(outDf,Final2=="Diatom")
dino <- subset(outDf,Final2=="Dinoflagellate")
hapto <- subset(outDf,Final2=="Haptophyte")
mast <- subset(outDf,Final2=="MAST")
oeuk <- subset(outDf,Final2=="Other Eukaryote")
ostr <- subset(outDf,Final2=="Other Stramenopile")
rhiz <- subset(outDf,Final2=="Rhizaria")
syn <- subset(outDf,Final2=="Syndiniales")
unk <- subset(outDf,Final2=="Unknown Eukaryote")

v <- c(outSum$Final2,chlor$num,cil$num,cry$num,dia$num,dino$num,hapto$num,mast$num,oeuk$num,ostr$num,rhiz$num,syn$num,unk$num)

# Hierarchy df
hierarchy <- data.frame(from=outDf$Final2,to=outDf$Low)
tmp <- data.frame(from=c(rep("Origin",3)),to=c("Dinoflagellate","Syndiniales","Ciliate"))
hierarchy <- rbind(tmp,hierarchy)

# Vertices df
vertices <- data.frame(name=outDf$Low,value=1,group=outDf$Final2)
tmp2 <- data.frame(name=c("Origin","Dinoflagellate","Syndiniales","Ciliate"),value=c(1,1,1,1),group=c(NA,"Origin","Origin","Origin"))
vertices <- rbind(tmp2,vertices)

# Mygraph
out <- as.data.frame(outPara)
join <- data.frame(V1=c(outDf$V2),Low1=c(outDf$Low))
out <- left_join(out,join)
colnames(join) <- c("V2","Low2")
out <- left_join(out,join)

mygraph <- graph_from_data_frame( hierarchy, vertices=vertices)

connect <- data.frame(from=c(out$Low1),to=c(out$Low2),value=c(1))
from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)

#Let's add information concerning the label we are going to add: angle, horizontal adjustement and potential flip
#calculate the ANGLE of the labels
vertices$id <- NA
myleaves <- which(is.na( match(vertices$name, hierarchy$from) ))
nleaves <- length(myleaves)
vertices$id[ myleaves ] <- seq(1:nleaves)
vertices$angle <- 90 - 360 * vertices$id / nleaves

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)

# flip angle BY to make them readable
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)

# flip angle BY to make them readable
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

mygraph <- igraph::graph_from_data_frame( hierarchy, vertices=vertices )

# The connection object must refer to the ids of the leaves:
from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)


# Plot network
# v <- randomcoloR::distinctColorPalette(15)

ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05)) +
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.5, colour="grey", width=0.7) +
  #geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, label=name, angle = angle, hjust=hjust), size=1.5, alpha=1) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))+theme(legend.position = "right")+geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group),   size=2) +
  scale_colour_manual(name="Taxonomic Groups",values= c(taxCols))+ggtitle("")
ggsave("../Parasitism.pdf",width=8,height=6)
