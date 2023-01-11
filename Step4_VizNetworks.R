### SPOT Network Analysis ###
### Visualize Networks + Calculate Network Statistics ###
### By: Samantha Gleich ###
### Last Updated: 12/28/22 ###

# Load in network outputs
surf <- read.csv(file.choose(),header=TRUE,row.names=1)
dcm <- read.csv(file.choose(),header=TRUE,row.names=1)

surfName <- colnames(surf)
surf <- as.matrix(surf)
gSurf <- graph_from_adjacency_matrix(surf,mode="undirected")
surfEdge <- get.edgelist(gSurf)
surfEdge <- as.data.frame(surfEdge)

dcmName <- colnames(dcm)
dcm <- as.matrix(dcm)
gDCM <- graph_from_adjacency_matrix(dcm,mode="undirected")
dcmEdge <- get.edgelist(gDCM)
dcmEdge <- as.data.frame(dcmEdge)

# Find edges in surface and DCM network 
inter <- graph.intersection(gSurf,gDCM)
interEdge <- get.edgelist(inter)
interEdge <- as.data.frame(interEdge)

# Find edges that are only in surface or only in DCM network 
`%ni%` <- Negate(`%in%`)
surfEdge$fin <- paste(surfEdge$V1,surfEdge$V2,sep="_")
dcmEdge$fin <- paste(dcmEdge$V1,dcmEdge$V2,sep="_")
interEdge$fin1 <- paste(interEdge$V1,interEdge$V2,sep="_")
interEdge$fin2 <- paste(interEdge$V2,interEdge$V1,sep="_")

surfEdgeOnly <- subset(surfEdge,fin %ni% interEdge$fin1 & fin %ni% interEdge$fin2)
dcmEdgeOnly <- subset(dcmEdge,fin %ni% interEdge$fin1 & fin %ni% interEdge$fin2)

surfEdgeOnly$fin <- NULL
dcmEdgeOnly$fin <- NULL
interEdge$fin1 <- NULL
interEdge$fin2 <- NULL

surfEdgeOnly <- as.matrix(surfEdgeOnly)
dcmEdgeOnly <- as.matrix(dcmEdgeOnly)
interEdge<- as.matrix(interEdge)

gSurfOnly <- graph_from_edgelist(surfEdgeOnly,directed=FALSE)
gDcmOnly <- graph_from_edgelist(dcmEdgeOnly,directed=FALSE)
gInterOnly <- graph_from_edgelist(interEdge,directed=FALSE)

# Visualize networks
pdf("Surface_Network.pdf")
plot(gSurfOnly,vertex.size=4,vertex.label=NA,edge.width=0.5,layout=layout_with_kk(gSurfOnly),vertex.color="grey",vertex.color="grey20",edge.color="indianred",main="Associations only found at surface (5 m)")
dev.off()

pdf("DCM_Network.pdf")
plot(gDcmOnly,vertex.size=4,vertex.label=NA,edge.width=0.5,layout=layout_with_kk(gDcmOnly),vertex.color="grey",vertex.color="grey20",edge.color="deepskyblue4",main="Associations only found at DCM")
dev.off()

pdf("Intersect_Network.pdf")
plot(gInterOnly,vertex.size=4,vertex.label=NA,edge.width=0.5,layout=layout_with_kk(gInterOnly),vertex.color="grey",vertex.color="grey20",edge.color="olivedrab4",main="Associations found at surface and DCM (merged network)")
dev.off()

# Network statistics
length(V(gInterOnly))
length(E(gInterOnly))
V(gInterOnly)$deg <- degree(gInterOnly)
mean(V(gInterOnly)$deg)
mean_distance(gInterOnly)
transitivity(gInterOnly,type="global")

dfDeg <- data.frame(name=c(V(gInterOnly)$name),deg=c(V(gInterOnly)$deg))
dfDeg <- dfDeg %>% arrange(desc(deg)) %>% as.data.frame()
dfDeg$name <- str_replace_all(dfDeg$name,"S_","")
t <- read.delim("taxonomy_90.tsv",header=TRUE)
t$Confidence <- NULL
colnames(t)[1] <- "name"
dfDeg <- left_join(dfDeg,t)

# Save edgelists 
interEdge <- as.data.frame(interEdge)
surfEdgeOnly <- as.data.frame(surfEdgeOnly)
dcmEdgeOnly <- as.data.frame(dcmEdgeOnly)
write.csv(interEdge,"Surf_DCM_edges.csv")
write.csv(surfEdgeOnly,"Surf_edges.csv")
write.csv(dcmEdgeOnly,"DCM_edges.csv")
