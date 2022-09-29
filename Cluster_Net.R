### Load eLSA Output and Filter ###
# Load data
surf <- read.delim(file.choose(),header=TRUE)
dcm <- read.delim(file.choose(),header=TRUE)

# Subset based on p and q-value
surfSig <- subset(surf,Psscc<0.01 & Qsscc < 0.01)
dcmSig <- subset(dcm,Psscc<0.01 & Qsscc < 0.01)

# Subset lagged and unlagged surface
surfLag <- subset(surfSig,Delay!=0)
surfNoLag <- subset(surfSig,Delay==0)
surfLagG <- data.frame(x=c(surfLag$X),y=c(surfLag$Y),scc=c(surfLag$SCC),delay=c(surfLag$Delay))
surfNoLagG <- data.frame(x=c(surfNoLag$X),y=c(surfNoLag$Y),scc=c(surfNoLag$SCC),delay=c(surfNoLag$Delay))

### Plot lagged and unlagged networks surface ###
# Unlagged surface
g <- graph_from_data_frame(surfNoLagG,directed=FALSE)

plot(g,vertex.size=4,vertex.label=NA,edge.width=0.1,layout=layout_with_lgl(g))

# Lagged surface
surfLagG$new <- ifelse(surfLagG$delay==-1,surfLagG$y,NA)
surfLagG$y <- ifelse(surfLagG$delay==-1,surfLagG$x,surfLagG$y)
surfLagG$x <- ifelse(surfLagG$delay==-1,surfLagG$new,surfLagG$x)
surfLagG$new <- NULL

gLag <- graph_from_data_frame(surfLagG, directed=TRUE)
plot(gLag,vertex.size=4,vertex.label=NA,edge.width=0.1,layout=layout_with_lgl(gLag),edge.arrow.size=0.2)

# Subset lagged and unlagged DCM
DCMLag <- subset(dcmSig,Delay!=0)
DCMNoLag <- subset(dcmSig,Delay==0)
DCMLagG <- data.frame(x=c(DCMLag$X),y=c(DCMLag$Y),scc=c(DCMLag$SCC),delay=c(DCMLag$Delay))
DCMNoLagG <- data.frame(x=c(DCMNoLag$X),y=c(DCMNoLag$Y),scc=c(DCMNoLag$SCC),delay=c(DCMNoLag$Delay))

### Plot lagged and unlagged networks DCM ###
# Unlagged DCM
g2 <- graph_from_data_frame(DCMNoLagG,directed=FALSE)

plot(g2,vertex.size=4,vertex.label=NA,edge.width=0.1,layout=layout_with_lgl(g2))

# Lagged DCM
DCMLagG$new <- ifelse(DCMLagG$delay==-1,DCMLagG$y,NA)
DCMLagG$y <- ifelse(DCMLagG$delay==-1,DCMLagG$x,DCMLagG$y)
DCMLagG$x <- ifelse(DCMLagG$delay==-1,DCMLagG$new,DCMLagG$x)
DCMLagG$new <- NULL

gLag2 <- graph_from_data_frame(DCMLagG, directed=TRUE)
plot(gLag2,vertex.size=4,vertex.label=NA,edge.width=0.1,layout=layout_with_lgl(gLag2),edge.arrow.size=0.2)

### Find edges that are in surface and DCM networks ###
tax <- read.csv(file.choose(),header=TRUE,row.names=1)
int <- intersection(g,g2) # Only edges that are in both
un <- graph.union(g,g2) # All edges in both

E(un)$col <- ifelse(E(un) %in% E(int),"darkgoldenrod4","indianred")
E(un)$col <- ifelse(E(un) %in% E(g2) & E(un)$col=="indianred","navy",E(un)$col)
unEdge <- data.frame(get.edgelist(un))
unEdge$col <- E(un)$col

tax2 <- data.frame(X1=tax$num,Final=tax$Final)
unEdge <- left_join(unEdge,tax2)
colnames(tax2) <- c("X2","Final2")
unEdge <- left_join(unEdge,tax2)

taxSum <- unEdge %>% group_by(Final,Final2,col) %>% tally()%>%as.data.frame()
taxSum$Assoc <- paste(taxSum$Final,taxSum$Final2,sep="-")

ggplot(taxSum,aes(x=reorder(Assoc, -n),y=n,fill=as.factor(col)))+geom_bar(stat="identity")+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12))+scale_fill_manual(name="Depth",values=c("darkgoldenrod3","indianred","navy"),labels=c("Surface & DCM","Surface Only","DCM Only"))+ylab("Number of Associations")+xlab("Association Type")
# ggsave("Figure5_Prelim.pdf",width=18,height=8)


### Cluster network ###
# Delete verticies without edges
int <- delete.vertices(igraph::simplify(int), degree(int)==0)
length(E(int))
length(V(int))

clus <- cluster_louvain(int)
length(clus)

V(int)$community <- clus$membership

pdf("FigureX_Prelim.pdf")  
plot(int,vertex.label=NA,layout=layout_with_fr(int),vertex.size=3,mark.border="black",edge.width=0.3,mark.groups=communities(clus),vertex.color="grey",edge.color="black",mark.col=adjustcolor(c("coral1","deepskyblue2","darkseagreen1","mediumpurple"),alpha=0.4))
dev.off()
