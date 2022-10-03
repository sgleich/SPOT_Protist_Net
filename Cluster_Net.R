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

# plot(g,vertex.size=4,vertex.label=NA,edge.width=0.1,layout=layout_with_lgl(g))

# Lagged surface
surfLagG$new <- ifelse(surfLagG$delay==-1,surfLagG$y,NA)
surfLagG$y <- ifelse(surfLagG$delay==-1,surfLagG$x,surfLagG$y)
surfLagG$x <- ifelse(surfLagG$delay==-1,surfLagG$new,surfLagG$x)
surfLagG$new <- NULL

gLag <- graph_from_data_frame(surfLagG, directed=TRUE)
# plot(gLag,vertex.size=4,vertex.label=NA,edge.width=0.1,layout=layout_with_lgl(gLag),edge.arrow.size=0.2)

# Subset lagged and unlagged DCM
DCMLag <- subset(dcmSig,Delay!=0)
DCMNoLag <- subset(dcmSig,Delay==0)
DCMLagG <- data.frame(x=c(DCMLag$X),y=c(DCMLag$Y),scc=c(DCMLag$SCC),delay=c(DCMLag$Delay))
DCMNoLagG <- data.frame(x=c(DCMNoLag$X),y=c(DCMNoLag$Y),scc=c(DCMNoLag$SCC),delay=c(DCMNoLag$Delay))

### Plot lagged and unlagged networks DCM ###
# Unlagged DCM
g2 <- graph_from_data_frame(DCMNoLagG,directed=FALSE)

# plot(g2,vertex.size=4,vertex.label=NA,edge.width=0.1,layout=layout_with_lgl(g2))

# Lagged DCM
DCMLagG$new <- ifelse(DCMLagG$delay==-1,DCMLagG$y,NA)
DCMLagG$y <- ifelse(DCMLagG$delay==-1,DCMLagG$x,DCMLagG$y)
DCMLagG$x <- ifelse(DCMLagG$delay==-1,DCMLagG$new,DCMLagG$x)
DCMLagG$new <- NULL

gLag2 <- graph_from_data_frame(DCMLagG, directed=TRUE)
# plot(gLag2,vertex.size=4,vertex.label=NA,edge.width=0.1,layout=layout_with_lgl(gLag2),edge.arrow.size=0.2)

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


for (i in 1:nrow(taxSum)){
  t <- taxSum[i,]
  t2 <- taxSum[-i,]
  s <- subset(t2,Final==t$Final2 & Final2==t$Final & col==t$col)
  if (nrow(s)!=0){
    num <- as.numeric(rownames(s))
    x <- s$Final2
    s$Final2 <- s$Final
    s$Final <- x
    taxSum[num,] <- s
  }}

taxSum2 <- taxSum %>% group_by(Final,Final2,col) %>% tally(n)%>%as.data.frame()




taxSum2$Assoc <- paste(taxSum2$Final,taxSum2$Final2,sep="-")



ggplot(taxSum2,aes(x=reorder(Assoc, -n),y=n,fill=as.factor(col)))+geom_bar(stat="identity")+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12))+scale_fill_manual(name="Depth",values=c("darkgoldenrod3","indianred","navy"),labels=c("Surface & DCM","Surface Only","DCM Only"))+ylab("Number of Associations")+xlab("Association Type")
# ggsave("Figure5_Prelim.pdf",width=18,height=8)


### Cluster network ###
# Delete verticies without edges
int <- delete.vertices(igraph::simplify(int), degree(int)==0)
length(E(int))
length(V(int))

clus <- cluster_louvain(int)
length(clus)

V(int)$community <- clus$membership

# pdf("FigureX_Prelim.pdf")  
plot(int,vertex.label=NA,layout=layout_with_fr(int),vertex.size=3,mark.border="black",edge.width=0.3,mark.groups=communities(clus),vertex.color="grey",edge.color="black",mark.col=adjustcolor(c("coral1","deepskyblue2","darkseagreen1","mediumpurple"),alpha=0.4))
# dev.off()

### Community analysis ###
df <- read.csv(file.choose(),header=TRUE,row.names=1)

# Remove dups 
df <- subset(df,rownames(df)!="SPOTRe.DNA.115.02.16.12.5m_S31_L001_R1_trimmed.fastq" & rownames(df)!="SPOTRe.DNA.115.02.16.12.DCM_S32_L001_R1_trimmed.fastq")

df <- df[1:242,]

edgesClus <- data.frame(name=V(int)$name, clus=V(int)$community)

df <- mutate_all(df, function(x) as.numeric(as.character(x)))

df <- df/rowSums(df)
df <- as.data.frame(t(df))

df$name <- rownames(df)
edgesClus <- left_join(edgesClus,df)

# Fxn
dfClus <- edgesClus
rownames(dfClus) <- dfClus$name
dfClus$name <- NULL
clust <- dfClus$clus
dfClus$clus <- NULL

dfClus <- as.data.frame(t(dfClus))
dfClus$meanz <- rowMeans(dfClus)
namez <- colnames(dfClus)
dfClus <- as.matrix(dfClus)
sdz <- rowSds(dfClus)
dfClus <- as.data.frame(dfClus)
colnames(dfClus) <- namez
dfClus$sdz <- sdz

dfZ <- (dfClus - dfClus$meanz)/dfClus$sdz
dfZ$meanz <- NULL
dfZ$sdz <- NULL

dfZ <- as.data.frame(t(dfZ))
dfZ$clus <- clust
dfZ <- as.data.frame(t(dfZ))

# Environmental data
env <- read.csv(file.choose(),header=TRUE) #ctd
env2 <- read.csv(file.choose(),header=TRUE) # nut


colz <- colsplit(rownames(dfZ),"\\.",c("spot","dna","num","month","day","year","fin"))
colz <- colz[1:(nrow(colz)-1),]
l <- c(3,4,5,6,7,8,9)
colz$year <- ifelse(colz$year %in% l, paste(200,colz$year,sep=""), paste(20,colz$year,sep=""))
colz <- data.frame(Month=colz$month,Day=colz$day,Year=colz$year,Depth=colz$fin)
colz$Year <- as.integer(colz$Year)
colz$Day <- as.integer(colz$Day)
colz$Month <- as.integer(colz$Month)
colz$Depth <- ifelse(grepl("5m",colz$Depth),"5m","DCM")

env2 <- left_join(colz,env2)

colz2 <- colsplit(env$date,"/",c("month","day","year"))
env$Year <- colz2$year
env$Day <- colz2$day
l <- c(3,4,5,6,7,8,9)
env$Year <- ifelse(env$Year %in% l, paste(200,env$Year,sep=""), paste(20,env$Year,sep=""))
env$Year <- as.numeric(env$Year)
env <- env[c(4,8:12,14,19:20)]

env <- left_join(colz,env)
env <- env %>% distinct(Month,Day,Year,Depth,.keep_all = TRUE) %>% as.data.frame()

env2 <- env2[c(1:4,6:11,14)]

envTotal <- left_join(env,env2)
envTotal <- envTotal[c(5:16)]
envTotal <- mutate_all(envTotal, function(x) as.numeric(as.character(x)))
envTotal <- missForest(envTotal)
envTotal <- envTotal$ximp

dfZ <- as.data.frame(t(dfZ))

plotFxn <- function(df,cluster,env,col){
  sub <- subset(df, clus==cluster)
  sub$clus <- NULL
  sub <- as.data.frame(t(sub))
  sub$mean <- rowMeans(sub)
  newDf <- data.frame(z=sub$mean, envTotal[env])
  colnames(newDf) <- c("z","env")
  p <- ggplot(newDf, aes(x=env,y=z))+geom_point(fill=col,shape=21)+theme_classic()+geom_smooth(method="lm",formula='y~x',se=FALSE,color="black")
  t <- lm(newDf$z~newDf$env)
  summary(t)
  return(p)}

clus1Temp <- plotFxn(dfZ,1,4,"coral1")+xlab("Temperature (°C)") +ylab("Mean community z-score")+ggtitle("Community #1 vs. Temperature")+xlim(32.5,34)
clus1Temp

clus1Chl <- plotFxn(dfZ,1,12,"coral1")+xlab(expression("Chlorophyll"~italic(a)~"("*mu*"g L"^{"-1"}*")")) +ylab("Mean community z-score")+ggtitle(expression("Community #1 vs. Chlorophyll"~italic(a)))
clus1Chl

clus2Temp <- plotFxn(dfZ,2,4,"deepskyblue2")+xlab("Temperature (°C)") +ylab("Mean community z-score")+ggtitle("Community #2 vs. Temperature")+xlim(32.5,34)
clus2Temp

clus2Chl <- plotFxn(dfZ,2,12,"deepskyblue2")+xlab(expression("Chlorophyll"~italic(a)~"("*mu*"g L"^{"-1"}*")")) +ylab("Mean community z-score")+ggtitle(expression("Community #2 vs. Chlorophyll"~italic(a)))
clus2Chl

clus3Temp <- plotFxn(dfZ,3,4,"darkseagreen1")+xlab("Temperature (°C)") +ylab("Mean community z-score")+ggtitle("Community #3 vs. Temperature")+xlim(32.5,34)
clus3Temp

clus3Chl <- plotFxn(dfZ,3,12,"darkseagreen1")+xlab(expression("Chlorophyll"~italic(a)~"("*mu*"g L"^{"-1"}*")")) +ylab("Mean community z-score")+ggtitle(expression("Community #3 vs. Chlorophyll"~italic(a)))
clus3Chl

clus4Temp <- plotFxn(dfZ,4,4,"mediumpurple")+xlab("Temperature (°C)") +ylab("Mean community z-score")+ggtitle("Community #4 vs. Temperature")+xlim(32.5,34)
clus4Temp

clus4Chl <- plotFxn(dfZ,4,12,"mediumpurple")+xlab(expression("Chlorophyll"~italic(a)~"("*mu*"g L"^{"-1"}*")")) +ylab("Mean community z-score")+ggtitle(expression("Community #4 vs. Chlorophyll"~italic(a)))
clus4Chl
ggarrange(clus1Temp,clus2Temp,clus3Temp,clus4Temp,clus1Chl,clus2Chl,clus3Chl,clus4Chl,nrow=2,ncol=4)
ggsave("Community_Analysis_SurfDCM.pdf",width=13,height=8)
