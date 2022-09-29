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

### Unite Surface and DCM Delay 0 graphs ###
int <- intersection(g,g2)

plot(int,vertex.size=4,vertex.label=NA,edge.width=0.3,layout=layout_with_fr(int),vertex.color="grey")

### Union graphs with same SCC direction ###
surfNoLagG$fin <- paste(surfNoLagG$x,surfNoLagG$y,sep="_") 
DCMNoLagG$fin1 <-paste(DCMNoLagG$x,DCMNoLagG$y,sep="_") 
DCMNoLagG$fin2 <-paste(DCMNoLagG$y,DCMNoLagG$x,sep="_") 

combine <- subset(surfNoLagG, fin %in% DCMNoLagG$fin1 | fin %in% DCMNoLagG$fin2)

out <- NULL
for (i in 1:nrow(combine)){
  tmp <- combine[i, ]
  tmp2 <- subset(DCMNoLagG, fin1==tmp$fin | fin2==tmp$fin)
  if (tmp$scc > 0 & tmp2$scc > 0){
    out <- rbind(out,tmp)
  }
  if (tmp$scc < 0 & tmp2$scc < 0){
    out <- rbind(out,tmp)
  }
}

### Make network graph with species level labels ###
tax <- read.csv(file.choose(),header=TRUE)

assignTax <- function(df){
  df$Kingdom <- as.character(df$Kingdom)
  df$Supergroup <- as.character(df$Supergroup)
  df$Division<- as.character(df$Division)
  df$Class <- as.character(df$Class)
  df$Order <- as.character(df$Order)
  df$Family <- as.character(df$Family)
  df$Genus <- as.character(df$Genus)
  df$Species <- as.character(df$Species)
  df[is.na(df)] <- "Unknown"
  df$Final <- "Other Eukaryote"
  df$Final <- ifelse(df$Supergroup=="Rhizaria","Rhizaria",df$Final)
  df$Final <- ifelse(df$Division=="Chlorophyta","Chlorophyte",df$Final)
  df$Final <- ifelse(df$Division=="Ciliophora","Ciliate",df$Final)
  df$Final <- ifelse(df$Division=="Haptophyta","Haptophyte",df$Final)
  df$Final <- ifelse(df$Division=="Cryptophyta","Cryptophyte",df$Final)
  df$Final <- ifelse(df$Division=="Fungi","Fungi",df$Final)
  df$Final <- ifelse(df$Class=="Syndiniales","Syndiniales",df$Final)
  df$Final <- ifelse(df$Class=="Bacillariophyta","Diatom",df$Final)
  df$Final <- ifelse(df$Class!="Syndiniales" & df$Division=="Dinoflagellata","Dinoflagellate",df$Final)
  df$Final <- ifelse(grepl("MAST",df$Class),"MAST",df$Final)
  df$Final <- ifelse(df$Supergroup=="Stramenopiles" & df$Final !="Diatom" & df$Final !="MAST","Other Stramenopile",df$Final)
  df$Final <- ifelse(df$Supergroup=="Alveolata" & df$Final !="Dinoflagellate" & df$Final !="Syndiniales" & df$Final !="Ciliate","Other Alveolate",df$Final)
  df$Final <- ifelse(df$Supergroup=="Archaeplastida" & df$Final !="Chlorophyte","Other Archaeplastida",df$Final)
  df$Final <- ifelse(df$Kingdom=="Eukaryota" & df$Supergroup=="Unknown","Unknown Eukaryote",df$Final)
  return(df)
}

# tax <- assignTax(tax)

# taxOut <- NULL
# for (row in 1:nrow(tax)){
  # tmp <- tax[row, ]
  # tmp$low <- ifelse(tmp$Supergroup=="Unknown",tmp$Kingdom,NA)
  # tmp$low <- ifelse(tmp$Division=="Unknown" & is.na(tmp$low),tmp$Supergroup,tmp$low)
  # tmp$low <- ifelse(tmp$Class=="Unknown" & is.na(tmp$low),tmp$Division,tmp$low)
  # tmp$low <- ifelse(tmp$Order=="Unknown" & is.na(tmp$low),tmp$Class,tmp$low)
  # tmp$low <- ifelse(tmp$Family=="Unknown" & is.na(tmp$low),tmp$Order,tmp$low)
  # tmp$low <- ifelse(tmp$Genus=="Unknown" & is.na(tmp$low),tmp$Family,tmp$low)
  # tmp$low <- ifelse(tmp$Species=="Unknown" & is.na(tmp$low),tmp$Genus,tmp$low)
  # tmp$low <- ifelse(tmp$Species!="Unknown" & is.na(tmp$low),tmp$Species,tmp$low)
  # taxOut <- rbind(taxOut,tmp)}

# taxSum <- taxOut %>% group_by(low)%>%tally()

# final <- NULL
# for (row in 1:nrow(taxOut)){
  # t <- taxOut[row, ]
  # c <- subset(taxSum, low==t$low)
  # if (c$n >1){
    # y <- c$n
    # taxSum$n <- ifelse(taxSum$low==t$low,taxSum$n-1,taxSum$n)
    # t$low <- paste(t$low,y,sep="_")
  # }
  # final <- rtebind(final, t)
# }



### Crazy Circos-type plot ###
# Let's just plot it for Rhizaria assoc because otherwise it's uninterpretable
rhizDf <- subset(out,scc < 0)
# rhizDf <- subset(surfNoLagG,x %in% rhiz$num | y %in% rhiz$num)
rhizG <-graph_from_data_frame(rhizDf,directed=FALSE)
out <- graph_from_data_frame(rhizDf)
rhizG <- out
rhizDf2 <- data.frame(num=c(V(rhizG)$name))
rhizDf2 <- left_join(rhizDf2,tax)

rhizSum <- rhizDf2 %>% group_by(Final) %>% tally()

chlor <- subset(rhizDf2,Final=="Chlorophyte")
cil <- subset(rhizDf2,Final=="Ciliate")
cry <- subset(rhizDf2,Final=="Cryptophyte")
dia <- subset(rhizDf2,Final=="Diatom")
dino <- subset(rhizDf2,Final=="Dinoflagellate")
hapto <- subset(rhizDf2,Final=="Haptophyte")
mast <- subset(rhizDf2,Final=="MAST")
oeuk <- subset(rhizDf2,Final=="Other Eukaryote")
ostr <- subset(rhizDf2,Final=="Other Stramenopile")
rhiz <- subset(rhizDf2,Final=="Rhizaria")
syn <- subset(rhizDf2,Final=="Syndiniales")
unk <- subset(rhizDf2,Final=="Unknown Eukaryote")

v <- c(rhizSum$Final,chlor$num,cil$num,cry$num,dia$num,dino$num,hapto$num,mast$num,oeuk$num,ostr$num,rhiz$num,syn$num,unk$num)


# Hierarchy df
hierarchy <- data.frame(from=c(rep("Origin",10),rep("Chlorophyte",5),rep("Ciliate",8),rep("Cryptophyte",4),rep("Diatom",4),rep("Dinoflagellate",23),rep("Haptophyte",12),rep("MAST",1),rep("Other Eukaryote",9),rep("Other Stramenopile",6),rep("Syndiniales",20)),to=c(v))

# Vertices df
vertices <- data.frame(name=c("Origin",hierarchy$to),value=1, group=c(NA,rep("Origin",10),rep("Chlorophyte",5),rep("Ciliate",8),rep("Cryptophyte",4),rep("Diatom",4),rep("Dinoflagellate",23),rep("Haptophyte",12),rep("MAST",1),rep("Other Eukaryote",9),rep("Other Stramenopile",6),rep("Syndiniales",20))) 

colnames(vertices) <- c("num","value","group")
tmp <- data.frame(tax$num,tax$low)
colnames(tmp) <- c("num","low")
vertices <- left_join(vertices,tmp)
vertices$num <- ifelse(!is.na(vertices$low),vertices$low,vertices$num)
vertices$low <- NULL
colnames(vertices)[1] <- "name"

colnames(hierarchy)[2] <- "num"
hierarchy <- left_join(hierarchy,tmp)
hierarchy$low <- ifelse(is.na(hierarchy$low),hierarchy$num,hierarchy$low)
hierarchy$num <- NULL
colnames(hierarchy)[2] <- "to"


# Mygraph
mygraph <- graph_from_data_frame( hierarchy, vertices=vertices)

# pos <- subset(rhizDf, scc < 0)
t <- data.frame(get.edgelist(rhizG))
colnames(t)[1] <- "num"
t <- left_join(t,tmp)
colnames(tmp)[2] <- "low2"
colnames(t)<- c("var1","num","low")
t <- left_join(t,tmp)
t$scc <- E(rhizG)$scc
t$val <- ifelse(t$scc > 0, 1, -1)

connect <- data.frame(from=c(t$low),to=c(t$low2),value=c(t$val))
from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)

vertices$id <- NA
edges <- hierarchy
myleaves <- which(is.na( match(vertices$name, edges$from) ))
nleaves <- length(myleaves)
vertices$id[ myleaves ] <- seq(1:nleaves)
vertices$angle <- 90 - 360 * vertices$id / nleaves

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)

# flip angle BY to make them readable
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

mygraph <- igraph::graph_from_data_frame( edges, vertices=vertices )

# The connection object must refer to the ids of the leaves:
from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)


# Plot network
v <- randomcoloR::distinctColorPalette(12)

ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05)) +
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.5, colour="skyblue", width=0.9) +
  geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, label=name, angle = angle, hjust=hjust), size=2.5, alpha=1) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))+theme(legend.position = "right")+geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group),   size=5) +
  scale_colour_manual(name="Taxonomic Groups",values= c(v))+ggtitle("")
