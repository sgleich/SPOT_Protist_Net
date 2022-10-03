### PIDA ###
pida <- read.csv(file.choose(),header=TRUE)

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

tax2 <- tax[c(3:13)]
colnames(tax2)<- c("Kingdom","Supergroup","Divison","Class","Order","Family","Genus","Species","V1","Final","Low")

intEdge <- as.data.frame(get.edgelist(int))
intEdge <- left_join(intEdge,tax2)
colnames(tax2)<- c("Kingdom2","Supergroup2","Divison2","Class2","Order2","Family2","Genus2","Species2","V2","Final2","Low2")
intEdge <- left_join(intEdge,tax2)

intEdge2 <- subset(intEdge, Class !="Unknown")

out <- NULL
for (i in 1:nrow(pida)){
  tmp <- pida[i, ]
  tmpA <- paste(tmp$Taxonomic.level.3..org1,tmp$Genus.org1,tmp$Species.org1,sep="_")
  tmpB <- paste(tmp$Taxonomic.level.3..org2,tmp$Genus.org2,tmp$Species.org2,sep="_")
  
  # Genus                        
  sub <- intEdge2 %>% filter(str_detect(paste(tmpA),Genus) & str_detect(paste(tmpB),Genus2) | str_detect(paste(tmpA),Genus2) & str_detect(paste(tmpB),Genus))
  # Family
  if (nrow(sub)==0){
    sub <- intEdge2 %>% filter(str_detect(paste(tmpA),Family) & str_detect(paste(tmpB),Family2) | str_detect(paste(tmpA),Family2) & str_detect(paste(tmpB),Family))}
  # Order
  if (nrow(sub)==0){
    sub <- intEdge2 %>% filter(str_detect(paste(tmpA),Order) & str_detect(paste(tmpB),Order2) | str_detect(paste(tmpA),Order2) & str_detect(paste(tmpB),Order))
  }
  # Class
  if (nrow(sub)==0){
    sub <- intEdge2 %>% filter(str_detect(paste(tmpA),Class) & str_detect(paste(tmpB),Class2) | str_detect(paste(tmpA),Class2) & str_detect(paste(tmpB),Class))
  }
  if (nrow(sub)!=0){
    sub$num <- i
    out <- rbind(out,sub)
  }}
