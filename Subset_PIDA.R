### PIDA ###
pida <- read.csv(file.choose(),header=TRUE)

### Load eLSA Output and Filter ###
# Load data
surf <- read.delim(file.choose(),header=TRUE)
dcm <- read.delim(file.choose(),header=TRUE)

# Subset based on p and q-value
surfSig <- subset(surf,Psscc<0.01 & Qsscc < 0.01)
dcmSig <- subset(dcm,Psscc<0.01 & Qsscc < 0.01)
surfG <- data.frame(x=c(surfSig$X),y=c(surfSig$Y),scc=c(surfSig$SCC),delay=c(surfSig$Delay))
dcmG <- data.frame(x=c(dcmSig$X),y=c(dcmSig$Y),scc=c(dcmSig$SCC),delay=c(dcmSig$Delay))

surfG$s <- ifelse(surfG$delay==-1,surfG$x,NA)
surfG$x <- ifelse(!is.na(surfG$s),surfG$y,surfG$x)
surfG$y <- ifelse(!is.na(surfG$s),surfG$s,surfG$y)
surfG$delay <- ifelse(!is.na(surfG$s),1,surfG$delay)
unique(surfG$delay)
 
### Plot lagged and unlagged networks surface ###
# Unlagged surface
surfG$scc <- NULL
surfG$s <- NULL
g <- graph_from_data_frame(surfG,directed=TRUE)

# pdf("SurfaceNet.pdf")
plot(g,vertex.size=4,vertex.label=NA,edge.width=0.1,layout=layout_with_lgl(g),edge.arrow.mode = ifelse(E(g)$delay==0, "-", ">"),edge.arrow.size=0.2,vertex.color="grey20")
# dev.off()
          
dcmG$s <- ifelse(dcmG$delay==-1,dcmG$x,NA)
dcmG$s <- ifelse(dcmG$delay==-1,dcmG$x,NA)
dcmG$x <- ifelse(!is.na(dcmG$s),dcmG$y,dcmG$x)
dcmG$y <- ifelse(!is.na(dcmG$s),dcmG$s,dcmG$y)
dcmG$delay <- ifelse(!is.na(dcmG$s),1,dcmG$delay)

# Plot surface and DCM networks with lagged and unlagged interactions
surfG$scc <- NULL
# Plot surface and DCM networks with lagged and unlagged interactions
dcmG$scc <- NULL
dcmG$s <- NULL
g2 <- graph_from_data_frame(dcmG,directed=TRUE)

# pdf("DCMNet.pdf")
plot(g2,vertex.size=4,vertex.label=NA,edge.width=0.1,layout=layout_with_lgl(g2),edge.arrow.mode = ifelse(E(g2)$delay==0, "-", ">"),edge.arrow.size=0.2,vertex.color="grey20")
# dev.off()
           
### Unite Surface and DCM Delay 0 graphs ###
int <- intersection(g,g2)

plot(int,vertex.size=4,vertex.label=NA,edge.width=0.3,layout=layout_with_fr(int),vertex.color="grey",edge.arrow.mode = ifelse(E(g2)$delay==0, "-", ">"),edge.arrow.size=0.2,vertex.color="grey20")

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
