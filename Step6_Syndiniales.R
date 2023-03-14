s5 <- read.csv(file.choose(),header=TRUE,row.names=1)
sDCM <- read.csv(file.choose(),header=TRUE,row.names=1)
colnames(s5) <- str_remove_all(colnames(s5),"^X")
colnames(sDCM) <- str_remove_all(colnames(sDCM),"^X")

s5 <- as.matrix(s5)
sDCM <- as.matrix(sDCM)

s5[s5<0] <- -1
s5[s5>0] <- 1
sDCM[sDCM<0] <- -1
sDCM[sDCM>0] <- 1


g5m <- graph_from_adjacency_matrix(s5,mode="undirected",weighted=TRUE)
gdcm <- graph_from_adjacency_matrix(sDCM,mode="undirected",weighted=TRUE)

int <- graph.intersection(g5m,gdcm)

interEdge <- as.data.frame(get.edgelist(int))
g5mEdge <- as.data.frame(get.edgelist(g5m))
g5mEdge$weight5 <- E(g5m)$weight
gdcmEdge <- as.data.frame(get.edgelist(gdcm))
gdcmEdge$weightDCM <- E(gdcm)$weight

interEdge$fin1 <- paste(interEdge$V1,interEdge$V2,sep="_")
g5mEdge$fin1 <- paste(g5mEdge$V1,g5mEdge$V2,sep="_")
g5mEdge$fin2 <- paste(g5mEdge$V2,g5mEdge$V1,sep="_")
gdcmEdge$fin1 <- paste(gdcmEdge$V1,gdcmEdge$V2,sep="_")
gdcmEdge$fin2 <- paste(gdcmEdge$V2,gdcmEdge$V1,sep="_")

interEdge <- left_join(interEdge,g5mEdge)
interEdge <- left_join(interEdge,gdcmEdge)

interEdge <- subset(interEdge,weight5==weightDCM)

t <- read.delim("taxonomy_90.tsv",header=TRUE)
interEdge2 <- interEdge[c(1:2)]
t2 <- t[c(1:2)]
colnames(t2) <- c("V1","Tax1")
interEdge2 <- left_join(interEdge2,t2)
colnames(t2) <- c("V2","Tax2")
interEdge2 <- left_join(interEdge2,t2)

genera <- c("Prorocentrum","Gyrodinium","Gonyaulax","Tripos","Heterocapsa","Gymnodinium","Karlodinium","Lepidodinium")
syn1 <- subset(interEdge2,grepl("Group-II-",interEdge2$Tax1) & grepl("Prorocentrum|Gyrodinium|Gonyaulax|Tripos|Heterocapsa|Gymnodinium|Karlodinium|Lepidodinium",interEdge2$Tax2))
syn2 <- subset(interEdge2,grepl("Group-II-",interEdge2$Tax2) & grepl("Prorocentrum|Gyrodinium|Gonyaulax|Tripos|Heterocapsa|Gymnodinium|Karlodinium|Lepidodinium",interEdge2$Tax1))
syn <- rbind(syn1,syn2)
syn$p <- ifelse(grepl("Prorocentrum",syn$Tax1)|grepl("Prorocentrum",syn$Tax2),"Prorocentrum",NA)
syn$p <- ifelse(grepl("Gyrodinium",syn$Tax1)|grepl("Gyrodinium",syn$Tax2),"Gyrodinium",syn$p)
syn$p <- ifelse(grepl("Gonyaulax",syn$Tax1)|grepl("Gonyaulax",syn$Tax2),"Gonyaulax",syn$p)
syn$p <- ifelse(grepl("Tripos",syn$Tax1)|grepl("Tripos",syn$Tax2),"Tripos",syn$p)
syn$p <- ifelse(grepl("Heterocapsa",syn$Tax1)|grepl("Heterocapsa",syn$Tax2),"Heterocapsa",syn$p)
syn$p <- ifelse(grepl("Gymnodinium",syn$Tax1)|grepl("Gymnodinium",syn$Tax2),"Gymnodinium",syn$p)
syn$p <- ifelse(grepl("Karlodinium",syn$Tax1)|grepl("Karlodinium",syn$Tax2),"Karlodinium",syn$p)
syn$p <- ifelse(grepl("Lepidodinium",syn$Tax1)|grepl("Lepidodinium",syn$Tax2),"Lepidodinium",syn$p)

tmp <- syn %>% filter(p=="Prorocentrum")

new <-  ape::read.tree(file.choose())
tax <- data.frame(Feature.ID=c(new$tip.label))
t <- read.delim("taxonomy_90.tsv",header=TRUE)
tax <- left_join(tax,t)
colz <- colsplit(tax$Taxon,";",c("s","k","p","c","o","f","g","sp"))
colz$new <- colz$f
colz$new <- ifelse(colz$new=="",colz$o,colz$new)
tax$New <- colz$new

colnames(tax)[1] <- "Feature.ID"


# new$colr <- tax$col
# tax$New <- paste(tax$Feature.ID,tax$New,sep="_")
new$tip.label<-tax$New 
new$node.label <- as.character(new$node.label)
s <- seq(1:100)
s <- s[50:100]
s <- as.character(s)
new$node.label <- ifelse(new$node.label %in% s,"",new$node.label)


pdf("../../TRYb.pdf",width=20,height=25)
plot(new)
# tiplabels(new$tip.label)
nodelabels(new$node.label,adj=c(-0.2),frame="none",cex=0.8)
dev.off()


#####
# Plot
counts <- read.delim(file.choose(),header=TRUE,row.names=1)
counts <- as.data.frame(t(counts))
manifest <- read.delim(file.choose(),header=FALSE,sep=",")

# We need to rename our samples in the counts dataframe so that they contain meaningful information
namez <- colsplit(manifest$V2,"-",c("a","b","c","d","e","f","g"))
namez$depth <- ifelse(grepl("5m",namez$g),"5m","DCM")
# Making a "fin" column that contains the date and depth of each sample
namez$fin <- paste(namez$c,namez$d,namez$e,namez$f,namez$depth,sep="_")

# Make a new dataframe with old sample name info and new sample name info
fin <- data.frame(namez=c(manifest$V1),new=c(namez$fin))
fin$new <- paste("SPOT",fin$new,sep="_")
fin <- fin %>% distinct(new,.keep_all = TRUE)

# Let's add the new (meaningful) sample names to our counts dataframe counts <- as.data.frame(t(counts))
counts$namez <- rownames(counts)
dfNew <- left_join(counts,fin)
rownames(dfNew) <- dfNew$new
dfNew$namez <- NULL
dfNew$new <- NULL

dfNew <- subset(dfNew,rownames(dfNew)!="SPOT_115_2_16_12_5m"& rownames(dfNew)!="SPOT_115_2_16_12_DCM")
dfCLR <- clr(dfNew)
dfCLR <- as.data.frame(dfCLR)

# Make a "date" column in ASV counts dataframe
datez <- colsplit(rownames(dfCLR),"_",c("spot","Cruise","m","d","y","Depth"))
env <- read.csv(file.choose(),header=TRUE)
datez <- left_join(datez,env)
datez2 <- colsplit(datez$Date,"/",c("month","day","year"))
dfCLR$month <- datez2$month
dfCLR$day <- datez2$day
dfCLR$year <- datez2$year
dfCLR$date <- as.Date(with(dfCLR, paste(year, month, day,sep="-")), "%y-%m-%d")

# Remove month, day, year columns since we have a date column now
dfCLR$month <- NULL
dfCLR$day <- NULL
dfCLR$year <- NULL

# Subset 5m and DCM samples
asvSurf <- subset(dfCLR,grepl("5m",rownames(dfCLR)))
asvDCM <- subset(dfCLR,grepl("DCM",rownames(dfCLR)))

asvSurf <- asvSurf %>% arrange(date) %>% as.data.frame() 
asvDCM<- asvDCM %>% arrange(date) %>% as.data.frame() 

rownames(asvSurf) <- asvSurf$date
rownames(asvDCM) <- asvDCM$date
asvSurf$date <- NULL
asvDCM$date <- NULL

asvSurf <- as.data.frame(t(asvSurf))
asvDCM <- as.data.frame(t(asvDCM))

net <- read.csv(file.choose(),header=TRUE,row.names=1)
t <- read.delim("taxonomy_90.tsv",header=TRUE)
net <- data.frame(Feature.ID=rownames(net))
net <- left_join(net,t)

group2 <- subset(net,grepl("Dino-Group-II-",net$Taxon))
dinos <- subset(net,grepl("Dinophyceae",net$Taxon))

asvSurfSub <- subset(asvSurf,rownames(asvSurf) %in% group2$Feature.ID |rownames(asvSurf) %in% dinos$Feature.ID)
asvSurfSub$group <- ifelse(rownames(asvSurfSub) %in% group2$Feature.ID,"Group-II Syndiniales","Dinoflagellate")
asvSurfSum <- asvSurfSub %>% group_by(group)%>% summarize_all(mean) %>% as.data.frame()
asvSurfSum <- melt(asvSurfSum,id.vars=c("group"))
asvSurfSum <- as.data.frame(t(asvSurfSum))
p1 <- ggplot(asvSurfSum,aes(x=variable,y=value,color=group,group=group))+geom_point()+geom_line()+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=4.7))+scale_color_manual(name="Taxonomic Group",values=c("grey20","grey60"))+xlab("Date of Sampling")+ylab("Mean CLR-transformed Abundance")

datez <- colsplit(asvSurfSum$variable,"-",c("Year","Month","day"))
asvSurfSum$Month <- datez$Month
asvSurfSum2 <- asvSurfSum
asvSurfSum2 <- asvSurfSum2 %>% pivot_wider(id_cols=c("Month","variable"),names_from = c("group"),values_from = ("value"))

p2 <- ggplot(asvSurfSum2,aes(x=Dinoflagellate,y=`Group-II Syndiniales`))+geom_point(aes(fill=as.factor(Month),shape=as.factor(Month)),size=2,color="black")+scale_fill_manual(name="Month",values=c("dodgerblue","dodgerblue","darkolivegreen4","darkolivegreen4","darkolivegreen4","goldenrod1","goldenrod1","goldenrod1","firebrick2","firebrick2","firebrick2","dodgerblue"),labels=c("January","February","March","April","May","June","July","August","September","October","November","December"))+scale_shape_manual(name="Month",values=c(21,22,21,22,24,21,22,24,21,22,24,24),labels=c("January","February","March","April","May","June","July","August","September","October","November","December"))+theme_classic()+xlab("Mean CLR-transformed Abundance Dinoflagellates")+ylab("Mean CLR-transformed Abundance Group-II Syndiniales")+stat_poly_line(se=FALSE,color="black")+stat_poly_eq()

p1+p2+ plot_layout(ncol = 1)
ggsave("../../syndiniales_surf.pdf",width=13,height=16)

### DCM

asvDCMSub <- subset(asvDCM,rownames(asvDCM) %in% group2$Feature.ID |rownames(asvDCM) %in% dinos$Feature.ID)
asvDCMSub$group <- ifelse(rownames(asvDCMSub) %in% group2$Feature.ID,"Group-II Syndiniales","Dinoflagellate")
asvDCMSum <- asvDCMSub %>% group_by(group)%>% summarize_all(mean) %>% as.data.frame()
asvDCMSum <- melt(asvDCMSum,id.vars=c("group"))

p1 <- ggplot(asvDCMSum,aes(x=variable,y=value,color=group,group=group))+geom_point()+geom_line()+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=4.7))+scale_color_manual(name="Taxonomic Group",values=c("grey20","grey60"))+xlab("Date of Sampling")+ylab("Mean CLR-transformed Abundance")

datez <- colsplit(asvDCMSum$variable,"-",c("Year","Month","day"))
asvDCMSum$Month <- datez$Month
asvDCMSum2 <- asvDCMSum
asvDCMSum2 <- asvDCMSum2 %>% pivot_wider(id_cols=c("Month","variable"),names_from = c("group"),values_from = ("value"))

p2 <- ggplot(asvDCMSum2,aes(x=Dinoflagellate,y=`Group-II Syndiniales`))+geom_point(aes(fill=as.factor(Month),shape=as.factor(Month)),size=2,color="black")+scale_fill_manual(name="Month",values=c("dodgerblue","dodgerblue","darkolivegreen4","darkolivegreen4","darkolivegreen4","goldenrod1","goldenrod1","goldenrod1","firebrick2","firebrick2","firebrick2","dodgerblue"),labels=c("January","February","March","April","May","June","July","August","September","October","November","December"))+scale_shape_manual(name="Month",values=c(21,22,21,22,24,21,22,24,21,22,24,24),labels=c("January","February","March","April","May","June","July","August","September","October","November","December"))+theme_classic()+xlab("Mean CLR-transformed Abundance Dinoflagellates")+ylab("Mean CLR-transformed Abundance Group-II Syndiniales")+stat_poly_line(se=FALSE,color="black")+stat_poly_eq()

p1+p2+ plot_layout(ncol = 1)
ggsave("../../syndiniales_DCM2.pdf",width=13,height=16)
