######################################################
t2 <- read.csv(file.choose(),header=TRUE,row.names=1)
t2 <- data.frame(name=c(rownames(t2)))
t2$name <- str_remove_all(t2$name,"S_")

colnames(tax) <- c("name","tax")

t2 <- left_join(t2,tax)

syn <- subset(t2,grepl("Syndinia",t2$tax))

for (line in 1:nrow(full)){
  l <- full[line,]
  if (grepl("Syndinial",l$Tax1)){
    s <- subset(full,V1==l$V1)
    print(s$V1)
    print(nrow(s))
  }
  if (grepl("Syndinial",l$Tax2)){
    s <- subset(full,V2==l$V2)
    print(s$V2)
    print(nrow(s))
  }
  
}

try <- NULL
for (i in 1:100){
  t <- erdos.renyi.game(389, 2237,"gnm")
  V(t)$deg <- degree(t)
  mDeg <- mean(V(t)$deg)
  mPath <- mean_distance(t)
  mClus <- transitivity(t,type="global")
  vec <- c(i,mDeg,mPath,mClus) 
  try <- rbind(try,vec)}
try <- as.data.frame(try)
######################################################
counts <- read.delim("./SPOT/export_dir/feature-table.tsv",header=TRUE,row.names=1)
manifest <- read.delim("./SPOT/export_dir/manifest.txt",header=FALSE,sep=",")

# We need to rename our samples in the counts dataframe so that they contain meaningful information
namez <- colsplit(manifest$V2,"-",c("a","b","c","d","e","f","g"))
namez$depth <- ifelse(grepl("5m",namez$g),"5m","DCM")
# Making a "fin" column that contains the date and depth of each sample
namez$fin <- paste(namez$c,namez$d,namez$e,namez$f,namez$depth,sep="_")

# Make a new dataframe with old sample name info and new sample name info
fin <- data.frame(namez=c(manifest$V1),new=c(namez$fin))
fin$new <- paste("SPOT",fin$new,sep="_")
fin <- fin %>% distinct(new,.keep_all = TRUE)

# Let's add the new (meaningful) sample names to our counts dataframe 
counts <- as.data.frame(t(counts))
counts$namez <- rownames(counts)
dfNew <- left_join(counts,fin)
rownames(dfNew) <- dfNew$new
dfNew$namez <- NULL
dfNew$new <- NULL

# Network ASVs
surf <- read.csv(file.choose(),header=TRUE,row.names=1)
df <- data.frame(Feature.ID=c(rownames(surf)))
df$Feature.ID <- str_remove_all(df$Feature.ID,"S_")
tax <- read.delim("./SPOT/export_dir/taxonomy_90.tsv",header=TRUE)
tax$Confidence <- NULL

df <- left_join(df,tax)

syn2 <- subset(df,grepl("Group-II;",df$Taxon))
syn1 <- subset(df,grepl("Group-I;",df$Taxon))
dino <- subset(df,grepl("Dinophyceae",df$Taxon))
cil <- subset(df,grepl("Ciliopho",df$Taxon))

dfNew <- as.data.frame(clr(dfNew))
dfNew <- as.data.frame(t(dfNew))

syn2C <- subset(dfNew,rownames(dfNew) %in% syn2$Feature.ID)
syn1C <- subset(dfNew,rownames(dfNew) %in% syn1$Feature.ID)
dinoC <- subset(dfNew,rownames(dfNew) %in% dino$Feature.ID)
cilC <- subset(dfNew,rownames(dfNew) %in% cil$Feature.ID)

final <- data.frame(dino=c(colMeans(dinoC)),syn1=c(colMeans(syn1C)),syn2=c(colMeans(syn2C)),cil=c(colMeans(cilC)))

colz <- colsplit(rownames(final),"_",c("SPOT","Cruise","m","y","d","Depth"))
env <- read.csv(file.choose(),header=TRUE,row.names=1)
colz <- left_join(colz,env)

final$month <- colz$Month
final$depth <- colz$Depth
final$month <- ifelse(is.na(final$month),2,final$month)
final$cruise <- colz$Cruise
final$date <- colz$Date


fun_color_range <- colorRampPalette(c("darkgreen","white","darkgoldenrod3")) 
my_colors <- fun_color_range(12)   


dino <- final %>% ggplot(aes(x=dino,y=syn2))+geom_point(aes(fill=as.factor(month)),shape=21)+theme_classic()+xlab("Mean CLR-transformed Abundance Dinoflagellate ASVs")+ylab("Mean CLR-transformed Abundance Group-II Syndiniales ASVs")+scale_fill_manual(name="Month",values=c(my_colors),labels=c("January","February","March","April","May","June","July","August","September","October","November","December"))+geom_smooth(formula="y~x",method="lm",color="black")+ggtitle("Group-II Syndiniales vs. Dinoflagellates")

dino+cil+plot_layout(guides = "collect")

new <- data.frame(cil=c(final$cil),syn=c(final$syn1),cruise=c(final$cruise),date=c(final$date),depth=c(final$depth))
new <- melt(new,id.vars=c("cruise","date","depth"))
new5 <- subset(new,depth=="5m")
new5 <- subset(new5,cruise!=115)
new5 <- new5 %>% arrange(cruise)%>%as.data.frame()
new5$cruise2 <- rep(1:122,each=2)
new5$cruise2 <- rep(1:120,each=2)

p5a <- new5 %>% filter(depth=="5m") %>% ggplot(aes(x=cruise2,y=value,color=variable))+geom_point()+geom_line()+scale_x_continuous(breaks=new5$cruise2, labels=new5$date,expand = expansion(add = 5))+theme_classic(base_size=14)+xlab("Date of Sampling")+ylab("Mean CLR-transformed Abundance")+scale_color_manual(name="Taxonomic Group",values=c("grey20","grey80"),,labels=c("Ciliate","Group-I Syndiniales"))+theme(axis.text.x = element_text(angle = 90, vjust=0.5,hjust=1,size=6))+ggtitle("Surface (5 m)")
ggsave("see.pdf",width=12,height=6)    

p5a <- p5a+ggtitle("Surface (5 m)")

ggarrange(p5a,p5b,ncol=1,nrow=2,common.legend = TRUE)
#######################################################

### SPOT Network Analysis ###
### Process Random Forest Output ###
### By: Samantha Gleich ###
### Last Updated: 1/12/23 ###

# Colors to use for all taxonomic analyses
taxCols <- c("#E1746D","#76C3D7","#DE5AB1","#D5E0AF","#DED3DC","#87EB58","#D4DC60","#88E5D3","#88AAE1","#DBA85C","#8B7DDA","#9A8D87","#D99CD1","#B649E3","#7EDD90")
names(taxCols) <- c("Chlorophyte","Ciliate","Cryptophyte","Diatom","Dinoflagellate","Fungi","Haptophyte","MAST","Other Alveolate","Other Archaeplastida","Other Eukaryote","Other Stramenopiles","Rhizaria","Syndiniales","Unknown Eukaryote")

rfOut <- read.csv("./SPOT/randomForest_Results_January2023.csv",header=TRUE,row.names=1)
y <- unique(rfOut$Y)
asv <- 0
env <- 0
both <- 0
asv_namez <- NULL
env_namez <- NULL
both_namez <- NULL

# Loop through results to see if each ASV was best predicted (lowest MSE) by ASVs, Environmental Variables, or ASVs + Environmental Variables
for(i in 1:length(unique(y))){
  subs <- subset(rfOut,Y==y[i])
  m <- min(subs$estimate)
  subs2 <- subset(subs,estimate==m)
  if (subs2[1,5]=="ASVs"){
    asv <- asv+1
    asv_namez <- c(asv_namez,subs[1,7])}
  if (subs2[1,5]=="Environment"){
    env <- env+1
    env_namez <- c(env_namez,subs[1,7])}
  if (subs2[1,5]=="ASV+ Environment"){
    both_namez <- c(both_namez,subs[1,7])
    both <- both+1}
}

# Let's clean up the output
fin <- data.frame(name=c(asv_namez)) 
fin$category <- "ASVs"
fin2 <- data.frame(name=c(env_namez)) 
fin2$category <- "Environment"
fin3 <- data.frame(name=c(both_namez)) 
fin3$category <- "ASV + Environment"
fin <- rbind(fin,fin2,fin3)

# Take a look at the taxonomic breakdown of the ASVs that were best predicted by each group of variables
tax <- read.delim("./SPOT/export_dir/taxonomy_90.tsv",header=TRUE)
tax$Confidence <- NULL
colnames(tax)[1] <- "name" 
fin <- left_join(fin,tax)

# Classify taxonomy
taxz <- colsplit(fin$Taxon,";",c("Supergroup","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxz$fin <- ifelse(taxz$Class=="Syndiniales","Syndiniales",NA)
taxz$fin <- ifelse(taxz$Class=="Dinophyceae","Dinoflagellate",taxz$fin)
taxz$fin <- ifelse(taxz$Class=="Bacillariophyta","Diatom",taxz$fin)
taxz$fin <- ifelse(grepl("MAST",taxz$Class),"MAST",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Rhizaria","Rhizaria",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Chlorophyta","Chlorophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Cryptophyta","Cryptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Haptophyta","Haptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Ciliophora","Ciliate",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Stramenopiles" & is.na(taxz$fin),"Other Stramenopiles",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Archaeplastida" & is.na(taxz$fin),"Other Archaeplastids",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Alveolata" & is.na(taxz$fin),"Other Alveolate",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin),"Other Eukaryote",taxz$fin)
unique(taxz$fin)

# Find # of ASVs in each taxonomic group that were best predicted by each predictor group
fin$tax <- taxz$fin
fin$name <- NULL
fin$Taxon <- NULL
finSum <- fin %>% group_by(tax,category) %>% tally()

# Get ready to plot
# colrs <- randomcoloR::distinctColorPalette(length(unique(fin$tax)))
finSum$category <- factor(finSum$category,levels=c("ASVs","ASV + Environment","Environment"))

ggplot(finSum,aes(x=category,y=n,fill=tax))+geom_bar(stat="identity",color="grey20")+scale_fill_manual(name="Taxonomic Group",values=c(taxCols))+theme_classic()+ylab("Number of ASVs")+xlab("Random Forest Predictors")+scale_x_discrete(guide = guide_axis(angle = 45),labels=c("ASVs Only","ASVs + Environmental Variables","Environmental Variables Only"))
ggsave("../RF_Output.pdf",width=7,height=6)

# Now we identify "clusters" of organisms that were best predicted by environmental variables
rfPerVar <- read.csv("./SPOT/randomForest_Results_Env_January2023.csv",header=TRUE,row.names=1)

subs <- subset(rfPerVar,name %in% env_namez)
subs <- left_join(subs,tax)

taxz <- colsplit(subs$Taxon,";",c("Supergroup","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxz$fin <- ifelse(taxz$Class=="Syndiniales","Syndiniales",NA)
taxz$fin <- ifelse(taxz$Class=="Dinophyceae","Dinoflagellate",taxz$fin)
taxz$fin <- ifelse(taxz$Class=="Bacillariophyta","Diatom",taxz$fin)
taxz$fin <- ifelse(grepl("MAST",taxz$Class),"MAST",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Rhizaria","Rhizaria",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Chlorophyta","Chlorophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Cryptophyta","Cryptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Haptophyta","Haptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Ciliophora","Ciliate",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Stramenopiles" & is.na(taxz$fin),"Other Stramenopiles",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Archaeplastida" & is.na(taxz$fin),"Other Archaeplastids",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Alveolata" & is.na(taxz$fin),"Other Alveolate",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin),"Other Eukaryote",taxz$fin)
unique(taxz$fin)

subs$taxFin <- taxz$fin
subsSum <- subs %>% group_by(taxFin,Param) %>% tally()

subsSum$Param <- factor(subsSum$Param,levels=c("MEI","DayOfYear","NO2.NO3","cyanos_PA","Chla","cyanos","O2Wink","sar11_PA","CTDBEAM","DayLength","NH4","PP","SST"))

ggplot(subsSum,aes(x=Param,y=n,fill=taxFin))+geom_bar(stat="identity",color="black")+theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust=1))+xlab("Environmental Predictors")+ylab("Number of ASVs")+scale_fill_manual(name="Taxonomic Group",values=c(taxCols))+scale_x_discrete(labels=c("MEI","Day of Year","Nitrate","Particle-associated cyanobacteria",expression("Chlorophyll"~italic(a)),"Free-living cyanobacteria","Oxygen (Winkler)","Particle-associated SAR11","Beam attenuation","Day Length","Ammonium","Average monthly PP (satellite)","Average monthly SST (satellite)"))+ylim(0,15)

ggsave("../FigureS1b.pdf",width=10,height=8)

### CLUSTER ANALYSIS??!! ###
# We will limit our cluster analysis to ASVs whose variability in abundances could be explained by enviornmental at > 30%
dfEnv <- subset(rfOut,Predictors=="Environment")
dfEnv <- subset(dfEnv, PercentVar > 30)
rfPerVar <- data.frame(ASV_ID=c(rfPerVar$name),Param=c(rfPerVar$Param))
j <- left_join(dfEnv,rfPerVar)

# Now let's load in raw ASV abundances so we can visualize the abundances of the different clusters as a function of environmental covariates
counts <- read.delim("feature-table.tsv",header=TRUE,row.names=1)
manifest <- read.delim("manifest.txt",header=FALSE,sep=",")

# We need to rename our samples in the counts dataframe so that they contain meaningful information
namez <- colsplit(manifest$V2,"-",c("a","b","c","d","e","f","g"))
namez$depth <- ifelse(grepl("5m",namez$g),"5m","DCM")
# Making a "fin" column that contains the date and depth of each sample
namez$fin <- paste(namez$c,namez$d,namez$e,namez$f,namez$depth,sep="_")

# Make a new dataframe with old sample name info and new sample name info
fin <- data.frame(namez=c(manifest$V1),new=c(namez$fin))
fin$new <- paste("SPOT",fin$new,sep="_")
fin <- fin %>% distinct(new,.keep_all = TRUE)

# Let's add the new (meaningful) sample names to our counts dataframe 
counts <- as.data.frame(t(counts))
counts$namez <- rownames(counts)
dfNew <- left_join(counts,fin)
rownames(dfNew) <- dfNew$new
dfNew$namez <- NULL
dfNew$new <- NULL
dfNew<- subset(dfNew,rownames(dfNew)!="SPOT_115_2_16_12_5m"& rownames(dfNew)!="SPOT_115_2_16_12_DCM")

# Calculate relative abundances
dfNew <- dfNew/rowSums(dfNew)
dfNew <- as.data.frame(t(dfNew))

# Remove ASVs that are not in environmental clusters
namez <- j$ASV_ID
dfKeep <- subset(dfNew,rownames(dfNew) %in% namez)

# Calculate z-scores
dfKeep <- as.data.frame(t(dfKeep))
dfZscore <- data.frame(x=rep(1,242))
for (col in 1:ncol(dfKeep)){
  x <- dfKeep[c(col)]
  name <- colnames(x)
  colnames(x) <- "V1"
  meanz <- mean(x$V1)
  sdz <- sd(x$V1)
  new <- ((x$V1)-meanz)/sdz
  new <- as.data.frame(new)
  colnames(new) <- name
  dfZscore <- cbind(dfZscore,new)
}
dfZscore$x <- NULL
rownames(dfZscore) <- rownames(dfKeep)

# Add information about the community each ASV is associated with
dfZscore <- as.data.frame(t(dfZscore))
dfZscore$name <- rownames(dfZscore)
communities <- data.frame(name=c(j$ASV_ID),community=c(j$Param))
dfZscore <- left_join(dfZscore,communities)
rownames(dfZscore) <- dfZscore$name
dfZscore$name <- NULL
rownames(dfZscore) <- paste(dfZscore$community,rownames(dfZscore),sep="_")
dfZscore$community <- NULL

# Add environmental data to z-score transformed ASV abundances 
env <- read.csv("../SPOT_Env_NewJan11.csv",header=TRUE)
env$Month <- NULL
env$Date <- NULL

dfZscore <- as.data.frame(t(dfZscore))
colz <- colsplit(rownames(dfZscore),"_",c("spot","Cruise","Month","Day","Year","Depth"))
dfZscore$Depth <- colz$Depth
dfZscore$Cruise <- colz$Cruise
all <- left_join(dfZscore,env)

# Plot communities

j2 <- j %>% group_by(Param) %>% tally() %>% arrange(desc(n))
# colz <- randomcoloR::distinctColorPalette(10)
ggplot(j2,aes(x=reorder(Param,-n),y=n))+geom_bar(stat="identity",color="black")+theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust=1))+geom_text(aes(label=n),vjust=-0.5)+scale_x_discrete(labels=c("MEI","Temperature","Free-living cyanobacteria","Particle-associated SAR11","Sea surface height","Day of year","Average monthly SST (satellite)",expression("Average monthly Chlorophyll"~italic(a)~"(satellite)"),"Day length","Ammonium",expression("Chlorophyll"~italic(a)),"Salinity","Particle-associated cyanobacteria","Oxygen (Winkler)","Average monthly primary production (satellite)"))+ylim(0,21)+xlab("Environmental Predictors")+ylab("Number of ASVs")
ggsave("../FigureSX2.pdf",width=8,height=6)

colrs <- distinctColorPalette(10)

# TEMPERATURE
subsTmp <- subset(all,select=c(grepl("CTDTMP",colnames(all))))
subsTmp$CTDTMP <- NULL
subsTmp <- as.data.frame(t(subsTmp))
subsTmp <- subsTmp %>% summarize_all(mean) %>% as.data.frame()
subsTmp <- data.frame(t(subsTmp))
envB <- all[c(62:ncol(all))]
finTmp <- cbind(subsTmp,envB)
finTmp <- subset(finTmp,!is.na(CTDTMP))
pTmp <- ggplot(finTmp,aes(x=CTDTMP,y=t.subsTmp.))+geom_point(shape=21,fill=colrs[1],color="black")+theme_classic()+ggtitle("Temperature Community (# ASVs = 3)")+xlab("Temperature (°C)")+ylab("Mean Community Z-score")+xlim(0,25)
pTmp

# CYANOS
subsCyan <- subset(all,select=c(grepl("cyano",colnames(all))))
subsCyan <- subset(subsCyan,select=c(!grepl("cyanos_PA",colnames(subsCyan))))
subsCyan <- as.data.frame(t(subsCyan))
subsCyan <- subsCyan %>% summarize_all(mean) %>% as.data.frame()
subsCyan <- data.frame(t(subsCyan))
envB <- all[c(24:ncol(all))]
finCyan <- cbind(subsCyan,envB)
finCyan<- subset(finCyan,!is.na(cyanos))
pCyan <- ggplot(finCyan,aes(x=cyanos,y=t.subsCyan.))+geom_point(shape=21,fill=colz[2],color="black")+theme_classic()+ggtitle("Cyanobacterial Community (# ASVs = 5)")+xlab("Relative Abundance of Free-living Cyanobacteria")+ylab("Mean Community Z-score")
pCyan

# Day of Year
subsDay <- subset(all,select=c(grepl("DayOfYear",colnames(all))))
subsDay <- as.data.frame(t(subsDay))
subsDay <- subsDay %>% summarize_all(mean) %>% as.data.frame()
subsDay <- data.frame(t(subsDay))
envB <- all[c(24:ncol(all))]
finDay <- cbind(subsDay,envB)
finDay<- subset(finDay,!is.na(DayOfYear))
pDay <- ggplot(finDay,aes(x=DayOfYear,y=t.subsDay.))+geom_point(shape=21,fill=colz[3],color="black")+theme_classic()+ggtitle("Day of Year Community (# ASVs = 4)")+xlab("Day of Year")+ylab("Mean Community Z-score")
pDay

# SSH
subsSsh <- subset(all,select=c(grepl("SLA",colnames(all))))
subsSsh<- as.data.frame(t(subsSsh))
subsSsh <- subsSsh %>% summarize_all(mean) %>% as.data.frame()
subsSsh <- data.frame(t(subsSsh))
envB <- all[c(24:ncol(all))]
finSsh <- cbind(subsSsh,envB)
finSsh<- subset(finSsh,!is.na(SLA))
pSsh <- ggplot(finSsh,aes(x=SLA,y=t.subsSsh.))+geom_point(shape=21,fill=colz[4],color="black")+theme_classic()+ggtitle("Sea Surface Height Community (# ASVs = 4)")+xlab("Sea Surface Height (m)")+ylab("Mean Community Z-score")
pSsh

# Fluorescence
subsFlu <- subset(all,select=c(grepl("CTDFLUOR",colnames(all))))
subsFlu<- as.data.frame(t(subsFlu))
subsFlu <- subsFlu %>% summarize_all(mean) %>% as.data.frame()
subsFlu <- data.frame(t(subsFlu))
envB <- all[c(24:ncol(all))]
finFlu <- cbind(subsFlu,envB)
finFlu<- subset(finFlu,!is.na(CTDFLUOR))
pFlu <- ggplot(finFlu,aes(x=CTDFLUOR,y=t.subsFlu.))+geom_point(shape=21,fill=colz[5],color="black")+theme_classic()+ggtitle(expression("Chlorophyll"~italic(a)~"Community (# ASVs = 2)"))+xlab(expression("Chlorophyll"~italic(a)~"("*mu*"g L"^{-1}*")"))+ylab("Mean Community Z-score")
pFlu

# Nitrate
subsNo3 <- subset(all,select=c(grepl("NO3",colnames(all))))
subsNo3<- as.data.frame(t(subsNo3))
subsNo3 <- subsNo3 %>% summarize_all(mean) %>% as.data.frame()
subsNo3 <- data.frame(t(subsNo3))
envB <- all[c(23:ncol(all))]
finNo3 <- cbind(subsNo3,envB)
finNo3<- subset(finNo3,!is.na(NO3))
pNo3 <- ggplot(finNo3,aes(x=NO3,y=t.subsNo3.))+geom_point(shape=21,fill=colz[6],color="black")+theme_classic()+ggtitle("Nitrate (# ASVs = 2)")+xlab(expression("NO"[3]^{"-"}~"("*mu*"M)"))+ylab("Mean Community Z-score")
pNo3

# PA Sar
subsSar <- subset(all,select=c(grepl("sar11_PA",colnames(all))))
subsSar<- as.data.frame(t(subsSar))
subsSar<- subsSar %>% summarize_all(mean) %>% as.data.frame()
subsSar <- data.frame(t(subsSar))
envB <- all[c(23:ncol(all))]
finSar <- cbind(subsSar,envB)
finSar<- subset(finSar,!is.na(sar11_PA))
pSar <- ggplot(finSar,aes(x=sar11_PA,y=t.subsSar.))+geom_point(shape=21,fill=colz[7],color="black")+theme_classic()+ggtitle("Particle-associated SAR11 (# ASVs = 2)")+xlab("Relative Abundance of Particle-associated SAR11")+ylab("Mean Community Z-score")
pSar

# CTDBEAM
subsBeam <- subset(all,select=c(grepl("CTDBEAM",colnames(all))))
subsBeam<- as.data.frame(t(subsBeam))
subsBeam<- subsBeam %>% summarize_all(mean) %>% as.data.frame()
subsBeam <- data.frame(t(subsBeam))
envB <- all[c(23:ncol(all))]
finBeam <- cbind(subsBeam,envB)
finBeam<- subset(finBeam,!is.na(CTDBEAM))
pBeam <- ggplot(finBeam,aes(x=CTDBEAM,y=t.subsBeam.))+geom_point(shape=21,fill=colz[8],color="black")+theme_classic()+ggtitle("Beam attenuation (# ASVs = 1)")+xlab(expression("Beam attenuation (m"^{"-1"}*")"))+ylab("Mean Community Z-score")
pBeam

# CTDOXY
subsOxy <- subset(all,select=c(grepl("CTDOXY",colnames(all))))
subsOxy<- as.data.frame(t(subsOxy))
subsOxy<- subsOxy %>% summarize_all(mean) %>% as.data.frame()
subsOxy <- data.frame(t(subsOxy))
envB <- all[c(23:ncol(all))]
finOxy <- cbind(subsOxy,envB)
finOxy<- subset(finOxy,!is.na(CTDOXY))
pOxy <- ggplot(finOxy,aes(x=CTDOXY,y=t.subsOxy.))+geom_point(shape=21,fill=colz[9],color="black")+theme_classic()+ggtitle("Oxygen (# ASVs = 1)")+xlab(expression("Oxygen (mg L"^{"-1"}*")"))+ylab("Mean Community Z-score")
pOxy

# PA cyano
subsCyanPA <- subset(all,select=c(grepl("cyanos_PA",colnames(all))))
subsCyanPA<- as.data.frame(t(subsCyanPA))
subsCyanPA<- subsCyanPA %>% summarize_all(mean) %>% as.data.frame()
subsCyanPA <- data.frame(t(subsCyanPA))
envB <- all[c(23:ncol(all))]
finCyanPA <- cbind(subsCyanPA,envB)
finCyanPA<- subset(finCyanPA,!is.na(cyanos_PA))
pCyanPA <- ggplot(finCyanPA,aes(x=cyanos_PA,y=t.subsCyanPA.))+geom_point(shape=21,fill=colz[10],color="black")+theme_classic()+ggtitle("Particle-associated cyanobacteria (# ASVs = 1)")+xlab("Relative abundance of particle-associated cyanobacteria")+ylab("Mean Community Z-score")
pCyanPA

# MEI
subsMEI<- subset(all,select=c(grepl("MEI",colnames(all))))
subsMEI$MEI <- NULL
subsMEI<- as.data.frame(t(subsMEI))
subsMEI<- subsMEI %>% summarize_all(mean) %>% as.data.frame()
subsMEI <- data.frame(t(subsMEI))
envB <- all[c(62:ncol(all))]
finMEI <- cbind(subsMEI,envB)
finMEI<- subset(finMEI,!is.na(MEI))
pMEI <- ggplot(finMEI,aes(x=MEI,y=t.subsMEI.))+geom_point(shape=21,fill=colrs[10],color="black")+theme_classic()+ggtitle("MEI")+xlab("MEI")+ylab("Mean Community Z-score")
pMEI

pCyan+pDay+pSsh+pTmp+pFlu+plot_layout(ncol = 3)
ggsave("../try.pdf",width=13,height=7)

# Let's visualize these clusters
interEdge <- read.csv("Surf_DCM_edges.csv",header=TRUE,row.names=1)
interEdge$V1 <- str_remove_all(interEdge$V1,"S_")
interEdge$V2 <- str_remove_all(interEdge$V2,"S_")
interEdge <- as.matrix(interEdge)

g <- graph_from_edgelist(interEdge,directed=FALSE)
names <- data.frame(ASV_ID=c(V(g)$name))
rfPerVar <- subset(rfPerVar,ASV_ID %in% namez)
names <- left_join(names,rfPerVar)
names$col <- ifelse(names$Param=="CTDTMP",colz[1],"grey")
names$col <- ifelse(names$Param=="cyanos",colz[2],names$col)
names$col <- ifelse(names$Param=="DayOfYear",colz[3],names$col)
names$col <- ifelse(names$Param=="SLA",colz[4],names$col)
names$col <- ifelse(names$Param=="CTDFLUOR",colz[5],names$col)
names$col <- ifelse(names$Param=="NO3",colz[6],names$col)
names$col <- ifelse(names$Param=="sar11_PA",colz[7],names$col)
names$col <- ifelse(names$Param=="CTDBEAM",colz[8],names$col)
names$col <- ifelse(names$Param=="CTDOXY",colz[9],names$col)
names$col <- ifelse(names$Param=="cyanos_PA",colz[10],names$col)
names$col <- ifelse(is.na(names$col),"grey",names$col)
unique(names$col)
V(g)$col <- names$col
names$level <- as.numeric(as.factor(names$col))
V(g)$level <- names$level

weight.community <- function(row,membership,weight.within,weight.between){
  if(as.numeric(membership[which(names(membership)==row[1])])==as.numeric(membership[which(names(membership)==row[2])])){
    weight=weight.within
  }else{
    weight=weight.between
  }
  return(weight)
}

membership <- c(V(g)$level)
names(membership) <-c(V(g)$name)
E(g)$weight=apply(get.edgelist(g),1,weight.community,membership,1000,1)
g$layout=layout.fruchterman.reingold(g,weights=E(g)$weight)

pdf("../Plot.pdf")
plot(g,vertex.size=7,vertex.label=NA,edge.width=0.5,vertex.color=V(g)$col,edge.color="grey",main="")
dev.off()

# Make legend
names$x <- 1
names$y <- 2
colorz <- c(colz,"grey")
ggplot(names,aes(x=x,y=y,fill=Param))+geom_point(shape=21)+scale_fill_manual(values=c(colorz),breaks=c("CTDTMP","cyanos","DayOfYear","SLA","CTDFLUOR","NO3","sar11_PA","CTDBEAM","CTDOXY","cyanos_PA",NA),labels=c("Temperature","Relative abundance of free-living cyanobacteria","Day of Year","Sea surface height","Chlorophyll a","Nitrate","Relative abundance of particle-associated SAR11","Beam attenuation","Oxygen","Relative abundance of particle-associated cyanobacteria","< 30% variability explained by environmental variables"))+theme_classic()
ggsave("../legend.pdf")
################################################
### SPOT Network Analysis ###
### Visualize Networks + Calculate Network Statistics ###
### By: Samantha Gleich ###
### Last Updated: 12/28/22 ###

taxCols <- c("#E1746D","#76C3D7","#DE5AB1","#D5E0AF","#DED3DC","#87EB58","#D4DC60","#88E5D3","#88AAE1","#DBA85C","#8B7DDA","#9A8D87","#D99CD1","#B649E3","#7EDD90")
names(taxCols) <- c("Chlorophyte","Ciliate","Cryptophyte","Diatom","Dinoflagellate","Fungi","Haptophyte","MAST","Other Alveolate","Other Archaeplastida","Other Eukaryote","Other Stramenopiles","Rhizaria","Syndiniales","Unknown Eukaryote")



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
t <- left_join(t,taxColsDf)
dfSurf <- data.frame(name=c(V(gSurfOnly)$name))
dfSurf$name <- str_remove_all(dfSurf$name,"S_")
dfSurf <- left_join(dfSurf,t)
V(gSurfOnly)$col <- dfSurf$taxCols

pdf("../../Surface_Network.pdf")
plot(gSurfOnly,vertex.size=4,vertex.label=NA,edge.width=0.5,layout=layout_with_kk(gSurfOnly),vertex.color=V(gSurfOnly)$col,edge.color="grey30",main="Associations only found at surface (5 m)")
dev.off()

dfDCM <- data.frame(name=c(V(gDcmOnly)$name))
dfDCM$name <- str_remove_all(dfDCM$name,"S_")
dfDCM <- left_join(dfDCM,t)
V(gDcmOnly)$col <- dfDCM$taxCols

pdf("../../DCM_Network.pdf")
plot(gDcmOnly,vertex.size=4,vertex.label=NA,edge.width=0.5,layout=layout_with_kk(gDcmOnly),vertex.color=V(gDcmOnly)$col,edge.color="grey30",main="Associations only found at DCM")
dev.off()

pdf("../../Intersect_Network.pdf")
plot(gInterOnly,vertex.size=V(gInterOnly)$siz,vertex.label=NA,edge.width=0.5,layout=layout_with_kk(gInterOnly),vertex.color=V(gInterOnly)$col,edge.color="grey30",main="Associations found at surface and DCM (merged network)")
dev.off()

# Network statistics
length(V(gInterOnly))
length(E(gInterOnly))
V(gInterOnly)$deg <- degree(gInterOnly)
mean(V(gInterOnly)$deg)
mean_distance(gInterOnly)
transitivity(gInterOnly,type="global")

# Size by degree
dfDeg <- data.frame(name=c(V(gInterOnly)$name),deg=c(V(gInterOnly)$deg))
dfDeg <- dfDeg %>% arrange(desc(deg)) %>% as.data.frame()
dfDeg$name <- str_replace_all(dfDeg$name,"S_","")
t <- read.delim("taxonomy_90.tsv",header=TRUE)
t$Confidence <- NULL
colnames(t)[1] <- "name"
dfDeg <- left_join(dfDeg,t)

# Color by tax
taxz <- colsplit(t$Taxon,";",c("Supergroup","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxz$fin <- ifelse(taxz$Class=="Syndiniales","Syndiniales",NA)
taxz$fin <- ifelse(taxz$Class=="Dinophyceae","Dinoflagellate",taxz$fin)
taxz$fin <- ifelse(taxz$Class=="Bacillariophyta","Diatom",taxz$fin)
taxz$fin <- ifelse(grepl("MAST",taxz$Class),"MAST",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Rhizaria","Rhizaria",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Chlorophyta","Chlorophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Cryptophyta","Cryptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Haptophyta","Haptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$Phylum=="Ciliophora","Ciliate",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Stramenopiles" & is.na(taxz$fin),"Other Stramenopiles",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Archaeplastida" & is.na(taxz$fin),"Other Archaeplastids",taxz$fin)
taxz$fin <- ifelse(taxz$Kingdom=="Alveolata" & is.na(taxz$fin),"Other Alveolate",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin),"Other Eukaryote",taxz$fin)
unique(taxz$fin)

dfDeg$tax <- taxz$fin
V(gInterOnly)$tax <- dfDeg$tax
taxColsDf <- data.frame(taxCols)
taxColsDf$tax <- rownames(taxColsDf)
dfDeg <- left_join(dfDeg,taxColsDf)
V(gInterOnly)$col <- dfDeg$taxCols

pdf("../../Degree.pdf")
plot(gInterOnly,vertex.size=V(gInterOnly)$deg*0.5,vertex.label=NA,edge.width=0.5,layout=layout_with_kk(gInterOnly),vertex.color=V(gInterOnly)$col,edge.color="grey30",main="Associations found at surface and DCM (merged network)")
dev.off()

counts <- read.delim("feature-table.tsv",header=TRUE,row.names=1)
counts <- as.data.frame(t(counts))
counts <- counts/rowSums(counts)
counts <- as.data.frame(t(counts))
counts <- subset(counts,rownames(counts) %in% dfDeg$name)
counts$meanz <- rowMeans(counts)
counts$median = apply(counts, 1, median, na.rm=T)
c <- data.frame(name=c(rownames(counts)),meanz=c(counts$meanz))
c <- data.frame(name=c(rownames(counts)),median=c(counts$median))
dfDeg <- left_join(dfDeg,c)
V(gInterOnly)$ab2 <- dfDeg$median

pdf("../../Abundance2.pdf")
plot(gInterOnly,vertex.size=V(gInterOnly)$ab2*1000,vertex.label=NA,edge.width=0.5,layout=layout_with_kk(gInterOnly),vertex.color=V(gInterOnly)$col,edge.color="grey30",main="Associations found at surface and DCM (merged network)")
dev.off()

ggplot(dfDeg,aes(x=deg,y=meanz,fill=tax))+geom_point(shape=21,size=4)+scale_fill_manual(name=c("Taxonomic Group"),values=c(taxCols))+theme_classic(base_size = 14)
ggsave("../../legend.pdf",width=7,height=10)






# Save edgelists 
interEdge <- as.data.frame(interEdge)
surfEdgeOnly <- as.data.frame(surfEdgeOnly)
dcmEdgeOnly <- as.data.frame(dcmEdgeOnly)
write.csv(interEdge,"Surf_DCM_edges.csv")
write.csv(surfEdgeOnly,"Surf_edges.csv")
write.csv(dcmEdgeOnly,"DCM_edges.csv")


df <- read.csv(file.choose(),header=TRUE)
################################################



