# Figure 7 and Figure SX

# Set working directory
# setwd("/Users/samanthagleich/Desktop/export_dir_feb2024")

# Libraries
library(reshape2)
library(tidyverse)
library(ggplot2)
library(compositions)

# Load data
counts <- read.delim("feature-table.tsv",header=TRUE,skip=1)
manifest <- read.delim("manifest.txt",header=FALSE)
manifest_cols <- colsplit(manifest$V1,",",c("name","sample","direction"))
manifest_cols <- manifest_cols %>% distinct(name,.keep_all = TRUE) %>% as.data.frame()
rownames(counts) <- counts$X.OTU.ID
counts$X.OTU.ID <- NULL
d <- data.frame(name=colnames(counts))
d <- left_join(d,manifest_cols)
sample_cols <- colsplit(d$sample,"_",c("name","other"))
colnames(counts) <- sample_cols$name
counts$`SPOTRe-DNA-115-02-16-12-5m` <- NULL
counts$`SPOTRe-DNA-115-02-16-12-DCM` <- NULL

# CLR transform data
counts <- as.data.frame(t(counts))
counts <- compositions::clr(counts)
counts <- as.data.frame(counts)

# Load in network outputs
networkSurf <- read.csv("Surf_SPOT_Feb12_2024.csv",header=TRUE)
networkSurf <- melt(networkSurf,id.vars="X")
networkSurf <- subset(networkSurf,value==1)
networkSurf$X <- str_remove_all(networkSurf$X,"S_")
networkSurf$variable <- str_remove_all(networkSurf$variable,"S_")

# Remove duplicate rows (ASV1-ASV2 == ASV2-ASV1)
networkSurf2 <- networkSurf[c(1:2)]
networkSurf2 <- data.frame(t(apply(networkSurf2,1,sort)))
networkSurf2 <- networkSurf2[!duplicated(networkSurf2),]

# Select only surface samples
countsSurf <- subset(counts,grepl("5m",rownames(counts)))
countsSurf <- as.data.frame(t(countsSurf))
countsSurf <- subset(countsSurf,rownames(countsSurf) %in% networkSurf2$X1 | rownames(countsSurf) %in% networkSurf2$X2)

# Calculate spearman correlation
countsSurf <- as.data.frame(t(countsSurf))
corSurf <- cor(countsSurf,method="spearman")
corSurf <- as.data.frame(corSurf)
corSurf$x <- rownames(corSurf)

# Find interactions with strong correlation
corSurfMelt <- melt(corSurf,id.vars="x")
corSurfMelt <- subset(corSurfMelt,x!=variable)
corSurfMelt$abs <- abs(corSurfMelt$value)
corSurfMelt <- corSurfMelt %>% arrange(desc(abs))
corSurfMelt$variable <- as.character(corSurfMelt$variable)

# Get taxonomic info
tax <- read.delim("taxonomy.tsv",header=TRUE)
tax$Confidence <- NULL
colnames(tax) <- c("x","xTax")
corSurfMelt <- left_join(corSurfMelt,tax)
colnames(tax) <- c("variable","variableTax")
corSurfMelt <- left_join(corSurfMelt,tax)

# Remove interactions that are likely due to intra-species diversity (same species; multiple ASVs)
t <- subset(corSurfMelt,xTax!=variableTax)
t2 <- t[c(1:2)]
# Remove duplicates (ASV1 - ASV2 == ASV2 - ASV1)
t2 <- data.frame(t(apply(t2,1,sort)))
t2 <- t2[!duplicated(t2),]
colnames(t)[1:2] <- c("X1","X2")
t2 <- left_join(t2,t)
colnames(t2) <- c("X1","X2","Value","Abs","X1Tax","X2Tax")


t2$YN <- "No"
for(i in 1:200){
  tmp <- t2[i,]
  v1 <- tmp$X1
  v2 <- tmp$X2
  d <- subset(countsSurf,select=c(v1,v2))
  colnames(d) <- c("t1","t2")
  
  # Find interactions where both ASVs are non-zero > 20 times
  if(sum(d$t1==0) < 20 & sum(d$t2==0) < 20){
    tmp$YN <- "Yes"
    t2[i,] <- tmp
  }
}

t2 <- subset(t2,YN=="Yes")
t2 <- subset(t2,X1!=X2)

surfRows <- c(2,3,5,7,9,10,12,13)
#dcmRows <- c(1,2,3,4)

t2[2,]

plotFxn1 <- function(num){
  v1 <- t2$X1[num]
  v2 <- t2$X2[num]
  d <- subset(countsSurf,select=c(v1,v2))
  colnames(d) <- c("t1","t2")
  
  scc_val <- round(t2$Value[num],2)
  
 
  ggplot(d,aes(t1,t2))+geom_point()+theme_bw()+xlab("ASV")+ylab("ASV")+ggtitle(paste("SCC =",scc_val,sep=" "))
}

#surfPlots <- p1+p2+p3+p4+p5+p6+p7+p8+plot_annotation(tag_levels="a",title="Surface")+plot_layout(nrow=2)
#ggsave("../surfPlots.pdf",width=10,height=3)


plotFxn2 <- function(num,dep,asv1,asv2){
  v1 <- t2$X1[num]
  v2 <- t2$X2[num]
  d <- subset(countsSurf,select=c(v1,v2))
  env <- read.csv("SPOT_Env_Feb2024_Final.csv",header=TRUE)
  c <- colsplit(rownames(d),"-",c("SPOT","DNA","Cruise","Other"))
  c <- left_join(c,env)
  c <- subset(c,Depth==dep)
  c$Date <- as.Date.character(c$Date,format="%m/%d/%y")
  d$Date <- c$Date
  colnames(d) <- c("ASV1","ASV2","Date")
  dMelt <- melt(d,id.vars="Date")
  ggplot(dMelt,aes(x=Date,y=value,color=variable))+geom_point()+geom_line()+theme_classic()+scale_x_date(breaks="1 year",date_labels="%Y")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("Date of Sample Collection")+ylab("CLR-transformed Abundance")+scale_color_manual(name="ASVs",values=c("black","grey"),labels=c("ASV1","ASV2"))
}
