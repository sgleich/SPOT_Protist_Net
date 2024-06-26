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
counts <- as.data.frame(t(counts))
counts <- compositions::clr(counts)
counts <- as.data.frame(counts)

networkSurf <- read.csv("Surf_SPOT_Feb12_2024.csv",header=TRUE)
networkSurf <- melt(networkSurf,id.vars="X")
networkSurf <- subset(networkSurf,value==1)
networkSurf$X <- str_remove_all(networkSurf$X,"S_")
networkSurf$variable <- str_remove_all(networkSurf$variable,"S_")

networkSurf2 <- networkSurf[c(1:2)]
networkSurf2 <- data.frame(t(apply(networkSurf2,1,sort)))
networkSurf2 <- networkSurf2[!duplicated(networkSurf2),]


countsSurf <- subset(counts,grepl("5m",rownames(counts)))
countsSurf <- as.data.frame(t(countsSurf))
countsSurf <- subset(countsSurf,rownames(countsSurf) %in% networkSurf2$X1 | rownames(countsSurf) %in% networkSurf2$X2)

countsSurf <- as.data.frame(t(countsSurf))
corSurf <- cor(countsSurf,method="spearman")
corSurf <- as.data.frame(corSurf)
corSurf$x <- rownames(corSurf)
corSurfMelt <- melt(corSurf,id.vars="x")
corSurfMelt <- subset(corSurfMelt,x!=variable)
corSurfMelt$abs <- abs(corSurfMelt$value)
corSurfMelt <- corSurfMelt %>% arrange(desc(abs))
corSurfMelt$variable <- as.character(corSurfMelt$variable)
tax <- read.delim("taxonomy.tsv",header=TRUE)
tax$Confidence <- NULL
colnames(tax) <- c("x","xTax")
corSurfMelt <- left_join(corSurfMelt,tax)
colnames(tax) <- c("variable","variableTax")
corSurfMelt <- left_join(corSurfMelt,tax)


t <- subset(corSurfMelt,xTax!=variableTax)
t <- t[1:200,]
t2 <- t[c(1:2)]
t2 <- data.frame(t(apply(t2,1,sort)))
t2 <- t2[!duplicated(t2),]
colnames(t)[1:2] <- c("X1","X2")
t2 <- left_join(t2,t)
colnames(t2) <- c("X1","X2","Value","Abs","X1Tax","X2Tax")

t2$YN <- "No"
for(i in 1:nrow(t2)){
  tmp <- t2[i,]
  v1 <- tmp$X1
  v2 <- tmp$X2
  d <- subset(countsSurf,select=c(v1,v2))
  colnames(d) <- c("t1","t2")
  if(sum(d$t1==0) < 20 & sum(d$t2==0) < 20){
    tmp$YN <- "Yes"
    t2[i,] <- tmp
  }
}
  
t2 <- subset(t2,YN=="Yes")
t2 <- subset(t2,X1!=X2)


# 2,3,5,7 surf
# 1,2,3,4 dcm
t2[2,]  
v1 <- t2$X1[4]
v2 <- t2$X2[4]
d <- subset(countsSurf,select=c(v1,v2))
colnames(d) <- c("t1","t2")

p4 <- ggplot(d,aes(t1,t2))+geom_point()+theme_bw()+xlab(expression(italic("Pelagomonas calceolata")~"ASV"))+ylab(expression(italic("Phaeocystis cordata")~"ASV"))+ggtitle("SCC = 0.70")

p4
p3
p2
p1

dcmPlots <- p1+p2+p3+p4+plot_annotation(tag_levels="e",title="DCM")+plot_layout(nrow=1)
ggsave("../dcmPlots.pdf",width=10,height=3)
