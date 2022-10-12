intEdgeSub <- subset(intEdge2,C1=="Syndiniales"|C2=="Syndiniales")
syn <- subset(pidaX,Taxonomic.level.3..org1=="Syndiniales" | Taxonomic.level.3..org2=="Syndiniales")
all <- data.frame(c=c(syn$Genus.org1,syn$Genus.org2))
all <- unique(all$c)
  
sOut <- NULL
for (item in 1:length(all)){
  s <- subset(intEdgeSub,grepl(paste(all[i]),intEdgeSub$G2)|grepl(paste(all[i]),intEdgeSub$G1))
  sOut <- rbind(sOut,s)
}

sub2 <-  intEdgeSub %>% filter(C1=="Syndiniales"& C2=="Syndiniales") 

total <- rbind(sOut,sub2)

write.csv(total,"Parasitism_Output.csv")
