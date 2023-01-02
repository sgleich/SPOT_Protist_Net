### SPOT Network Analysis ###
### Run Random Forest Analysis ###
### By: Samantha Gleich ###
### Last Updated: 1/2/23 ###

# Load libraries
library(randomForest)
library(missForest)
library(tidyverse)
library(reshape2)
library(caret)
library(yardstick)

# Set seed for reproducibility
set.seed(100)

# Load in ASV table and manifest
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

# CLR transform data (because they're relative abundances)
dfCLR <- as.data.frame(clr(dfNew))

# Let's only keep those ASVs that we included in our networks
netOut <- read.csv("Glasso_5m_SPOT_Dec28.csv",header=TRUE,row.names=1) 
namez <- colnames(netOut)
namez <- str_remove(namez,"S_")
dfCLR <- subset(dfCLR,select=c(namez))
namez <- colnames(dfCLR)
s <- c(paste("V",1:389,sep="_"))
colnames(dfCLR) <- s

# Load environmental data 
env <- read.csv("../SPOT_Env_Update_2022.csv",header=TRUE,row.names=1)

# Match up environmental data to order of ASV table samples (not necessarily ordered by time)
colz <- colsplit(rownames(dfCLR),"_",c("spot","Cruise","Month","Day","Year","Depth"))
dfCLR$Depth <- colz$Depth
dfCLR$Cruise <- colz$Cruise
all <- left_join(dfCLR,env)
dfCLR$Depth <- NULL
dfCLR$Cruise <- NULL

# Now we can select just the environmental variables we care about
envNew <- all[c(395,397:411,413:418)]
envNew$CTDPRS <- NULL
envNew$CTDCOND <- NULL
envNew <- missForest(envNew)
envNew <- envNew$ximp

# Rename CLR transformed dataframe
total <- dfCLR
total <- mutate_all(total, function(x) as.numeric(as.character(x)))

# Set up training and test dataset for ASV-ASV predictions
finalTest <- total[181:242,] # 25% of data used for test
finalTrain <- total[1:180,] # 75% of data used for train
horizon=12

out <- NULL
out2 <- NULL
# Set up loop 
for (i in 1:ncol(total)){
  # Run RF on ASV i with other ASVs as predictors
  y_train <- finalTrain[, i] # Target ASV
  X_train <- finalTrain[, -i] # Everything but target ASV
  y_test <- finalTest[,i] # Target ASV test
  X_test <- finalTest[,-i] # Everything but target ASV test
  
  tr_control <- caret::trainControl(
    method = 'timeslice',
    initialWindow = nrow(X_train) - horizon,
    fixedWindow = TRUE, savePredictions = "all")
  
  tune_grid <- expand.grid(
    mtry = c(2,5,10,15,19))
  
  res_asv <- caret::train(
    data.frame(X_train),
    y_train,
    method = 'rf',
    trControl = tr_control,
    tuneGrid = tune_grid,
    ntree=100)
  
  y_hat <- predict(res_asv, X_test)
  rf_scored_asv <- as_tibble(cbind(y_test, y_hat))
  RMSE_asv <- yardstick::rmse(rf_scored_asv, truth=y_test,estimate=y_hat)
  RMSE_asv[1,4] <- colnames(total)[i]
  RMSE_asv[1,5] <- "ASVs"
  x_asv <- 100 * (1-sum((res_asv$finalModel$y-res_asv$finalModel$pred)^2) /sum((res_asv$finalModel$y-mean(res_asv$finalModel$y))^2))
  RMSE_asv[1,6] <- x_asv
  RMSE_asv[1,7] <- namez[i]
  colnames(RMSE_asv) <- c("metric","s","estimate","Y","Predictors","PercentVar","ASV_ID")
  
  # Run RF on ASV i with environmental parameters as predictors
  y_train <- finalTrain[, i] 
  X_train <- envNew[1:180,] 
  y_test <- finalTest[,i]
  X_test <- envNew[181:242,]
  
  tr_control2 <- caret::trainControl(
    method = 'timeslice',
    initialWindow = nrow(X_train) - horizon,
    fixedWindow = TRUE, savePredictions = "all")
  
  tune_grid2 <- expand.grid(
    mtry = c(2,5,10,15,19))
  
  
  res_env <- caret::train(
    data.frame(X_train),
    y_train,
    method = 'rf',
    trControl = tr_control2,
    tuneGrid = tune_grid2,
    ntree=100)
  
  y_hat2 <- predict(res_env, X_test)
  rf_scored_env <- as_tibble(cbind(y_test, y_hat2))
  RMSE_env <- yardstick::rmse(rf_scored_env, truth=y_test,estimate=y_hat2)
  RMSE_env[1,4] <- colnames(total)[i] 
  RMSE_env[1,5] <- "Environment"
  x_env <- 100 * (1-sum((res_env$finalModel$y-res_env$finalModel$pred)^2) /sum((res_env$finalModel$y-mean(res_env$finalModel$y))^2))
  RMSE_env[1,6] <- x_env
  RMSE_env[1,7] <- namez[i]
  colnames(RMSE_env) <- c("metric","s","estimate","Y","Predictors","PercentVar","ASV_ID")
  varz <- res_env$finalModel$importance
  varz <- as.data.frame(varz)
  varz$Param <- rownames(varz)
  varz <- varz %>% arrange(desc(IncNodePurity)) %>% as.data.frame()
  varz$name <- namez[i]
  varz <- varz[1,]
  out2 <- rbind(out2,varz)
  
  # Run RF on ASV i with ASV AND environmental parameters as predictors
  allPred <- cbind(total,envNew)
  y_train <- finalTrain[, i] 
  X_train <- allPred[1:180,-i] 
  y_test <- finalTest[,i]
  X_test <- allPred[181:242,-i]
  
  tr_control3 <- caret::trainControl(
    method = 'timeslice',
    initialWindow = nrow(X_train) - horizon,
    fixedWindow = TRUE, savePredictions = "all")
  
  tune_grid3 <- expand.grid(
    mtry = c(2,5,10,15,19))
  
  
  res_all <- caret::train(
    data.frame(X_train),
    y_train,
    method = 'rf',
    trControl = tr_control3,
    tuneGrid = tune_grid3,
    ntree=100)
  
  y_hat3 <- predict(res_all, X_test)
  rf_scored_all <- as_tibble(cbind(y_test, y_hat3))
  RMSE_all <- yardstick::rmse(rf_scored_all, truth=y_test,estimate=y_hat3)
  RMSE_all[1,4] <- colnames(total)[i] 
  RMSE_all[1,5] <- "ASV+ Environment"
  x_all <- 100 * (1-sum((res_all$finalModel$y-res_all$finalModel$pred)^2) /sum((res_all$finalModel$y-mean(res_all$finalModel$y))^2))
  RMSE_all[1,6] <- x_all
  RMSE_all[1,7] <- namez[i]
  colnames(RMSE_all) <- c("metric","s","estimate","Y","Predictors","PercentVar","ASV_ID")
  
  out <- rbind(out,RMSE_asv,RMSE_env,RMSE_all)}

write.csv(out,"randomForest_Results_Dec2022b.csv")
write.csv(out2,"randomForest_Results_Env.csv")
