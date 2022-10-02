### SPOT Random Forest ###
# Load libraries
library(randomForest)
library(missForest)
library(tidyverse)
library(reshape2)
library(caret)
library(yardstick)

# Load in data (data are filtered and CLR transformed as was done for the network analysis)
surf <- read.csv("SPOT_5m_RF.csv",header=TRUE,row.names=1)
dcm <- read.csv("SPOT_DCM_RF.csv",header=TRUE,row.names=1)
surf <- as.data.frame(t(surf)) 
dcm <- as.data.frame(t(dcm)) 
# Combine surface and DCM data
total <- rbind(surf,dcm)

# Load environmental data 
env <- read.csv("Nut_env.csv",header=TRUE)
env2 <- read.csv("CTD_env.csv",header=TRUE)

# Clean up environmental data
colz <- colsplit(rownames(total),"\\.",c("Spot","dna","num","month","day","year","depth"))
l <- c(3,4,5,6,7,8,9)
colz$year <- ifelse(colz$year %in% l, paste(200,colz$year,sep=""), paste(20,colz$year,sep=""))
colz <- data.frame(Month=colz$month,Day=colz$day,Year=colz$year,Depth=colz$depth)
colz$Year <- as.integer(colz$Year)
colz$Day <- as.integer(colz$Day)
colz$Month <- as.integer(colz$Month)
colz$Depth <- ifelse(grepl("5m",colz$Depth),"5m","DCM")
env <- subset(env,Depth=="5m"|Depth=="DCM")
colz2 <- colsplit(env2$date,"/",c("month","day","year"))
env2$Year <- colz2$year
env2$Day <- colz2$day
l <- c(3,4,5,6,7,8,9)
env2$Year <- ifelse(env2$Year %in% l, paste(200,env2$Year,sep=""), paste(20,env2$Year,sep=""))
env2$Year <- as.numeric(env2$Year)
env2 <- env2[c(4,8:12,14,19:20)]

# Let's order the ASV data in time order (block CV)
total <- cbind(total,colz)
total <- total %>% arrange(Year,Month,Day)%>%as.data.frame()

# Now grab environmental data that aligns with ASV data
colz <- total[,399:402]
total <- total[-c(399:402)]
envNew <- left_join(env,env2)
envNew <- left_join(colz,envNew)
envNew2 <- envNew %>% distinct(Month,Day,Year,Depth, .keep_all = TRUE) %>% as.data.frame()
envNew2 <- envNew2[c(6:11,14,16:20)]

# Impute missing environmental data points
envNew2 <- mutate_all(envNew2, function(x) as.numeric(as.character(x)))
envNewImp <- missForest(envNew2)
envNewImp <- envNewImp$ximp

# Set up training and test dataset for ASV-ASV predictions
finalTest <- total[181:242,] # 25% of data used for test
finalTrain <- total[1:180,] # 75% of data used for train
horizon=12

out <- NULL
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
  mtry = c(20,80,120,180,220,280,320,350,380,398))
  
  
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
  colnames(RMSE_asv) <- c("metric","s","estimate","Y","Predictors")
  
  # Run RF on ASV i with environmental parameters as predictors
  y_train <- finalTrain[, i] 
  X_train <- envNewImp[1:180,] 
  y_test <- finalTest[,i]
  X_test <- envNewImp[181:242,]

  tr_control2 <- caret::trainControl(
  method = 'timeslice',
  initialWindow = nrow(X_train) - horizon,
  fixedWindow = TRUE, savePredictions = "all")
  
  tune_grid2 <- expand.grid(
  mtry = c(20,80,120,180,220,280,320,350,380,398))


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
  colnames(RMSE_env) <- c("metric","s","estimate","Y","Predictors")
  
  out <- rbind(out,RMSE_asv,RMSE_env)}

write.csv(out,"randomForest_Results.csv")

