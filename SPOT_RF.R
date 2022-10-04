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
env <- read.csv("environmental.csv",header=TRUE,row.names=1)
envNewImp <- env

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
  
  # Run RF on ASV i with ASV AND environmental parameters as predictors
  allPred <- cbind(total,envNewImp)
  y_train <- finalTrain[, i] 
  X_train <- allPred[1:180,-i] 
  y_test <- finalTest[,i]
  X_test <- allPred[181:242,-i]
  
  tr_control3 <- caret::trainControl(
    method = 'timeslice',
    initialWindow = nrow(X_train) - horizon,
    fixedWindow = TRUE, savePredictions = "all")
  
  tune_grid3 <- expand.grid(
    mtry = c(20,80,120,180,220,280,320,350,380,398,409))
  
  
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
  colnames(RMSE_all) <- c("metric","s","estimate","Y","Predictors")
  
  out <- rbind(out,RMSE_asv,RMSE_env,RMSE_all)}

write.csv(out,"randomForest_Results_Both.csv")
