horizon=12

y_train <- finalTrain[, 2] # the target
X_train <- finalTrain[, -2] # everything but the target
y_test <- finalTest[,2]
X_test <- finalTest[,-2]

tr_control <- caret::trainControl(
  method = 'timeslice',
  initialWindow = nrow(X_train) - horizon,
  fixedWindow = TRUE, savePredictions = "all")

tune_grid <- expand.grid(
  mtry = c(20,80,120,180,220,280,320,380,398))


holdout_result <- caret::train(
  data.frame(X_train),
  y_train,
  method = 'rf',
  trControl = tr_control,
  tuneGrid = tune_grid, ntree=100
)

holdout_result$results
holdout_result$bestTune

y_hat <- predict(holdout_result, X_test)
test.rf_scored <- as_tibble(cbind(y_test, y_hat))
RMSE_rf_TEST <- yardstick::rmse(test.rf_scored, truth=y_test, estimate=y_hat)
