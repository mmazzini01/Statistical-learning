## library
library(readr)
library(rmarkdown)
library(corrplot)
library(glmnet)
library(knitr)
library(kableExtra)
library(DescTools)
library(pROC)
library(class)

#Data Pre-Processing
data <- read_csv('heart_data.csv', show_col_types = FALSE)
summary(data)
str(data)

data$Sex <- as.factor(data$Sex)
data$ChestPainType <- as.factor(data$ChestPainType)
data$RestingECG <- as.factor(data$RestingECG)
data$ExerciseAngina <- as.factor(data$ExerciseAngina)
data$ST_Slope <- as.factor(data$ST_Slope)
data$HeartDisease <- as.factor(data$HeartDisease)
data$FastingBS <- as.factor(data$FastingBS)

# Filling missing values
data$Cholesterol[data$Cholesterol == 0] <- NA
data$Cholesterol[is.na(data$Cholesterol)]<- median(data$Cholesterol, na.rm = TRUE)
data$RestingBP[data$RestingBP == 0] <- NA
data$RestingBP[is.na(data$RestingBP)]<- median(data$RestingBP, na.rm = TRUE)
summary(data)

#train test splitting
#we want to respect the proportion of response variable in both train and test set
#we use a stratified sampling
# Proporzione del training set
set.seed(123)
train_proportion <- 0.75
train_set <- data.frame()
test_set <- data.frame()
for (level in unique(data$HeartDisease)) {
  subset_data <- data[data$HeartDisease == level, ]
  train_size <- round(nrow(subset_data) * train_proportion)
  train_indices <- sample(1:nrow(subset_data), train_size)
  train_subset <- subset_data[train_indices, ]
  test_subset <- subset_data[-train_indices, ]
  train_set <- rbind(train_set, train_subset)
  test_set <- rbind(test_set, test_subset)
}

#check proportion
train_prop <- round(prop.table(table(train_set$HeartDisease)),3)
test_prop <- round(prop.table(table(test_set$HeartDisease)),3)
data_prop<- round(prop.table(table(data$HeartDisease)),3)
cat('data: ',data_prop,'train: ',train_prop,'test: ',test_prop)

# X and y

X_train <- model.matrix(HeartDisease~., train_set)[,-1]
y_train <- as.numeric(as.character(train_set$HeartDisease))

X_test <- model.matrix(HeartDisease~., test_set)[,-1]
y_test <- as.numeric(as.character(test_set$HeartDisease))


#Ridge regression
ridge_cv <- cv.glmnet(X_train, y_train, alpha = 0, family = "binomial", type.measure = "deviance", nfolds = 10)
plot(ridge_cv)
lambda = ridge_cv$lambda.min
cat("The value for the minimum lambda is ", lambda)
ridge_coef <- coef(ridge_cv, s = "lambda.min")
print(ridge_coef)

pred_ridge_prob<- predict(ridge_cv, X_test, type = "response", s = lambda)
pred_ridge <- ifelse(pred_ridge_prob > 0.5, 1, 0)

# Create confusion matrix
confusion_matrix <- table(y_test, pred_ridge)
# Format confusion matrix as a data frame with 2 columns
confusion_df <- as.data.frame(matrix(0, nrow = 2, ncol = 2))
colnames(confusion_df) <- c("True Non-Disiased", "True Disiased")
rownames(confusion_df) <- c("Pred. Non-Disiased", "Pred. Disiased")
# Fill the confusion matrix with values
confusion_df[, 1:2] <- t(confusion_matrix)
# Add a total row to the data frame
confusion_df["Total",] <- colSums(confusion_df)
# Add a total column to the data frame
confusion_df <- cbind(confusion_df, Total = rowSums(confusion_df))

confusion_df

acc_ridge <- acc_knn <- round((confusion_matrix[1,1]+confusion_matrix[2,2])/sum(confusion_matrix),2)
cat("Accuracy:", acc_ridge, "\n")

prec_ridge <- round(confusion_matrix[2,2]/(confusion_matrix[2,2]+confusion_matrix[1,2]),2)
cat("Precision:", prec_ridge, "\n")

rec_ridge <- round(confusion_matrix[2,2]/(confusion_matrix[2,2]+confusion_matrix[2,1]),2)
cat("Recall:", rec_ridge, "\n")

spec_ridge <- round(confusion_matrix[1,1]/(confusion_matrix[1,1]+confusion_matrix[1,2]),2)
cat("Specificity:", spec_ridge, "\n")

type_ridge <- round(confusion_matrix[1,2]/(confusion_matrix[1,2]+confusion_matrix[1,1]),2)
cat("Type 1 error:", type_ridge, "\n")

roc_out <- roc(y_test, as.numeric(pred_ridge_prob))
plot(roc_out, print.auc = TRUE, xlab="False positive rate(Type 1 error)", ylab="True positive rate(Recall)", legacy.axes = TRUE)

#Accuracy: 0.91
#Precision: 0.91
#Recall: 0.94
#Specificity: 0.88
#Type 1 error: 0.12
#AUC: 0.954

#Lasso regression
lasso_cv <- cv.glmnet(X_train, y_train, alpha = 1, family = "binomial", type.measure = "deviance", nfolds = 10)
plot(lasso_cv)
lambda = lasso_cv$lambda.min
cat("The value for the minimum lambda is ", lambda)
lasso_coef <- coef(lasso_cv, s = "lambda.min")
print(lasso_coef)

pred_lasso_prob<- predict(lasso_cv, X_test, type = "response", s = lambda)
pred_lasso <- ifelse(pred_lasso_prob > 0.5, 1, 0)

# Create confusion matrix
confusion_matrix <- table(y_test, pred_lasso)
# Format confusion matrix as a data frame with 2 columns
confusion_df <- as.data.frame(matrix(0, nrow = 2, ncol = 2))
colnames(confusion_df) <- c("True Non-Disiased", "True Disiased")
rownames(confusion_df) <- c("Pred. Non-Disiased", "Pred. Disiased")
# Fill the confusion matrix with values
confusion_df[, 1:2] <- t(confusion_matrix)
# Add a total row to the data frame
confusion_df["Total",] <- colSums(confusion_df)
# Add a total column to the data frame
confusion_df <- cbind(confusion_df, Total = rowSums(confusion_df))

confusion_df

acc_lasso <- round((confusion_matrix[1,1]+confusion_matrix[2,2])/sum(confusion_matrix),2)
cat("Accuracy:", acc_lasso, "\n")

prec_lasso <- round(confusion_matrix[2,2]/(confusion_matrix[2,2]+confusion_matrix[1,2]),2)
cat("Precision:", prec_lasso, "\n")

rec_lasso <- round(confusion_matrix[2,2]/(confusion_matrix[2,2]+confusion_matrix[2,1]),2)
cat("Recall:", rec_lasso, "\n")

spec_lasso <- round(confusion_matrix[1,1]/(confusion_matrix[1,1]+confusion_matrix[1,2]),2)
cat("Specificity:", spec_lasso, "\n")

type_lasso <- round(confusion_matrix[1,2]/(confusion_matrix[1,2]+confusion_matrix[1,1]),2)
cat("Type 1 error:", type_lasso, "\n")

roc_out <- roc(y_test, as.numeric(pred_lasso_prob))
plot(roc_out, print.auc = TRUE, xlab="False positive rate(Type 1 error)", ylab="True positive rate(Recall)", legacy.axes = TRUE)

#Accuracy: 0.91
#Precision: 0.91
#Recall: 0.93
#Specificity: 0.89
#Type 1 error: 0.11
#AUC: 0.952

# KNN
# Scaling dataset
X_train[,"Age"] <- scale(X_train[,"Age"])
X_train[, "RestingBP"] <- scale(X_train[, "RestingBP"])
X_train[,"Cholesterol"] <- scale(X_train[,"Cholesterol"])
X_train[,"MaxHR"] <- scale(X_train[,"MaxHR"])
X_train[,"Oldpeak"] <- scale(X_train[,"Oldpeak"])

X_test[, "Age"] <- scale(X_test[, "Age"])
X_test[, "RestingBP"] <- scale(X_test[, "RestingBP"])
X_test[, "Cholesterol"] <- scale(X_test[, "Cholesterol"])
X_test[, "MaxHR"] <- scale(X_test[, "MaxHR"])
X_test[, "Oldpeak"] <- scale(X_test[, "Oldpeak"])

k <- 5
pred_knn <- knn(train = X_train, test = X_test, cl = y_train, k = k)

# Create confusion matrix
confusion_matrix <- table(y_test, pred_knn)
# Format confusion matrix as a data frame with 2 columns
confusion_df <- as.data.frame(matrix(0, nrow = 2, ncol = 2))
colnames(confusion_df) <- c("True Non-Disiased", "True Disiased")
rownames(confusion_df) <- c("Pred. Non-Disiased", "Pred. Disiased")
# Fill the confusion matrix with values
confusion_df[, 1:2] <- t(confusion_matrix)
# Add a total row to the data frame
confusion_df["Total",] <- colSums(confusion_df)
# Add a total column to the data frame
confusion_df <- cbind(confusion_df, Total = rowSums(confusion_df))

confusion_df

acc_knn <- round((confusion_matrix[1,1]+confusion_matrix[2,2])/sum(confusion_matrix),2)
cat("Accuracy:", acc_knn, "\n")

prec_knn <- round(confusion_matrix[2,2]/(confusion_matrix[2,2]+confusion_matrix[1,2]),2)
cat("Precision:", prec_knn, "\n")

rec_knn <- round(confusion_matrix[2,2]/(confusion_matrix[2,2]+confusion_matrix[2,1]),2)
cat("Recall:", rec_knn, "\n")

spec_knn <- round(confusion_matrix[1,1]/(confusion_matrix[1,1]+confusion_matrix[1,2]),2)
cat("Specificity:", spec_knn, "\n")

type_knn <- round(confusion_matrix[1,2]/(confusion_matrix[1,2]+confusion_matrix[1,1]),2)
cat("Type 1 error:", type_knn, "\n")

#Accuracy: 0.89
#Precision: 0.89
#Recall: 0.91
#Specificity: 0.86
#Type 1 error: 0.14

# Conclusion
# The three models have similar performances, with the Ridge regression model we achive a little bit better results in terms of recall.
