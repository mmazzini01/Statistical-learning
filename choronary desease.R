## library
library(readr)
library(rmarkdown)
library(corrplot)
library(MASS)
library(car)
data <- read.csv('heart_data.csv', header = T, sep = ',')
#####Pre-processing
#change categorical into factor
data$Sex <- as.factor(data$Sex)
data$ChestPainType <- as.factor(data$ChestPainType)
data$RestingECG <- as.factor(data$RestingECG)
data$ExerciseAngina <- as.factor(data$ExerciseAngina)
data$ST_Slope <- as.factor(data$ST_Slope)
data$FastingBS <- as.factor(data$FastingBS)
data$HeartDisease <- as.factor(data$HeartDisease)
attach(data)
any(is.na(data))
summary(data)
#from the summary we can see that there in non missing value but in the variable cholesterol and restingBp there values 
#that are equal to zero (172 and 1 respectively), this could be missing value or errors in the recording since in the first
#case the sierum cholesterol can't be qual to zero, especially for 172 patients. in the second case blood pressure in 
#human alive cant be zero, and the patient is alivr since heartj records are different from zero. we procced by substituing thi
#zero values with the median that is less sensitive than the mean from extrem observations
data[RestingBP == 0, c("MaxHR", "RestingECG",'HeartDisease')]
data[Cholesterol == 0, c("MaxHR", "RestingECG",'HeartDisease')]
detach(data)
data$Cholesterol[data$Cholesterol==0] <- NA
data$RestingBP[data$RestingBP==0] <- NA
data$Cholesterol[is.na(data$Cholesterol)] <- median(data$Cholesterol, na.rm = TRUE)
data$RestingBP[is.na(data$RestingBP)] <- median(data$RestingBP, na.rm = TRUE)
any(is.na(data))
summary(data)# data are ready to be processed
#stand  log
numeric_vars <- data[, sapply(data, is.numeric)]
cat_vars <- data[, sapply(data, is.factor)]
data_num_scaled <- scale(numeric_vars)
data_num_scaled_df <- as.data.frame(data_num_scaled)
data <- cbind(data_num_scaled_df, cat_vars)

#####Explorative Data Analysis
##Univariate Analysis
plot_histogram <- function(data, variable) {
  var_data <- data[[variable]]
  hist(var_data, main = paste("Histogram of", variable),
       xlab = variable, col = "lightblue")
}
plot_boxplot <- function(data, variable) {
  var_data <- data[[variable]]
  boxplot(var_data, main = paste("Boxplot of", variable),
       xlab = variable, col = "lightblue")
}
#age
par(mfrow = c(1, 2))
plot_histogram(data,'Age')
plot_boxplot(data,'Age')
#restingBP
par(mfrow = c(1, 2))
plot_histogram(data,'RestingBP')
plot_boxplot(data,'RestingBP')
#Cholesterol
par(mfrow = c(1, 2))
plot_histogram(data,'Cholesterol')
plot_boxplot(data,'Cholesterol')
#MaxHR
par(mfrow = c(1, 2))
plot_histogram(data,'MaxHR')
plot_boxplot(data,'MaxHR')
#oldpeack
par(mfrow = c(1, 2))
plot_histogram(data,'Oldpeak')
plot_boxplot(data,'Oldpeak')
##numerical variables
par(mfrow = c(2, 3))
for (i in 1:ncol(numeric_vars)) {
  hist(numeric_vars[, i], main = colnames(numeric_vars)[i],
       xlab = colnames(numeric_vars)[i],col = 'lightblue', freq = FALSE)
  dens <- density(numeric_vars[, i], na.rm=TRUE, adjust=1.25)
  lines(dens, col = "black", lwd = 1)
}

par(mfrow = c(2, 3))
for (i in 1:ncol(numeric_vars)) {
  boxplot(numeric_vars[, i], main = colnames(numeric_vars)[i],
       xlab = colnames(numeric_vars)[i],col = 'lightblue', freq = FALSE)
}
##categorical variables
par(mfrow = c(3, 3))
for (i in 1:ncol(cat_vars)) {
  barplot(table(cat_vars[,i]), main = colnames(cat_vars)[i],
          xlab = colnames(cat_vars)[i],col = c("skyblue", "salmon", "lightgreen",'yellow'))
}
###Bivariate Analysis
#num_var
par(mfrow = c(2, 3))
for (i in 1:ncol(numeric_vars)) {
  boxplot(numeric_vars[, i] ~ data$HeartDisease,
          main = paste("Boxplot of", colnames(numeric_vars)[i], "by Target"),
          xlab = "Target", ylab = colnames(numeric_vars)[i],
          col = c("lightblue", "lightgreen"))
}
par(mfrow = c(2, 3))
for (i in 1:ncol(numeric_vars)) {
  plot(density(numeric_vars[data$HeartDisease == 0, i]), 
       main = paste("Density Plot of", colnames(numeric_vars)[i], "for Target 0"),
       col = "lightblue", lwd = 2, xlim = range(numeric_vars[, i]))
  
  lines(density(numeric_vars[data$HeartDisease == 1, i]), col = "lightgreen", lwd = 2)
  
  legend("topright", legend = c("Target 0", "Target 1"), 
         col = c("lightblue", "lightgreen"), lwd = 2)
}
##categorical variables
#non runnare
par(mfrow = c(3, 2))
target_colors <- c("skyblue", "salmon")
for (i in 1:ncol(cat_vars[,-7])) { 
  contingency_table <- table(cat_vars[, i], cat_vars$HeartDisease)
  barplot(contingency_table,
          levels <- unique(cat_vars[, i]),
          colors <- rep(target_colors, length(levels) / length(target_colors)),
          main = paste("Barplot of", colnames(cat_vars)[i], "by Target"),
          xlab = colnames(cat_vars)[i], col =  c("skyblue", "salmon"),
          legend.text = TRUE, beside = TRUE)
}
##correlation
par(mfrow = c(1, 1))
pairs(numeric_vars, pch = 16)
cormat <- round(cor(numeric_vars),2)
corrplot(cormat, method = "number", type = "lower", 
         tl.col = "black", tl.srt = 45)
#partial correlation
library(igraph)
par(mfrow = c(1, 1))
S <- var(numeric_vars)
R <- -cov2cor(solve(S))
thr <- 0.1
G <- abs(R)>thr
diag(G) <- 0
Gi <- as(G, "igraph")
plot(Gi)
#contingecy table
calculate_chi_square <- function(data, target_var, categorical_var) {
  contingency_table <- table(data[[categorical_var]], data[[target_var]])
  chi_square_test <- chisq.test(contingency_table)
  return(list(contingency_table = contingency_table, chi_square_test = chi_square_test))
}
par(mfrow = c(2, 2)) 
for (col in colnames(cat_vars[,-7])) {
  result <- calculate_chi_square(data, "HeartDisease",col)
  cat("\n")
  cat(paste("Variable:", col, "\n"))
  print(result$contingency_table)
  print(result$chi_square_test)
}
par(mfrow = c(1, 1)) 
#data splitting
set.seed(123)
train_indices <- sample(1:nrow(data), 0.8 * nrow(data))
test_indices <- setdiff(1:nrow(data), train_indices)  # Indici rimanenti per il set di test

# Crea il set di training e il set di test utilizzando gli indici
train_set <- data[train_indices, ]
test_set <- data[test_indices, ]
##Logistic Regression
lr_model <- glm(HeartDisease ~ . , data=train_set, family=binomial)
lr_model_null <- glm(HeartDisease ~ +1 , data=train_set, family=binomial)
anova(lr_model_null, lr_model, test="Chisq")
summary(lr_model)
sort(vif(lr_model))
#from chi-test we can say that the difference of deviance between the null model and the dull model is significative,
# so the predictors brings information. However, in the full model lot of variables seems to be not significative
# and we procede with a feature selection. no collinearity problems seem to be raised
step_model_lr <- stepAIC(lr_model, direction = 'both')
summary(step_model_lr)
sort(vif(step_model_lr))
anova(lr_model, step_model_lr, test="Chisq")
#high values of p-value, therefore the reduced model and the full model are not signifcantly different, so the
#removed predictors are not significant
#prediction
fit_step_lr<- predict(step_model_lr, test_set, type = "response")
pred<- ifelse(fit_step_lr > 0.4, 1, 0)
#we tried different levels of threshold. reducing it would rise recall, that is the most usefull in our case.
#at the same time we want to keep the other metrics at decent levels
conf_matrix <- table(Predicted = pred, Actual = test_set$HeartDisease)

conf_matrix

precision <- round(conf_matrix[2, 2] / sum(conf_matrix[2, ]),2)
cat("Precision: ", precision, "\n")
accuracy <- round(sum(diag(conf_matrix)) / sum(conf_matrix),2)
cat("Accuracy: ", accuracy, "\n")
recall <- round(conf_matrix[2, 2] / sum(conf_matrix[, 2]),2)
cat("Recall: ", recall, "\n")
f1_score <- round(2 * (precision * recall) / (precision + recall),2)
cat("F1 Score: ", f1_score, "\n")
# Residuals plots
plot(step_model_lr, 
     which = c(2,3,4,5),
     col = as.numeric(train_set$HeartDisease),
     pch = as.numeric(train_set$HeartDisease),
)
######LDA
#shapiro test
response_0 <- train_set[train_set$HeartDisease == 0, ]
response_1 <- train_set[train_set$HeartDisease == 1, ]
shapiro_results_0 <- lapply(response_0[, numeric_vars, drop = FALSE], shapiro.test)
shapiro_results_1 <- lapply(response_1[, numeric_vars, drop = FALSE], shapiro.test)
# Funzione per stampare i risultati del test di Shapiro-Wilk
print_shapiro_results <- function(results, group_label) {
  for (name in names(results)) {
    cat("Shapiro-Wilk test for", name, "in group", group_label, ":\n")
    print(results[[name]])
    cat("\n")
  }
}
cat("Shapiro-Wilk Test Results for Response = 0:\n")
print_shapiro_results(shapiro_results_0, "Response = 0")
cat("Shapiro-Wilk Test Results for Response = 1:\n")
print_shapiro_results(shapiro_results_1, "Response = 1")
#all null hypothesis are rejected, no normal distribution in our variabile, we prooceed with LDA but we keep in mindi
#of thath
lda_model <- lda(HeartDisease ~ ., data = train_set)
lda_pred <-predict(lda_model, test_set[,-12])
predicted_classes <- lda_pred$class

# Risultati delle previsioni
actual_classes <- test_set$HeartDisease

# Calcolo delle metriche
confusion_matrix <- table(Predicted = predicted_classes, Actual = actual_classes)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)

# Precision, Recall, and F1-score
precision <- diag(confusion_matrix) / rowSums(confusion_matrix)
recall <- diag(confusion_matrix) / colSums(confusion_matrix)
f1_score <- 2 * (precision * recall) / (precision + recall)

# Stampa dei risultati
cat("Confusion Matrix:\n")
print(confusion_matrix)
cat("\n")

cat("Accuracy:", accuracy, "\n")

cat("Precision per classe:\n")
print(precision)
cat("\n")

cat("Recall per classe:\n")
print(recall)
cat("\n")

cat("F1-score per classe:\n")
print(f1_score)
cat("\n")

##QDA
qda_model <- qda(HeartDisease ~ ., data = train_set)
qda_pred <-predict(qda_model, test_set[,-12])
predicted_classes_qda <- qda_pred$class

# Calcolo delle metriche
confusion_matrix <- table(Predicted = predicted_classes_qda, Actual = actual_classes)

# Accuratezza
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)

# Precision, Recall, and F1-score
precision <- diag(confusion_matrix) / rowSums(confusion_matrix)
recall <- diag(confusion_matrix) / colSums(confusion_matrix)
f1_score <- 2 * (precision * recall) / (precision + recall)

# Stampa dei risultati
cat("Confusion Matrix:\n")
print(confusion_matrix)
cat("\n")

cat("Accuracy:", accuracy, "\n")

cat("Precision per classe:\n")
print(precision)
cat("\n")

cat("Recall per classe:\n")
print(recall)
cat("\n")

cat("F1-score per classe:\n")
print(f1_score)
cat("\n")

