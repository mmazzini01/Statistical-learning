###### Library
library(rmarkdown)
library(corrplot)
library(MASS)
library(car)
library(pROC)
library(glmnet)
library(class)

##### Pre-processing
data <- read.csv('heart_data.csv', header = T, sep = ',')

# Change categorical into factor
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
numeric_vars <- data[, sapply(data, is.numeric)]
cat_vars <- data[, sapply(data, is.factor)]
# from the summary we can see that there in non missing value but in the variable cholesterol and restingBp there values 
# that are equal to zero (172 and 1 respectively), this could be missing value or errors in the recording since in the first
# case the sierum cholesterol can't be qual to zero, especially for 172 patients. in the second case blood pressure in 
# human alive cant be zero, and the patient is alivr since heartj records are different from zero. we procced by substituing thi
# zero values with the median that is less sensitive than the mean from extrem observations
data[RestingBP == 0, c("MaxHR", "RestingECG",'HeartDisease')]
data[Cholesterol == 0, c("MaxHR", "RestingECG",'HeartDisease')]
detach(data)
data$Cholesterol[data$Cholesterol==0] <- NA
data$RestingBP[data$RestingBP==0] <- NA
data$Cholesterol[is.na(data$Cholesterol)] <- median(data$Cholesterol, na.rm = TRUE)
data$RestingBP[is.na(data$RestingBP)] <- median(data$RestingBP, na.rm = TRUE)
any(is.na(data))
summary(data)# data are ready to be processed

##### Explorative Data Analysis
### Univariate Analysis
## Numerical variables
par(mfrow = c(2, 3))
for (i in 1:ncol(numeric_vars)) {
  hist(numeric_vars[, i], main = colnames(numeric_vars)[i],
       xlab = colnames(numeric_vars)[i],col = 'lightblue', freq = FALSE)
  dens <- density(numeric_vars[, i], na.rm=TRUE, adjust=1.25)
  lines(dens, col = "black", lwd = 1)
}
# age and maxHR are nearly normal, resting BP ha alcuni valori motlo alti. cholesterol seems to be bimodal. and oldpeak
# has concentration between zero values with lot of observaztion with higher values

par(mfrow = c(2, 3))
for (i in 1:ncol(numeric_vars)) {
  boxplot(numeric_vars[, i], main = colnames(numeric_vars)[i],
       xlab = colnames(numeric_vars)[i],col = 'lightblue', freq = FALSE)
}
# oldpeack and cholesterol seems to have outliers in the right tail as well as restingBp

## Categorical variables
par(mfrow = c(3, 3))
for (i in 1:ncol(cat_vars)) {
  barplot(table(cat_vars[,i]), main = colnames(cat_vars)[i],
          xlab = colnames(cat_vars)[i],col = c("skyblue", "salmon", "lightgreen",'yellow'))
}
# male and female seemst to be quite unbalanced like fastingbs. asynthomatic patient are more numerosi than chestpaiin patient
# down st slope are much way less than other two levels. our response variable is quite balanced.

### Bivariate Analysis
## Numerical variables
par(mfrow = c(2, 3))
for (i in 1:ncol(numeric_vars)) {
  boxplot(numeric_vars[, i] ~ data$HeartDisease,
          main = paste("Boxplot of", colnames(numeric_vars)[i], "by Target"),
          xlab = "Target", ylab = colnames(numeric_vars)[i],
          col = c("skyblue", "salmon"))
}
# Age: Individuals with heart disease tend to be slightly older.
# RestingBP*: Slightly higher resting blood pressure is observed in individuals with heart disease.
# Cholesterol*: Higher cholesterol levels are more common among those with heart disease.
# MaxHR*: Lower maximum heart rate is seen in individuals with heart disease.
# Oldpeak*: Higher oldpeak values are associated with heart disease.


par(mfrow = c(2, 3))
for (i in 1:ncol(numeric_vars)) {
  plot(density(numeric_vars[data$HeartDisease == 0, i]), 
       main = paste("Density Plot ", colnames(numeric_vars)[i], "vs Target"),
       col = "skyblue", lwd = 2, xlim = range(numeric_vars[, i]))
  
  lines(density(numeric_vars[data$HeartDisease == 1, i]), col = "salmon", lwd = 2)
}
plot.new() 
legend("center", legend = c("Non-Disiased", " Disiased"), 
       col = c("skyblue", "salmon"), lwd = 2)
# from the plot we can see that diseased people tend to have higher age. resting Bp seems to have no difference.
# maxhr is higher in patient with disease. Patients with Disease have higher values of oldpeak and it seems that in the 
# last one group there bimodal distribution with two subgroups with picchi of 0 old peack and 2 respectively. there could 
# be a subgroup of cholesterol people with higher values. althougth, when we analyze cholesterol we have to remind that 
# a bunch of values are replaced with the median

## Categorical variables
par(mfrow = c(3, 2), mar = c(2, 2, 4,8))
for (i in 1:(ncol(cat_vars) - 1)) { 
  contingency_table <- table(cat_vars[, i], cat_vars$HeartDisease)
  num_levels <- nlevels(cat_vars[, i])
  if (num_levels == 2) {
    colors <- c('skyblue', 'salmon')
  } else if (num_levels == 3) {
    colors <- c('skyblue', 'salmon', 'lightgreen')
  } else if (num_levels == 4) {
    colors <- c('skyblue', 'salmon', 'lightgreen', 'yellow')
  } else {
    colors <- rainbow(num_levels) 
  }
  
  barplot(contingency_table,
          main = paste("Barplot of", colnames(cat_vars)[i], "by HeartDisease"),
          col = colors,
          xlab = colnames(cat_vars)[i],
          beside = TRUE)
  legend("topright", legend = levels(cat_vars[, i]), fill = colors, cex = 0.9, inset = c(0, 0))
}
# More males have heart disease compared to females, and more females do not have heart disease compared to males.
# most individuals with heart disease have ASY chest pain, while ATA and NAP chest pains are more common in individuals without heart disease. this seems to be controintuitivo but we will discuss later on our analysis
# Individuals with heart disease have a higher proportion of fasting blood sugar ≥ 120 mg/dl compared to those without heart disease.
# Individuals without heart disease have a higher proportion of normal resting ECG results, while individuals with heart disease show a higher proportion of ST results. but it seems to be not much significant
# Exercise-induced angina is more common in individuals with heart disease compared to those without it.
# Individuals with heart disease have a higher proportion of flat ST slopes, while those without heart disease have a higher proportion of up ST slopes.


# Contingecy table
calculate_chi_square <- function(data, target_var, categorical_var) {
  contingency_table <- table(data[[categorical_var]], data[[target_var]])
  chi_square_test <- chisq.test(contingency_table)
  return(list(contingency_table = contingency_table, chi_square_test = chi_square_test))
}
for (col in colnames(cat_vars[,-7])) {
  result <- calculate_chi_square(data, "HeartDisease",col)
  cat("\n")
  cat(paste("Variable:", col, "\n"))
  print(result$contingency_table)
  print(result$chi_square_test)
}
# null hypothesis are rejected. chi-squared test confirms associations between our categorical variables and response variable

## Correlation
par(mfrow = c(1, 1))
pairs(numeric_vars, pch = 16)
cormat <- round(cor(numeric_vars),2)
corrplot(cormat, method = "number", type = "lower", 
         tl.col = "black", tl.srt = 45)
# age seems to be 'central' from the other numeric variables, they are correlated with age, the only exception
# is given by cholesterol 
par(mfrow = c(1, 1))

## Data splitting
set.seed(123)
train_indices <- sample(1:nrow(data), 0.8 * nrow(data))
test_indices <- setdiff(1:nrow(data), train_indices) 

train_set <- data[train_indices, ]
test_set <- data[test_indices, ]

# since our numeric data are represented in different units of measures we proceed with a standardization 
# Scaling
numeric_vars_train <- train_set[, sapply(train_set, is.numeric)]
cat_vars_train <- train_set[, sapply(train_set, is.factor)]

numeric_vars_test <- test_set[, sapply(test_set, is.numeric)]
cat_vars_test <- test_set[, sapply(test_set, is.factor)]

data_num_scaled_train <- scale(numeric_vars_train)
data_num_scaled_df_train <- as.data.frame(data_num_scaled_train)
train_set <- cbind(data_num_scaled_df_train, cat_vars_train)

data_num_scaled_test <- scale(numeric_vars_test)
data_num_scaled_df_test <- as.data.frame(data_num_scaled_test)
test_set <- cbind(data_num_scaled_df_test, cat_vars_test)

# Logistic Regression
lr_model <- glm(HeartDisease ~ . , data=train_set, family=binomial)
lr_model_null <- glm(HeartDisease ~ +1 , data=train_set, family=binomial)
anova(lr_model_null, lr_model, test="Chisq")
summary(lr_model)
sort(vif(lr_model))
# from chi-test we can say that the difference of deviance between the null model and the dull model is significative,
# so the predictors brings information. However, in the full model lot of variables seems to be not significative
# and we procede with a feature selection. no collinearity problems seem to be raised
step_model_lr <- stepAIC(lr_model, direction = 'both')
summary(step_model_lr)
sort(vif(step_model_lr))
anova(lr_model, step_model_lr, test="Chisq")
# high values of p-value, therefore the reduced model and the full model are not signifcantly different, so the
# removed predictors are not significant

# Prediction
pred_lr_prob <- predict(step_model_lr, test_set, type = "response")
pred_lr <- ifelse(pred_lr_prob > 0.5, 1, 0)


# Confusion matrix
compute_confusion_matrix <- function(Predicted, Actual) {
  conf_matrix <- table(Predicted, Actual)
  conf_df <- as.data.frame(matrix(0, nrow = 2, ncol = 2))
  colnames(conf_df) <- c("True Non-Disiased", "True Disiased")
  rownames(conf_df) <- c("Pred. Non-Disiased", "Pred. Disiased")
  conf_df[, 1:2] <- t(conf_matrix)
  conf_df["Total",] <- colSums(conf_df)
  conf_df <- cbind(conf_df, Total = rowSums(conf_df))
  return(conf_df)
}

conf_matrix_lr <- compute_confusion_matrix(test_set$HeartDisease, pred_lr)
conf_matrix_lr

# Metrics
compute_metrics <- function(conf_matrix, Actual, Predicted_prob) {
  acc <- round((conf_matrix[1,1]+conf_matrix[2,2])/conf_matrix[3,3],3)
  prec <- round(conf_matrix[2,2]/conf_matrix[2,3],3)
  rec <- round(conf_matrix[2,2]/conf_matrix[3,2],3)
  spec <- round(conf_matrix[1,1]/conf_matrix[3,1],3)
  type_1 <- round(conf_matrix[2,1]/conf_matrix[3,1],3)
  f1_score <- round(2 * (prec * rec) / (prec + rec),3)
  
  cat("Accuracy:", acc, "\n")
  cat("Precision:", prec, "\n")
  cat("Recall:", rec, "\n")
  cat("Specificity:", spec, "\n")
  cat("Type 1 error:", type_1, "\n")
  cat("F1 Score: ", f1_score, "\n")
  
  roc_out <- roc(Actual, as.numeric(Predicted_prob))
  plot(roc_out, print.auc = TRUE, 
       xlab="False positive rate(Type 1 error)", 
       ylab="True positive rate(Recall)", 
       legacy.axes = TRUE)
}

compute_metrics(conf_matrix_lr, test_set$HeartDisease, pred_lr_prob)
# we tried different levels of threshold. reducing it would rise recall, that is the most usefull in our case.
# at the same time we want to keep the other metrics at decent levels. 0.5 is the best since lower levels give 
# negligible advantages in terms of recall but significant worst results in precision

##Logistic-Regression diagnostic
# Residuals plots
plot(step_model_lr, 
     which = c(4,5),
     col = as.numeric(train_set$HeartDisease),
     pch = as.numeric(train_set$HeartDisease),
)
# from cook distance and leverage plot we can see that there 3 influential points, we try to eliminate them 
# and build again our model to compare results

lev_points <- c(680,557,786)
data[lev_points,]
train_set_clean <- train_set[!rownames(train_set) %in% lev_points, ]
lr_model <- glm(HeartDisease ~ . , data=train_set_clean, family=binomial)
step_model_lr <- stepAIC(lr_model, direction = 'both')
pred_lr_prob <- predict(step_model_lr, test_set, type = "response")
pred_lr <- ifelse(pred_lr_prob > 0.5, 1, 0)
conf_matrix_lr <- compute_confusion_matrix(test_set$HeartDisease, pred_lr)
conf_matrix_lr
compute_metrics(conf_matrix_lr, test_set$HeartDisease, pred_lr_prob)
#model achieve better quality metrics, especially precision
# and specificity gain significant points with same recall's level
summary(step_model_lr)#now the model is composed by also the cholesterol variable
# Sex , ChestPainType, Fasting Blood Sugar, Exercise-Induced Angina, and ST Slope seems to be the most effective 
# predictors. the coefficients represent the variations of the logit by an additional increase of unit of that variable
# by keeping the others fixed.Cholesterol is not significant but we keep it since there might be interaction effects between cholesterol 
# and other variables that is significant even if the individual effect is not. let's check
updated_model <- update(step_model_lr, . ~ . - Cholesterol)
summary(updated_model)
# the AIC confirms that model with cholesterol is better 503.91 vs 504.11
# odds ratio
lr_coefficients <- coef(step_model_lr)
# Calculate the odds ratios
odds_ratios <- exp(lr_coefficients)
odds_ratios
# risk factors for cvd disease are: having one additional unit of age increment the prob of beeing diseased by 38%
# oldpeak: 44%
# sex: beeing male 482%
# asyntomathic chest pain type is considered a risk factor since the other categories bring a decrease of 85-6% for 
#ATA and NAP, and 62% for TA. this is controintuitive but it could be by the fact that our data are clearly unbalanced
#and most of the patients have asyntomatic chest pain type
# fastingangina: yes increase 261%
# st-slope_ flat increase 238% , slope up is not consiered a risk factor, prob decrease by 78%

## Ridge and Lasso Regression
# X and y
X_train <- model.matrix(HeartDisease~., train_set)[,-1]
y_train <- as.numeric(as.character(train_set$HeartDisease))

X_test <- model.matrix(HeartDisease~., test_set)[,-1]
y_test <- as.numeric(as.character(test_set$HeartDisease))

cumpute_AIC_R_L <- function(pred, coef, y){
  n <- length(y)
  rss <- sum((y - pred)^2)
  edf <- sum(coef != 0)
  aic_ridge <- n * log(rss/n) + 2 * edf
  cat("AIC: ", aic_ridge, "\n")
}

# Ridge regression
ridge_cv <- cv.glmnet(X_train, y_train, alpha = 0, family = "binomial", type.measure = "deviance", nfolds = 10)
plot(ridge_cv)
lambda = ridge_cv$lambda.min
cat("The value for the minimum lambda is ", lambda)


pred_ridge_prob <- predict(ridge_cv, X_test, type = "response", s = lambda)
pred_ridge <- ifelse(pred_ridge_prob > 0.5, 1, 0)

# Confusion Matrix
conf_matrix_ridge <- compute_confusion_matrix(y_test, pred_ridge)
conf_matrix_ridge

# Metrics
compute_metrics(conf_matrix_ridge, y_test, pred_ridge_prob)
ridge_coef <- coef(ridge_cv, s = "lambda.min")
ridge_coef
lr_coefficients
cumpute_AIC_R_L(pred_ridge, ridge_coef, y_test)
# predictors shrunk near zero: RestingBP, RestingEcg. other preditors with smaller absolute values are 
# age, cholesterol and maxHR

# Lasso regression
lasso_cv <- cv.glmnet(X_train, y_train, alpha = 1, family = "binomial", type.measure = "deviance", nfolds = 10)
plot(lasso_cv)
lambda = lasso_cv$lambda.min
cat("The value for the minimum lambda is ", lambda)
pred_lasso_prob <- predict(lasso_cv, X_test, type = "response", s = lambda)
pred_lasso <- ifelse(pred_lasso_prob > 0.5, 1, 0)

# Confusion Matrix
conf_matrix_lasso <- compute_confusion_matrix(y_test, pred_lasso)
conf_matrix_lasso

# Metrics
compute_metrics(conf_matrix_lasso, y_test, pred_lasso_prob)

lasso_coef <- coef(lasso_cv, s = "lambda.min")
lasso_coef
cumpute_AIC_R_L(pred_lasso, lasso_coef, y_test)
# best logistic regression model seen so far looking at the metrics results.
# parameters shrunk to zero: RestingBP, RestingECG
# Cholesterol and MaxHR more shrunked than in ridge, less importance to them.

# LDA
# null hypothesis are rejected, in almost all group of variable. no normal distribution in our variable, 
# we proced with LDA but we keep in mind of that

cumpute_AIC_L_Q <- function(model_type, model, data, X, y){
  pred <- predict(model, data)$posterior
  log_likelihood <- sum(log(rowSums(pred * (model.matrix(~ HeartDisease - 1, data = data)))))
  if(model_type == 'LDA') {
    n_params <- ncol(X) + (ncol(X) * (ncol(X) + 1)) / 2 + length(unique(y)) - 1
  }
  else if (model_type == 'QDA'){
    n_classes <- length(unique(y))
    n_params <- n_classes * (ncol(X) + (ncol(X) * (ncol(X) + 1)) / 2) + n_classes - 1
  }
  aic <- -2 * log_likelihood + 2 * n_params
  cat("AIC: ", aic, "\n")
}

lda_model <- lda(HeartDisease ~ ., data = train_set)
pred_lda_prob <- predict(lda_model, test_set,type ='response')$posterior[,2]
pred_lda <- as.factor(ifelse(pred_lda_prob > 0.5, 1, 0))
lda_model$scaling

# Each coefficient indicates the importance and direction of the influence of a predictor variable on class discrimination.
# A positive coefficient indicates that an increase in the variable is associated with an increase in the probability of belonging to a particular class,
#while a negative coefficient indicates the opposite.
#Strongly Influencing Variables:
# SexM, ExerciseAnginaY, and ST_SlopeFlat have significant positive coefficients
#Weakly Influencing Variables:
#RestingBP, RestingECGNormal, and RestingECGST .they have coefficients very close to zero, indicating that their influence
#on class discrimination is very weak

# Confusion Matrix
conf_matrix_lda <- compute_confusion_matrix(test_set$HeartDisease, pred_lda)
conf_matrix_lda

# Metrics
compute_metrics(conf_matrix_lda, test_set$HeartDisease, pred_lda_prob)
# model very valid same metrics score of lasso

cumpute_AIC_L_Q('LDA',lda_model, train_set, train_set[, 1:11], train_set[, 12])

# QDA
qda_model <- qda(HeartDisease ~ ., data = train_set)
pred_qda_prob <- predict(qda_model, test_set,type ='response')$posterior[,2]
pred_qda <- as.factor(ifelse(pred_qda_prob>0.5,1,0))

# Confusion Matrix
conf_matrix_qda <- compute_confusion_matrix(test_set$HeartDisease, pred_qda)
conf_matrix_qda

# Metrics
compute_metrics(conf_matrix_qda, test_set$HeartDisease, pred_qda_prob)

cumpute_AIC_L_Q('QDA',qda_model, train_set, train_set[, 1:11], train_set[, 12])
# KNN
#k <- 
#pred_knn <- knn(train = X_train, test = X_test, cl = y_train, k = k)

# Confusion Matrix
#conf_matrix_knn <- compute_confusion_matrix(y_test, pred_knn)
#conf_matrix_knn

# Metrics
#compute_metrics(conf_matrix_knn, y_test, as.numeric(pred_knn))

# our data seems to be ambigous since it seems that having chest pain reduce the probability to have hearth disease, so be asyntomatic 
# is considered a risk factor that is not very likely ro see in reality. so we want to explore better our data to undartand this phenomena

# first we divide chestPainType into two group. the fisrt is composed by asyntomatic patient, while the second  by patient that have a chest pain
# with no distinction of its type. we want to investigate how this two groups interact with others risk factor variables

patient_with_chest_pain <- data[data$ChestPainType != 'ASY', ]
patient_without_chest_pain <- data[data$ChestPainType == 'ASY', ] 

par(mfrow = c(1, 2))
boxplot(patient_with_chest_pain$Age, patient_without_chest_pain$Age,
        names = c("Pain", "No-Pain"),
        main = "Age",
        col = c("skyblue", "salmon"))

boxplot(patient_with_chest_pain$Oldpeak, patient_without_chest_pain$Oldpeak,
        names = c("Pain", "No-Pain"),
        main = "Oldpeak ",
        col = c("skyblue", "salmon"))
# contegency table for cat. variable
#ST_slope
ST_table_with_chest_pain <- table(patient_with_chest_pain$ST_Slope)
ST_table_without_chest_pain <- table(patient_without_chest_pain$ST_Slope)
ST_contingency_table <- rbind(ST_table_with_chest_pain, ST_table_without_chest_pain)
rownames(ST_contingency_table) <- c("With Chest Pain", "Without Chest Pain")
#Exercise_Angina
EA_table_with_chest_pain <- table(patient_with_chest_pain$ExerciseAngina)
EA_table_without_chest_pain <- table(patient_without_chest_pain$ExerciseAngina)
EA_contingency_table <- rbind(EA_table_with_chest_pain, EA_table_without_chest_pain)
rownames(EA_contingency_table) <- c("With Chest Pain", "Without Chest Pain")
ST_chi_square_test <- chisq.test(ST_contingency_table)
EA_chi_square_test <- chisq.test(EA_contingency_table)
# from the boxplot and contigecy table we can see that asynthomatic patient tend to have higher values in risk factor. this couold be
# brought by two reasons. the first is that our data are biased. the second is a medical one, for example most of hearth disease have no
# syntom or maybe there is others medical issue that we cannot grasp
# to check if our data are biased we try to remove this variable and we see if this variable really bring additional information or other
# variables themselves could do the work alone. we try our best model lasso regression
par(mfrow = c(1, 1))

# Lasso regression less parameter
X_train_reduced <- model.matrix(HeartDisease~., train_set[,-c(7)])[,-1]
y_train <- as.numeric(as.character(train_set$HeartDisease))

X_test_reduced <- model.matrix(HeartDisease~., test_set[,-c(7)])[,-1]
y_test <- as.numeric(as.character(test_set$HeartDisease))
lasso_red <- cv.glmnet(X_train_reduced, y_train, alpha = 1, family = "binomial", type.measure = "deviance", nfolds = 10)
plot(lasso_cv)
lambda = lasso_cv$lambda.min
cat("The value for the minimum lambda is ", lambda)
pred_lasso_red__prob <- predict(lasso_red, X_test_reduced, type = "response", s = lambda)
pred_lasso_red <- ifelse(pred_lasso_prob > 0.5, 1, 0)

# Confusion Matrix
conf_matrix_lasso_reduced <- compute_confusion_matrix(y_test, pred_lasso_red)
conf_matrix_lasso_reduced

# Metrics
compute_metrics(conf_matrix_lasso_reduced, y_test, pred_lasso_red_prob)

# it seems that the model achieve the same results also with that variable, this suggest to us that our data have some bais
# since at fisrt chestpaintipe is considere one of the most significative variable for our prediction. but in reality
# no additional information are brought since the model remains pretty much the same. (guardare confounding effect)

lasso_red_coef <- coef(lasso_cv, s = "lambda.min")
lasso_red_coef

#now we saw that most of the variables considered risk factors are the ones coming from a cardiac stress test. this could be very
#usefull  information in every day life because this test could be use alone for predicting patient heart disesase. now , taking into
#account this aspect we want to build a model with only these variable and compare the results with our full model. if the variation 
# is not significative this coul be use as a point of riferimento for meds anc other exams like analisi del sangue could be trascured.
# we try again our best model

# Lasso regression with demographic and stress test parameters
X_train_stress <- model.matrix(HeartDisease~., train_set[,-c(3,7,8)])[,-1]
y_train <- as.numeric(as.character(train_set$HeartDisease))

X_test_stress <- model.matrix(HeartDisease~., test_set[,-c(3,7,8)])[,-1]
y_test <- as.numeric(as.character(test_set$HeartDisease))
lasso_stress <- cv.glmnet(X_train_stress, y_train, alpha = 1, family = "binomial", type.measure = "deviance", nfolds = 10)
plot(lasso_cv)
lambda = lasso_cv$lambda.min
cat("The value for the minimum lambda is ", lambda)
pred_lasso_stress__prob <- predict(lasso_stress, X_test_stress, type = "response", s = lambda)
pred_lasso_stress <- ifelse(pred_lasso_stress__prob > 0.5, 1, 0)

# Confusion Matrix
conf_matrix_lasso_stress<- compute_confusion_matrix(y_test, pred_lasso_stress)
conf_matrix_lasso_stress

# Metrics
compute_metrics(conf_matrix_lasso_stress, y_test, pred_lasso_stress_prob)

lasso_stress_coef <- coef(lasso_cv, s = "lambda.min")
lasso_stress_coef
# valid metrics, this model could be use 

# LDA model with stress test parameters
lda_stress_model <- lda(HeartDisease ~ ., data = train_set[,-c(3,7,8)])
pred_stress_lda_prob <- predict(lda_stress_model, test_set,type ='response')$posterior[,2]
pred_stress_lda <- as.factor(ifelse(pred_stress_lda_prob > 0.5, 1, 0))
lda_stress_model$scaling

# Confusion Matrix
conf_matrix_lda_stress <- compute_confusion_matrix(test_set$HeartDisease, pred_stress_lda)
conf_matrix_lda_stress

# Metrics
compute_metrics(conf_matrix_lda_stress, test_set$HeartDisease, pred_stress_lda)
# also this model very valid, slightly better recall but worst precision
