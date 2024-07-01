## library
library(readr)
library(rmarkdown)
library(corrplot)
library(MASS)
library(car)

#####
##Data Pre-Processing
c_names <- c('Age', 'Gender','TB','DB','Alkphos','Sgpt','Sgot','TP','ALB','AG_ratio','Response')
data <- read_csv('Indian Liver Patient Dataset (ILPD).csv', col_names = c_names)
data$Response[data$Response == 2] <- 0
str(data)
#tranform gender and response variable into factor 
data$Gender <- as.factor(data$Gender)
data$Response <- as.factor(data$Response)
summary(data) #4 NA values in ag_ratio variable
#replacing NA's
data$AG_ratio[is.na(data$AG_ratio)]<- median(data$AG_ratio, na.rm = TRUE)
summary(data)
any(is.na(data))

#stand  log
num_vars <- sapply(data, is.numeric)
data_num <- data[, num_vars]
data_cat <- data[, !num_vars]
data_num_scaled <- scale(log(data_num))
data_num_scaled_df <- as.data.frame(data_num_scaled)
data <- cbind(data_num_scaled_df, data_cat)

##Logistic Regression
lr_model <- glm(Response ~ . , data=train_set, family=binomial)
lr_model_null <- glm(Response ~ +1 , data=train_set, family=binomial)
anova(lr_model_null, lr_model, test="Chisq")
summary(lr_model)
sort(vif(lr_model))
#from chi-test we can say that the difference of deviance between the null model and the dull model is significativr
#however, in the full model lot of variables seems to be not significative and the model suffers from multicollinearity
## stepwise selection
step_model_lr <- stepAIC(lr_model, direction = "both")
summary(step_model_lr)
sort(vif(step_model_lr))
anova(lr_model, step_model_lr, test="Chisq")
#high values of p-value, therefore the reduced model and the full model are not signifcantly different, so the
#removed predictors are not significant
#prediction
fit_step_lr<- predict(step_model_lr, test_set, type = "response")
pred<- ifelse(fit_step_lr > 0.4, 1, 0)
conf_matrix <- table(pred,test_set$Response)
conf_matrix
#
precision <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
cat("Precision: ", precision, "\n")
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
cat("Accuracy: ", accuracy, "\n")
recall <- conf_matrix[2, 2] / sum(conf_matrix[, 2])
cat("Recall: ", recall, "\n")
f1_score <- 2 * (precision * recall) / (precision + recall)
cat("F1 Score: ", f1_score, "\n")

# Residuals plots
plot(step_model_lr, 
     which = c(1,2,3,4,5),
     col = as.numeric(train_set$Response),
     pch = as.numeric(train_set$Response),
)

##non runnare

#leverage <- hatvalues(step_model_lr)
#high_leverage <- leverage > (3 * mean(leverage))
#which(high_leverage)
#train_set <- train_set[!(rownames(train_set) %in% c(115)), ]
