#####
## library
library(readr)
library(rmarkdown)
library(caret)

#####
##Data Pre-Processing
c_names <- c('Age', 'Gender','TB','DB','Alkphos','Sgpt','Sgot','TP','ALB','AVG_ratio','Response')
data <- read_csv('Indian Liver Patient Dataset (ILPD).csv', col_names = c_names)
data$Response[data$Response == 2] <- 0
str(data)
#tranform gender and response variable into factor 
data$Gender <- as.factor(data$Gender)
data$Response <- as.factor(data$Response)
summary(data) #4 NA values in avg_ratio variable
#replacing NA's
data$AVG_ratio[is.na(data$AVG_ratio)]<- median(data$AVG_ratio, na.rm = TRUE)
summary(data)
any(is.na(data))
#train test splitting
#we want to respect the proportion of response variable in both train and test set
#we use a stratified sampling
train_indices <- createDataPartition(data$Response, p = 0.75, list = FALSE)
data_train <- data[train_indices, ]
data_test <- data[-train_indices, ]
#check proportion
data_prop<-round(prop.table(table(data$Response)),3)
train_prop<-round(prop.table(table(data_train$Response)),3)
test_prop<-round(prop.table(table(data_test$Response)),3)
cat('data: ',data_prop,'train: ',train_prop,'test: ',test_prop)

#####
## Esplorative analysis


