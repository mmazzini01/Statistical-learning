#####
## library
library(readr)
library(rmarkdown)
library(corrplot)


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
#train test splitting
#we want to respect the proportion of response variable in both train and test set
#we use a stratified sampling
# Proporzione del training set
train_proportion <- 0.75
train_set <- data.frame()
test_set <- data.frame()
for (level in unique(data$Response)) {
  subset_data <- data[data$Response == level, ]
  train_size <- round(nrow(subset_data) * train_proportion)
  train_indices <- sample(1:nrow(subset_data), train_size)
  train_subset <- subset_data[train_indices, ]
  test_subset <- subset_data[-train_indices, ]
  train_set <- rbind(train_set, train_subset)
  test_set <- rbind(test_set, test_subset)
}
#check proportion
train_prop <- round(prop.table(table(train_set$Response)),3)
test_prop <- round(prop.table(table(test_set$Response)),3)
data_prop<- round(prop.table(table(data$Response)),3)
cat('data: ',data_prop,'train: ',train_prop,'test: ',test_prop)

#train_indices <- createDataPartition(data$Response, p = 0.75, list = FALSE)
#data_train <- data[train_indices, ]
#data_test <- data[-train_indices, ]
#check proportion
#data_prop<-round(prop.table(table(data$Response)),3)
#train_prop<-round(prop.table(table(data_train$Response)),3)
#test_prop<-round(prop.table(table(data_test$Response)),3)
#cat('data: ',data_prop,'train: ',train_prop,'test: ',test_prop)

#####
## Esplorative analysis
##Univariate Analysis
#Response
barplot(table(data$Response),
        main = "Response frequency",
        xlab = "Response",
        ylab = "Frequency",
        col = c("skyblue", "salmon"),
        ylim = c(0, max(table(data$Response)) + 1))
#response variable quite unbalanced as seen before
#Gender
barplot(table(data$Gender),
        main = "Gender Bar Plot",
        xlab = "Category",
        ylab = "Frequency",
        col = c("skyblue", "salmon", "lightgreen"))
#Age
par(mfrow = c(1, 2))
boxplot(data$Age,
        main = "Age Boxplot",
        ylab = "Age",
        col = "lightblue")
hist(data$Age,
     main = "Age Histogram",
     xlab = "Values",
     col = "lightblue",
     border = "black")
#age variable seems to be quite symmetric
par(mfrow = c(1, 1))
qqnorm(data$Age)
qqline(data$Age, col = "red", lwd = 2)
#TB
par(mfrow = c(1, 2))
boxplot(data$TB,
        main = "TB Boxplot",
        ylab = "TB",
        col = "lightblue")
hist(data$TB,
     main = "TB Histogram",
     xlab = "Values",
     col = "lightblue",
     border = "black")
#DB
par(mfrow = c(1, 2))
boxplot(data$DB,
        main = "DB Boxplot",
        ylab = "DB",
        col = "lightblue")
hist(data$DB,
     main = "DB Histogram",
     xlab = "Values",
     col = "lightblue",
     border = "black")
#Alkphos
par(mfrow = c(1, 2))
boxplot(data$Alkphos,
        main = "Alkphos Boxplot",
        ylab = "Alkphos",
        col = "lightblue")
hist(data$Alkphos,
     main = "Alkphos Histogram",
     xlab = "Values",
     col = "lightblue",
     border = "black")
#Sgpt
par(mfrow = c(1, 2))
boxplot(data$Sgpt,
        main = "Sgpt Boxplot",
        ylab = "Sgpt",
        col = "lightblue")
hist(data$Sgpt,
     main = "Sgpt Histogram",
     xlab = "Values",
     col = "lightblue",
     border = "black")
#sgot
par(mfrow = c(1, 2))
boxplot(data$Sgot,
        main = "Sgot Boxplot",
        ylab = "Sgot",
        col = "lightblue")
hist(data$Sgot,
     main = "Sgot Histogram",
     xlab = "Values",
     col = "lightblue",
     border = "black")
#TP
par(mfrow = c(1, 2))
boxplot(data$TP,
        main = "TP Boxplot",
        ylab = "TP",
        col = "lightblue")
hist(data$Sgpt,
     main = "TP Histogram",
     xlab = "Values",
     col = "lightblue",
     border = "black",
     log = "y")
#ALB
par(mfrow = c(1, 2))
boxplot(data$ALB,
         main = "ALB Boxplot",
         ylab = "ALB",
         col = "lightblue")
hist(data$ALB,
     main = "ALB Histogram",
     xlab = "Values",
     col = "lightblue",
     border = "black")
#ALB seems to be quite symmetric
par(mfrow = c(1, 1))
qqnorm(data$ALB)
qqline(data$ALB, col = "red", lwd = 2)
#AG_ratio
par(mfrow = c(1, 2))
boxplot(data$AG_ratio,
        main = "AG_ratio Boxplot",
        ylab = "AG_ratio",
        col = "lightblue")
hist(data$AG_ratio,
     main = "AG_ratio Histogram",
     xlab = "Values",
     col = "lightblue",
     border = "black")
par(mfrow = c(1, 1))
#Age and ALB quite symmetric, both likely a normal distribution with some deviation in the tails
#variable TB,DB, Alkphos, Sgpt, Sgot and TP have strong positive skewnwess
#AG_ratio positive skeweness
##Bivariate analysis
#Gender and Response
contingency_table <- table(data$Gender, data$Response)  
contingency_table
chi_squared_test <- chisq.test(contingency_table)
chi_squared_test


#Function to boxplot
boxplot_response <- function(data, variable) {
  data_0 <- data[[variable]][data$Response == 0]
  data_1 <- data[[variable]][data$Response == 1]
  boxplot(data_0, data_1,
          main = paste(variable, "Boxplot"),
          ylab = variable,
          xlab = "Response",
          names = c("0", "1"),
          col = c("lightblue", "red"))
}
#Function to density plot
density_response <- function(data, variable) {
  dens_0 <- density(data[[variable]][data$Response == 0])
  dens_1 <- density(data[[variable]][data$Response == 1])
  plot(dens_0, col = "lightblue", main = paste(variable, "Density Plot"), xlab = variable, ylab = "Density", lwd = 2)
  lines(dens_1, col = "red", lwd = 2)
  legend("topright", legend = c("0", "1"), fill = c("lightblue", "red"))
}

#Age and Response
par(mfrow = c(1, 2))
boxplot_response(data, "Age")
density_response(data, "Age")
#The age do not seem to be a good predictor of the response variable

#TB and Response
par(mfrow = c(1, 2))
boxplot_response(data, "TB")
density_response(data, "TB")
#The TB seems to be a good predictor of the response variable

#DB and Response
par(mfrow = c(1, 2))
boxplot_response(data, "DB")
density_response(data, "DB")
#The DB seems to be a good predictor of the response variable

#Alkphos and Response
par(mfrow = c(1, 2))
boxplot_response(data, "Alkphos")
density_response(data, "Alkphos")
#The Alkphos seems to be a good predictor of the response variable

#Sgpt and Response
par(mfrow = c(1, 2))
boxplot_response(data, "Sgpt")
density_response(data, "Sgpt")
#The Sgpt seems to be a good predictor of the response variable

#Sgot and Response
par(mfrow = c(1, 2))
boxplot_response(data, "Sgot")
density_response(data, "Sgot")
#The Sgot seems to be a good predictor of the response variable

#TP and Response
par(mfrow = c(1, 2))
boxplot_response(data, "TP")
density_response(data, "TP")
#The TP do not seem to be a good predictor of the response variable

#ALB and Response
par(mfrow = c(1, 2))
boxplot_response(data, "ALB")
density_response(data, "ALB")
#The ALB do not seem to be a good predictor of the response variable

#AVG_ratio and Response
par(mfrow = c(1, 2))
boxplot_response(data, "AG_ratio")
density_response(data, "AG_ratio")
#The AVG_ratio do not seem to be a good predictor of the response variable
##Correlation
num_vars <- c_names[-c(2,11)]
pairs(data[, num_vars], pch = 16)
cormat <- round(cor(data[,num_vars]),2)
corrplot(cormat, method = "circle", type = "lower", 
         tl.col = "black", tl.srt = 45)
#TB and DB are highly correlated
#Sgpt and Sgot are highly correlated
#TP and ALB are highly correlated
#AVG_ratio and ALB are highly correlated
