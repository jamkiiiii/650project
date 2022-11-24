######################################################################### 
#-----------------------------------------------------------------------#
# Created on:    2022/11/8                                              #
# Updated on:    2022/11/20                                             #
# EMAIL:         dengfy@umich.edu                                       #
#                Umich  Feiyang Deng.                                   #
#########################################################################

#### Data Analysis Script ####

rm(list=ls())
library('readxl')
library('ggplot2')
library('tidyverse')
library('mice')
library('VIM')


### import data
setwd("E:/Biostat Study/BIOSTAT 650/Group Project")
set.seed(123)

inputFile <- "Depression Data.xls"
raw_data <- read_excel(inputFile, sheet=1, na='NA')
raw_data <- as.data.frame(raw_data)
summary(raw_data)

## Visualize missing values
md.pattern(raw_data)
aggr(raw_data, prop=F, numbers=T)


#########################################################
## Descriptive analysis for complete and incomplete data
#########################################################


## add indicator to identify complete rows
na_rows <- which(rowSums(is.na(raw_data)) > 0)
raw_data$complete <- 1
raw_data$complete[na_rows] <- 0
raw_data$complete <- as.factor(raw_data$complete)


### Descriptive analysis of AGE4, TOTMMSE4, TOTIADL4, TOTADL4, EE46
## All continuous variables (except outcome) and one categorical predictor of interest
## Boxplot
ggplot(raw_data, aes(x=complete, y=AGE4, fill=complete)) +
  geom_boxplot()
ggplot(raw_data, aes(x=complete, y=TOTMMSE4, fill=complete)) +
  geom_boxplot()
ggplot(raw_data, aes(x=complete, y=TOTIADL4, fill=complete)) +
  geom_boxplot()
ggplot(raw_data, aes(x=complete, y=TOTADL4, fill=complete)) +
  geom_boxplot()
ggplot(raw_data, aes(x=complete, y=EE46, fill=complete)) +
  geom_boxplot()

## Density
ggplot(raw_data, aes(x=AGE4, fill=complete)) +
  geom_density(alpha=.25)
ggplot(raw_data, aes(x=TOTMMSE4, fill=complete)) +
  geom_density(alpha=.25)
ggplot(raw_data, aes(x=TOTIADL4, fill=complete)) +
  geom_density(alpha=.25)
ggplot(raw_data, aes(x=TOTADL4, fill=complete)) +
  geom_density(alpha=.25)
## Histogram
ggplot() +
  geom_histogram(data=raw_data[-na_rows,], 
                 binwidth = 0.5,
                 aes(x=EE46, y=..count../(1682-length(na_rows))),
                 fill = 13, alpha = 0.7) +
  geom_histogram(data=raw_data[na_rows,], 
                 binwidth = 0.5,
                 aes(x=EE46, y=-..count../sum(1-is.na(raw_data$EE46[na_rows]))),
                 fill = 10, alpha = 0.7) +
  ylab('Incomplete data VS Complete data') +
  ggtitle('Histogram of EE46') +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate(geom = "label", 
           x = 4, y = 0.4, 
           label = "EE46 of Complete data",
           fill = 12, 
           color = "white", hjust = 1) +
  annotate(geom = "label", 
           x = 4, y = -0.3, 
           label = "EE46 of Inomplete data",
           fill = 14, 
           color = "white", hjust = 1)


## test if same distribution AGE4, TOTMMSE4, TOTIADL4, TOTADL4, EE46
ks.test(raw_data$AGE4[-na_rows], raw_data$AGE4[na_rows], na.omit=TRUE)
ks.test(raw_data$TOTMMSE4[-na_rows], raw_data$TOTMMSE4[na_rows], na.omit=TRUE)
ks.test(raw_data$TOTIADL4[-na_rows], raw_data$TOTIADL4[na_rows], na.omit=TRUE)
ks.test(raw_data$TOTADL4[-na_rows], raw_data$TOTADL4[na_rows], na.omit=TRUE)
chisq.test(raw_data$EE46, raw_data$complete)




#################
## Imputation
#################
library("DMwR2")
library("mice")
library("Hmisc")
# column index
# Categorical variables : GRADE, USBORN, marstat, nkids4, health4, khyper41, mdiab41, nfrac41, u43s, cc43, ee46, hha4, oo49lang, Male
# indeces : GRADE(2), USBORN(3), MARSTAT4(5), NKIDS4(6), HEALTH4(7), KHYPER41(8),
contin_idx <- c(4,12,14,15)
outcome_idx <- c(13)
cate_na_idx <- c(2,6,8,9,10,11,16,17,18)# index of categorical with missing values
cate_full_idx <- c(3,5,7,19,20)
other_idx <- c(1,21)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

Imputation <- function(raw_data, contin_idx, outcome_idx, cate_na_idx, cate_full_idx){
  ## this function only work for our raw_data ##
  
  ### Continuous data imputation
  ## 1. KNN
  contin_data_knn <- knnImputation(raw_data[-c(1,21)])[,contin_idx-1]
  # summary(contin_data_knn)
  
  ## 2. Mice
  tmp <- mice(raw_data[-c(1,21)],m=5,maxit=50,meth='pmm',seed=500)
  contin_data_mice <- complete(tmp,1)[,contin_idx-1]
  # summary(contin_data_mice)
  
  ## 3. Mean
  contin_data_mean <- raw_data[,contin_idx]
  for(i in c(1:4)){
    contin_data_mean[,i] <- impute(contin_data_mean[,i], mean)  # 均值替代
  }
  # summary(contin_data_mean)
  
  ## 4. Median
  contin_data_median <- raw_data[,contin_idx]
  for(i in c(1:4)){
    contin_data_median[,i] <- impute(contin_data_median[,i], median)  # 中位数替代
  }
  # summary(contin_data_median)
  
  ### Categorical data imputation
  # Mode
  cate_na_data <- raw_data[,cate_na_idx]
  
  mode_list <- c()
  for( i in cate_na_idx){
    mode_list <- c(mode_list, getmode(raw_data[,i]))
  }
  mode_list
  for( i in 1:ncol(cate_na_data)){
    ind=which(is.na(cate_na_data[,i]))
    cate_na_data[ind,i]=mode_list[i]
  }
  
  # head(cate_na_data,10)
  
  ### Outcome NA
  outcome_na_idx <- which(is.na(raw_data[outcome_idx]))
  
  return(list(outcome=raw_data[outcome_idx], 
              cate_full=raw_data[cate_full_idx], 
              cate_na=cate_na_data, 
              contin_data_knn=contin_data_knn,
              contin_data_mice=contin_data_mice,
              contin_data_mean=contin_data_mean,
              contin_data_median=contin_data_median,
              outcome_na_idx=outcome_na_idx))
}

impute_result <- Imputation(raw_data, contin_idx, outcome_idx, cate_na_idx, cate_full_idx)

summary(impute_result$cate_na)
summary(impute_result$contin_data_knn)
#### combine data after imputing
new_data <- cbind(impute_result$outcome, 
                  impute_result$cate_full, 
                  impute_result$cate_na, 
                  impute_result$contin_data_knn)[-impute_result$outcome_na_idx]


summary(new_data)
### Scatter plot for Depresion vs. predictor of interest
ggplot(new_data, aes(x = TOTIADL4, y = CESDTOT4)) + 
  geom_point(color = "blue")+
  geom_smooth(method = lm, color = "red", fill="#69b3a2", se = TRUE)

ggplot(new_data, aes(x = EE46, y = CESDTOT4)) + 
  geom_point(color = "blue")+
  geom_smooth(method = lm, color = "red", fill="#69b3a2", se = TRUE)












