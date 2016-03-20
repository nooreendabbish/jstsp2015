library(glmnet)
library(lattice)
library(igraph)
library (corrplot)
library(reshape)
library(reshape2)
library(randomForest)
library(e1071)
library(pROC)

setwd("/Users/Brendan/Dropbox/Kaggle 2014/Data")

load("labels_train.Rda")
load("FNC_train.Rda")
load("SBM_train.Rda")



##Original Data
SBM_labels <- merge(labels_train, SBM_train, by="Id")
labeledSBM <- SBM_labels[,c(-1, -35)]
labeledSBM$Class <- factor(labeledSBM$Class)


###LP Score Function####
LP.Score.fun <- function(x,m){ 
  u <- (rank(x,ties.method = c("average")) - .5)/length(x) 
  m <- min(length(unique(u ))-1, m )
  S.mat <- as.matrix(poly(u ,df=m)) 
  return(as.matrix(scale(S.mat)))
}



###LP transformations of SBM Data###
#Combines LP transformations and labels into one data set#
LP_SBM <- function(x, y, m) {
  SBM_mat <- as.matrix(x)
  SBM.LP <- as.data.frame(y)
  for (i in 1:ncol(SBM_mat)) {   
    LP_mat <- as.data.frame(LP.Score.fun(SBM_mat[, i], m))
    names(LP_mat)[1:m] <- paste0(colnames(SBM_mat)[i], "_", 1:m)
    SBM.LP <- cbind(SBM.LP, LP_mat)
  }
  names(SBM.LP)[1] <- "Class"
  SBM.LP$Class <- as.factor(SBM.LP$Class)
  return(SBM.LP)
}

#Matrices of LP Transformations
sbmLP1 <- LP_SBM(SBM_train[,-1], labels_train$Class, 1)
sbmLP2 <- LP_SBM(SBM_train[,-1], labels_train$Class, 2)
sbmLP3 <- LP_SBM(SBM_train[,-1], labels_train$Class, 3)

###Lasso-Logistic Regression Function###

##Label Vector Must be First Column in Data Matrix and a Factor### 
##Returns Lasso-Logistic Accuracies over Given iteration and Resampling###
##Percent is the percentage of data to use as training data, the rest will
#be used for test
#iteration tell the funciton how many times to perform lasso logistic
LasLog <- function(X, percent, iteration) {
  bound <- floor((nrow( X )/ 100)*percent)
  acc <- NULL
  for (i in 1:iteration){
    Sample.tr <- X[sample(nrow(X)), ]
    train <- Sample.tr[1:bound, ] 
    test <- Sample.tr[(bound+1):nrow(Sample.tr), ]  
    log.cv.fit <- cv.glmnet(as.matrix(train[,-1]), train$Class, 
                       family = "binomial", alpha=1)
    log.predict <- as.factor(predict(log.cv.fit, as.matrix(test[,-1]),
                         type="class", s="lambda.min"))
    test[,1] <- as.numeric(levels(test[,1]))[test[,1]]
    log.predict <- as.numeric(levels(log.predict))[log.predict]
    acc[i] <- auc(test[,1], log.predict)
    }
acc <- as.data.frame(acc)
return(acc)
}

OrigLog <- LasLog(labeledSBM, 80, 200)
LPLog1 <- LasLog(sbmLP1, 80, 200)
LPLog2 <- LasLog(sbmLP2, 80, 200)
LPLog3 <- LasLog(sbmLP3, 80, 200)


####Random Forest Accuracy Function####
##Returns Accuracies of performing iterations of resampled data##


rf <- function(X, percent, iteration) {
  bound <- floor((nrow( X )/ 100)*percent)
  acc <- NULL
  for (i in 1:iteration){
    Sample.tr <- X[sample(nrow(X)), ]
    train <- Sample.tr[1:bound, ] 
    test <- Sample.tr[(bound+1):nrow(Sample.tr), ]  
    rf.fit <- randomForest(as.matrix(train[,-1]), train$Class, 
                            family = "binomial", alpha=1)
    rf.predict <- as.factor(predict(rf.fit, as.matrix(test[,-1]),
                                     type="class"))
    test[,1] <- as.numeric(levels(test[,1]))[test[,1]]
    rf.predict <- as.numeric(levels(rf.predict))[rf.predict]
    acc[i] <- auc(test[,1], rf.predict)
  }
  acc <- as.data.frame(acc)
  return(acc)
}

Orig.rf <- rf(labeledSBM, 80, 200)
LP.rf1 <- rf(sbmLP1, 80, 200)
LP.rf2 <- rf(sbmLP2, 80, 200)
LP.rf3 <- rf(sbmLP3, 80, 200)

##SVM Accuracy Function####
#X must be data matrix with factored labels as first column
#k determines value for k-fold cross validation
SVMcross <- function(X, percent, k, iteration) {
  accuracy <- NULL
  for (i in 1:iteration) {
    X.tr.samp <- X[sample(nrow(X)), ]
    bound <- floor((nrow(X.tr.samp)/100)*80)
    data.tr <- X.tr.samp[1:bound, ]
    data.test <- X.tr.samp[(bound+1):nrow(X.tr.samp),]
    data.tune <- tune.svm(data.tr[,-1], data.tr[,1],
                          gamma = 2^(-15:3), cost = 2^(-3:15), 
                          tune.control(cross = k))
    bestGamma <- data.tune$best.parameters[[1]]
    bestC <- data.tune$best.parameters[[2]]
    data.svm <- svm(data.tr[,-1], data.tr[,1], gamma = bestGamma,
                    cost = bestC)
    data.predict <- as.factor(predict(data.svm, data.test[,-1],
                            type="class"))
    data.predict <- as.numeric(levels(data.predict))[data.predict]
    data.test[,1] <- as.numeric(levels(data.test[,1]))[data.test[,1]]
    accuracy[i] <- auc(data.test[,1], data.predict)
  }
  accuracy <- as.data.frame(accuracy)
  return(accuracy)
}

Orig.SVM <- SVMcross(labeledSBM, 80, 10, 100)
LP.svm3 <- SVMcross(sbmLP3, 80, 10, 100)


##Building Boxplot Comparisons

library(ggplot2)
acc.rf <- cbind(LP.rf3, Orig.rf)
names(acc.rf) <- c("LP", "Original")
acc.rf <- melt(acc.rf)
acc.rf$Model <- rep("RF", length(acc.rf))

acc.log <- cbind(LPlog3, OrigLog)
names(acc.log) <- c("LP", "Original")
acc.log <- melt(acc.log)
acc.log$Model <- rep("Logistic", length(acc.log))


acc.SVM <- cbind(LP.svm3, Orig.SVM)
names(acc.SVM) <- c("LP", "Original")
acc.SVM <- melt(acc.SVM)
acc.SVM$Model <- rep("SVM", length(acc.SVM))

acc.box <- rbind(acc.rf, acc.log, acc.SVM)


ggplot(data = acc.box, aes(x=Model, y=value)) +  geom_boxplot(aes(fill=variable))


###Boxplots Comparing different Lasso-Logistics###

logaccs <- cbind(OrigLog, LPLog1, LPlog2, LPlog3)
names(logaccs) <- c("Original", "m=1", "m=2", "m=3")
logaccs <- melt(logaccs)

ggplot(data = logaccs, aes(x=variable, y=value)) +  geom_boxplot(aes(fill=variable))

####Boxplots Comparing Random Forests###
rfaccs <- cbind(Orig.rf, LP.rf1, LP.rf2, LP.rf3)
names(rfaccs) <- c("Original", "m=1", "m=2", "m=3")
rfaccs <- melt(rfaccs)

ggplot(data = rfaccs, aes(x=variable, y=value)) +  geom_boxplot(aes(fill=variable))