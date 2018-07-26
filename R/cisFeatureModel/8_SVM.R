# SVM

load('train.RData')
library(e1071)
library(class)
library(ROCR)
library(caret)

SVMCV <- function(data, formula, n, cost, gamma, n.site=1750) {
  # the return value - a list of prediction tables
  pred.list <- list()
  
  for (i in 1:n) {
    print(paste("Running CV:", i))
    
    all.index <- 1 : nrow(data)
    train.index <- ((i - 1) * n.site + 1) : (i * n.site);
    test.index <- all.index[-train.index]  # all sites except train sites
    train <- data[train.index,]
    test <- data[test.index,]
    
    # run svm
    svm <- svm(formula, train, cost = cost, gamma = gamma, probability=TRUE)
    pred <- predict(svm, test, decision.values=TRUE, probability=TRUE )
    prob<-data.frame(attr(pred, "probabilities")[,1])
    colnames(prob)<-c("prob")
    pred.table <- data.frame(true = test$class, beta = test$beta, pred = pred, prob=prob)
    pred.list[[i]] <- pred.table
  }
  
  pred.list
}

cost.v <- c(2^-1, 2^1, 2^3, 2^5, 2^7, 2^9)
gamma.v <- c(2^1, 2^-1, 2^-3, 2^-5, 2^-7, 2^-9)

i <- 1
for (cost in cost.v) {
  for (gamma in gamma.v) {
    print(paste("trial:", i, "cost:", cost, "gamma:", gamma))
    pred.list <- SVMCV(features, formula.all, 10, cost, gamma, 1750)
    for (j in 1:10) {
      file.name <- paste("svm/pred", i, j, sep = ".")
      write.csv(pred.list[[j]], file.name)
    }
    i <- i + 1
  }
}

evalSVM <- function(pred.list, n) {
  # vectors for 10 folds of accuracy, specificity, sensitivity and MCC
  acc.v <- NULL
  spec.v <- NULL
  sens.v <- NULL
  mcc.v <- NULL
  
  for (i in 1:n) {
    
    # count tp, tn, fp, fn
    pred <- pred.list[[i]]$pred
    true <- pred.list[[i]]$true
    tp <- as.numeric(length(which(pred == 1 & true == 1)))  # to prevent integer onverflow
    tn <- as.numeric(length(which(pred == 0 & true == 0)))
    fp <- as.numeric(length(which(pred == 1 & true == 0)))
    fn <- as.numeric(length(which(pred == 0 & true == 1)))
    
    # specificity
    spec <- tn / (tn + fp)
    spec.v <- c(spec.v, spec)
    
    # sensitivity
    sens <- tp / (tp + fn)
    sens.v <- c(sens.v, sens)
    
    # accuracy
    acc <- (tp + tn) / (tp + fp + tn + fn)
    acc.v <- c(acc.v, acc)
    
    # MCC
    mcc <- (tp * tn - fp * fn) / sqrt((tp + fn) * (tp + fp) * (tn + fp) * (tn + fn))
    mcc.v <- c(mcc.v, mcc)
  }
 
  
  data.frame(acc.v, spec.v, sens.v, mcc.v)
}

i = 1
acc.v = NULL
for (cost in cost.v) {
  for (gamma in gamma.v) {
    pred.list = list()
    for (j in 1:10) {
      file.name = paste("svm/pred", i, j, sep = ".")
      pred.list[[j]] = read.csv(file.name, header = TRUE, row.names = 1,
                                colClasses = c("factor", "factor", "double", "factor", "double"))
    }
    eval.svm = evalSVM(pred.list, 10)
    acc = mean(eval.svm$acc.v)
    print(paste("trial:", i, "cost:", cost, "gamma:", gamma, "acc:", acc))
    acc.v = c(acc.v, acc)
    i = i + 1;
  }
}

acc.m = matrix(acc.v, nrow = 6, ncol = 6)
rownames(acc.m) = gamma.v
colnames(acc.m) = cost.v
write.csv(acc.m, "svm/acc.csv")
save.image('SVM.RData')


