# Random forest model

load('feature2803.RData')

if (!require('randomForest')) install.packages('randomForest', repos="http://cran.rstudio.com/")
if (!require('caret')) install.packages('caret', repos="http://cran.rstudio.com/")
if (!require('class')) install.packages('class', repos="http://cran.rstudio.com/")
if (!require('ROCR')) install.packages('ROCR', repos="http://cran.rstudio.com/")

library(randomForest)
library(caret)
library(class)
library(ROCR)

shuffle.idx <- sample.int(nrow(features))
features <- features[shuffle.idx, ]

RFCV <- function(data, formula, n.fold, train.number) {
    pred.list <- list()
    for (i in 1:n.fold) {
        print(paste("CV fold", i))

        all.idx <- 1:nrow(data)
        train.idx <- ((i - 1) * train.number + 1) : (i * train.number)
        test.idx <- all.idx[-train.idx]
        train <- data[train.idx,]
        test <- data[test.idx,]

        rf <- randomForest(formula, train, ntree=1000)
        prob <- predict(rf, test, type='prob')
        score <- prob[,2]
        label <- test$class

        pred <- score
        pred[pred >= 0.5] <- 1
        pred[pred < 0.5] <- 0
        pred <- as.factor(pred)

        pred.table <- data.frame(label=label, beta=test$beta, 
                                 pred=pred, score=score)
        pred.list[[i]] <- pred.table
    }
    pred.list
}

formula.all <- paste(colnames(features)[c(7, 9, seq(14, 141))], collapse = '+')
formula.all <- paste("as.factor(class) ~ ", formula.all)
formula.all <- formula(formula.all)

pred.RF <- RFCV(features, formula.all, 10, 1750)

# Evaluation metrics

evalRF <- function(pred.list, n.fold) {
    acc.v <- NULL
    auc.v <- NULL
    spec.v <- NULL
    sens.v <- NULL
    mcc.v <- NULL
    tp.v <- NULL
    fp.v <- NULL
    tn.v <- NULL
    fn.v <- NULL
    for (i in seq(n.fold)) {
        ROC <- prediction(pred.list[[i]][,4], pred.list[[i]][,1])

        auc.pref <- performance(ROC, "auc")
        auc <- auc.pref@y.values[[1]]
        auc.v <- c(auc.v, auc)

        pred <- pred.list[[i]]$pred
        label <- pred.list[[i]]$label
        tp <- as.numeric(length(which(pred==1 & label==1)))
        tn <- as.numeric(length(which(pred==0 & label==0)))
        fp <- as.numeric(length(which(pred==1 & label==0)))
        fn <- as.numeric(length(which(pred==0 & label==1)))
        tp.v <- c(tp.v, tp)
        tn.v <- c(tn.v, tn)
        fp.v <- c(fp.v, fp)
        fn.v <- c(fn.v, fn)

        spec <- tn / (tn + fp)
        spec.v <- c(spec.v, spec)

        sens <- tp / (tp + fn)
        sens.v <- c(sens.v, sens)

        acc <- (tp + tn) / (tp + fp + tn + fn)
        acc.v <- c(acc.v, acc)

        mcc <- (tp * tn - fp * fn) / sqrt((tp + fn) * (tp + fp) * (tn + fp) * (tn + fn))
        mcc.v <- c(mcc.v, mcc)
    }
    data.frame(acc.v, auc.v, spec.v, sens.v, mcc.v, tp.v, fp.v, tn.v, fn.v)
}

eval.RF <- evalRF(pred.RF, 10)
save(pred.RF, eval.RF, file='RFresults.RData')
write.table(eval.RF, file='RFres.csv', sep=',', row.names=T, col.names=F, quote=F)
save(features, formula.all, file='train.RData')

