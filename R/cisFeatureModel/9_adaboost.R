# AdaBoost

load('train.RData')
if (!require('adabag') install.packages('adabag', repos="http://cran.rstudio.com/")
library(adabag)
library(class)
library(ROCR)

AdaBoostCV <- function(data, formula, n, n.site=10000) {
  # the return value - a list of prediction tables
  pred.list <- list()

  for (i in 1:n) {
    print(paste("Running CV:", i))
    
    all.index <- 1 : nrow(data)
    train.index <- ((i - 1) * n.site + 1) : (i * n.site);
    test.index <- all.index[-train.index]  # all sites except train sites
    train <- data[train.index,]
    test <- data[test.index,]
    
    ada <- boosting(formula, train)
    prob <- predict(ada, test)
    score <- prob$prob[,2]
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
features$class <- as.factor(features$class)
formula.all <- paste(colnames(features)[c(7, 9, seq(14, 141))], collapse = '+')
formula.all <- paste("class ~ ", formula.all)
formula.all <- formula(formula.all)

pred.adaBoost <- AdaBoostCV(features, formula.all, 10, 1750)

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
        ROC <- prediction(pred.list[[i]]$score, pred.list[[i]]$label)

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


eval.adaBoost <- evalRF(pred.adaBoost, 10)
save(pred.adaBoost, eval.adaBoost, file='adaBoostRes.RData')
write.table(eval.adaBoost, file='adaBoostRes.csv', sep=',', row.names=T, col.names=F, quote=F)





