#' Use Lasso to do feature selection
#'
#' @param df.input Row is sample, column is feature. Required
#' @param target.vec y vector. Required
#' @param nfolds glmnet CV nfolds
#' @param logisticRegression doing logistic regression or linear regression.
#' @param nRun number of glmnet runs
#' @param alpha same as in glmnet
#' @import glmnet
#' @return signature
#' @export
#' @examples
#' data('iris')
#' getSignatureFromMultipleGlmnet(iris[,1:4],
#' c(rep(1,100), rep(0,50)), nfolds = 3, nRun = 10)

getSignatureFromMultipleGlmnet <- function(df.input,
    target.vec,
    nfolds = 10,
    logisticRegression = TRUE,
    nRun=100,
    alpha = 1){

    # make target numeric
    if (!is.numeric(target.vec)){
        target.vec[target.vec == target.vec[1]] <- 1
        target.vec[target.vec != 1] <- 0
    }

    # factorize catogorical variables
    factor.col.index <- c()
    factor.variable.name <- c()
    for (i in 1:ncol(df.input)){
        if (is.character(df.input[,i]) | is.factor(df.input[,i])){
            df.input[,i] <- as.factor(df.input[,i])
            factor.variable.name <- c(factor.variable.name,
    colnames(df.input)[i])
            factor.col.index <- c(factor.col.index, i)
        }
    }

    #return(head(df.input))

    if (is.null(factor.col.index)){
        x <- as.matrix(df.input)
    } else{
        df.input.factor <- df.input[,factor.col.index]
        df.input.factor.onehot <- model.matrix( ~ .-1,
        data.frame(df.input.factor))
        df.input.factor.numeric <- df.input[,-factor.col.index]
        x <- as.matrix(cbind(df.input.factor.numeric,
        df.input.factor.onehot))
    }

    #return(x)


    #body


    featureDict <- list()
    featureNum <- c()
    weights.vec <- rep(0,ncol(x))
    for (i in seq(1,nRun,1)) {
        if (logisticRegression == FALSE){
            fit2 <- suppressWarnings(cv.glmnet(x, target.vec,
    alpha = alpha,
    nfolds = nfolds))
        } else{
            target.vec <- as.factor(target.vec)
            fit2 <- suppressWarnings(cv.glmnet(x, target.vec,
    alpha = alpha,
    family = "binomial",
    type.measure = "auc",
    nfolds = nfolds))
        }
        weights.mat <- as.matrix(coef(fit2, s = "lambda.min"))
        weights.vec <- (weights.vec + weights.mat[-1,1])

        tmp_vec <- as.vector((coef(fit2, s="lambda.min") != 0))
        if (sum(tmp_vec) <= 1){
    next
        }
        featureFromGlmnet <- colnames(x)[tmp_vec[-1]]
        featureNum <- c(featureNum, length(featureFromGlmnet))
        for (k in seq(1,length(featureFromGlmnet),1)){
            gene <- featureFromGlmnet[k]
            if (gene %in% names(featureDict)){
                featureDict[[gene]] <- featureDict[[gene]] + 1
            }
            else{
                if (is.na(gene) == FALSE){
                    featureDict[[gene]] <- 1
                }

            }
        }
    }
    #print(featureDict)

    featureSelectionComplete <- names(featureDict)

    numFloor <- floor(mean(featureNum))

    featureDictInverse <- list()
    for (i in seq(1,length(featureDict),1)){
        numTmp <- featureDict[[i]]
        #print(numTmp)
        numTmpChr <- as.character(numTmp)
        if (numTmp %in% names(featureDictInverse)){
            featureDictInverse[[numTmpChr]] <-
    c(featureDictInverse[[numTmpChr]],
    names(featureDict)[i])
        }
        else {
            featureDictInverse[[numTmpChr]] <- c(names(featureDict)[i])
        }
    }

    numIndex <- sort(as.numeric(names(featureDictInverse)),
    decreasing = TRUE)
    featureSelectionFloor <- c()

    for (i in seq(1,length(numIndex),1)){
        numTmp <- numIndex[i]
        numTmpChr <- as.character(numTmp)
        featureSelectionFloor <- c(featureSelectionFloor,
    featureDictInverse[[numTmpChr]])
        if (length(featureSelectionFloor) > numFloor) {
            break
        }
    }

    return.list <- list()
    return.list[["feature"]] <- featureSelectionFloor
    return.list[["counts"]] <- featureDict
    return.list[["counts.inverse"]] <- featureDictInverse
    return.list[["weights"]] <- weights.vec

    selection.rate <- c()
    for (i in 1:length(featureSelectionFloor)){
    selection.rate <- c(selection.rate,
    round(featureDict[[
    paste(featureSelectionFloor[i])]]/nRun, 2))
    }
    percent <- function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
    }
    return.list[["selection_rate"]] <- percent(selection.rate)


    return.list[["feature_weights"]] <- weights.vec[match(featureSelectionFloor,
    colnames(x))]/nRun

    if (sum(selection.rate > 1) >0){
    feature.remove.index <- which(selection.rate > 1)
    return.list <- lapply(return.list, function(x)
        {x <- x[-feature.remove.index]} )
    }
    if (sum(startsWith(featureSelectionFloor,"others")) > 0){
    feature.remove.index <- which(startsWith(featureSelectionFloor,
    "others"))
    return.list <- lapply(return.list, function(x)
        {x <- x[-feature.remove.index]} )

    }

    # Here we get the features:
    return(return.list)

}

#' Do bootstrap and LOOCV
#'
#' @param df Row is sample, column is feature. Required
#' @param targetVec y vector. Required
#' @param nboot number of BOOTSTRAP
#' @importFrom gmodels ci
#' @return bootstrap loocv result dataframe
#' @export
#' @examples
#' data('iris')
#' Bootstrap_LOOCV_LR_AUC(iris[,1:4],
#' c(rep(1,100), rep(0,50)), nboot = 3)
Bootstrap_LOOCV_LR_AUC <- function(df, targetVec, nboot=50){
    output.auc.vec <- c()
    output.other.df <- NULL

    # make target numeric
    if (!is.numeric(targetVec)){
        targetVec[targetVec == targetVec[1]] <- 1
        targetVec[targetVec != 1] <- 0
    }
    targetVec <- as.numeric(targetVec)
    auc.vec <- c()
    for (i in 1:nboot){
        index.boot <- sample(1:nrow(df), nrow(df), replace = TRUE)
        df.tmp <- df[index.boot,]
        auc.vec <- c(auc.vec, LOOAUC_simple_multiple_noplot_one_df(df.tmp,
        targetVec[index.boot]))
    }

    result.type <- c("AUC Estimate","CI lower","CI upper","Std. Error")
    output.df <- data.frame(result.type, ci(auc.vec))
    colnames(output.df) <- c("Type", "Value")
    return(output.df)

}




#' LOOCV
#'
#' @param df Row is sample, column is feature. Required
#' @param targetVec y vector. Required
#' @importFrom ROCR prediction performance
#' @return mean auc
#' @export
#' @examples
#' data('iris')
#' LOOAUC_simple_multiple_noplot_one_df(iris[,1:4],
#' c(rep(1,100), rep(0,50)))

LOOAUC_simple_multiple_noplot_one_df <- function(df, targetVec){
    auc.vec <- c()
    nSample <- nrow(df)
    vecProbTmp <- c()
    testPredictionClassVec <- c()
    for (j in 1:nSample){
        train = as.matrix(df[-j,])
        test = matrix(as.numeric(df[j,]), nrow = 1)
        fit <- glmnet(train, targetVec[-j], family = "binomial")
        testProb <- predict(fit,type="response", newx = test, s = 0)
        vecProbTmp <- c(vecProbTmp, testProb)
        testPredictionClassVec[j] <- predict(fit,type="class",
        newx = test, s = 0)
    }
    loo.pred = prediction(vecProbTmp, targetVec)
    loo.perf = performance(loo.pred,"tpr","fpr")
    auc <- performance(loo.pred,"auc")
    auc <- unlist(slot(auc, "y.values"))
    aucRound <- round(auc,3)
    auc.vec <- c(auc.vec, aucRound)

    return(mean(auc.vec))
}



#' LOOCV with ROC curve
#'
#' @param df Row is sample, column is feature. Required
#' @param targetVec y vector. Required
#' @importFrom ROCR prediction performance
#' @return the ROC
#' @export
#' @examples
#' data('iris')
#' LOOAUC_simple_multiple_one_df(iris[,1:4],
#' c(rep(1,100), rep(0,50)))

LOOAUC_simple_multiple_one_df <- function(df, targetVec){

    auc.vec <- c()
    nSample <- nrow(df)
    vecProbTmp <- c()
    testPredictionClassVec <- c()

    # make target numeric
    if (!is.numeric(targetVec)){
        targetVec[targetVec == targetVec[1]] <- 1
        targetVec[targetVec != 1] <- 0
    }
    targetVec <- as.numeric(targetVec)

    for (j in 1:nSample){
        train = as.matrix(df[-j,])
        test = matrix(as.numeric(df[j,]), nrow = 1)
        fit <- glmnet(train, targetVec[-j], family = "binomial")
        testProb <- predict(fit,type="response", newx = test, s = 0)
        vecProbTmp <- c(vecProbTmp, testProb)
        testPredictionClassVec[j] <- predict(fit,
        type="class",
        newx = test,
        s = 0)
    }
    loo.pred.plot = prediction(vecProbTmp, targetVec)
    loo.perf.plot = performance(loo.pred.plot,"tpr","fpr")
    return(list(loo.pred.plot = loo.pred.plot,
    loo.perf.plot = loo.perf.plot,
    testPredictionClassVec = testPredictionClassVec))
}






