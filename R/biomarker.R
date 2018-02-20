#' Use Lasso to do feature selection
#'
#' @param df.input Row is sample, column is feature. Required
#' @param target.vec y vector. Required
#' @param nfolds glmnet CV nfolds
#' @param logisticRegression doing logistic regression or linear regression.
#' @param nRun number of glmnet runs
#' @param alpha same as in glmnet
#' @export
#' @examples
#' getSignatureFromMultipleGlmnet(df.input, target.vec)

getSignatureFromMultipleGlmnet <- function(df.input, target.vec, nfolds = 10, logisticRegression = TRUE, nRun=100, alpha = 1){

    #package
    require(glmnet)

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
            factor.variable.name <- c(factor.variable.name, colnames(df.input)[i])
            factor.col.index <- c(factor.col.index, i)
        }
    }

    #return(head(df.input))

    if (is.null(factor.col.index)){
      x <- as.matrix(df.input)
    } else{
      df.input.factor <- df.input[,factor.col.index]
      df.input.factor.onehot <- model.matrix( ~ .-1, data.frame(df.input.factor))
      df.input.factor.numeric <- df.input[,-factor.col.index]
      x <- as.matrix(cbind(df.input.factor.numeric, df.input.factor.onehot))
    }

    #return(x)


    #body


    featureDict <- list()
    featureNum <- c()
    weights.vec <- rep(0,ncol(x))
    for (i in seq(1,nRun,1)) {
        if (logisticRegression == FALSE){
            fit2 <- cv.glmnet(x, target.vec, alpha = alpha, nfolds = nfolds)
        } else{
            target.vec <- as.factor(target.vec)
            fit2 <- cv.glmnet(x, target.vec, alpha = alpha, family = "binomial", type.measure = "auc", nfolds = nfolds)
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
            featureDictInverse[[numTmpChr]] <- c(featureDictInverse[[numTmpChr]], names(featureDict)[i])
        }
        else {
            featureDictInverse[[numTmpChr]] <- c(names(featureDict)[i])
        }
    }

    numIndex <- sort(as.numeric(names(featureDictInverse)), decreasing = TRUE)
    featureSelectionFloor <- c()

    for (i in seq(1,length(numIndex),1)){
        numTmp <- numIndex[i]
        numTmpChr <- as.character(numTmp)
        featureSelectionFloor <- c(featureSelectionFloor, featureDictInverse[[numTmpChr]])
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
                          round(featureDict[[paste(featureSelectionFloor[i])]]/nRun, 2))
    }
    percent <- function(x, digits = 2, format = "f", ...) {
      paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
    }
    return.list[["selection_rate"]] <- percent(selection.rate)


    return.list[["feature_weights"]] <- weights.vec[match(featureSelectionFloor, colnames(x))]/nRun

    if (sum(selection.rate > 1) >0){
      feature.remove.index <- which(selection.rate > 1)
      return.list <- lapply(return.list, function(x) {x <- x[-feature.remove.index]} )
    }
    if (sum(startsWith(featureSelectionFloor,"others")) > 0){
      feature.remove.index <- which(startsWith(featureSelectionFloor,"others"))
      return.list <- lapply(return.list, function(x) {x <- x[-feature.remove.index]} )

    }

    # Here we get the features:
    return(return.list)

}
