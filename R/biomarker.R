#' Use Lasso to do feature selection
#'
#' @param df.input Row is feature, column is sample. Required
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
    
    if (!is.numeric(target.vec)){
        target.vec[target.vec == target.vec[1]] <- 1
        target.vec[target.vec != 1] <- 0
    }
    
    #body
    x <- t(as.matrix(df.input))
    
    featureDict <- list()
    featureNum <- c()
    weights.vec <- rep(0,nrow(df.input))
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
        featureFromGlmnet <- rownames(df.input)[tmp_vec]
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
    # Here we get the features:
    return(return.list)
    
}
