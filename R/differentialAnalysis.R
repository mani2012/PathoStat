################################################################################
#' Convert phyloseq OTU count data into DGEList for edgeR package
#'
#' Further details.
#'
#' @param physeq (Required).  A \code{\link{phyloseq-class}} or
#' an \code{\link{otu_table-class}} object.
#' The latter is only appropriate if \code{group} argument is also a
#' vector or factor with length equal to \code{nsamples(physeq)}.
#'
#' @param group (Required). A character vector or factor giving the experimental
#' group/condition for each sample/library. Alternatively, you may provide
#' the name of a sample variable. This name should be among the output of
#' \code{sample_variables(physeq)}, in which case
#' \code{get_variable(physeq, group)}
#' would return either a character vector or factor.
#' This is passed on to \code{\link[edgeR]{DGEList}},
#' and you may find further details or examples in its documentation.
#'
#' @param method (Optional). The label of the
#' edgeR-implemented normalization to use.
#' See \code{\link[edgeR]{calcNormFactors}} for supported options and details.
#' The default option is \code{"RLE"}, which is a scaling factor method
#' proposed by Anders and Huber (2010).
#' At time of writing, the \link[edgeR]{edgeR} package supported
#' the following options to the \code{method} argument:
#'
#' \code{c("TMM", "RLE", "upperquartile", "none")}.
#'
#' @param ... Additional arguments passed on to \code{\link[edgeR]{DGEList}}
#'
#' @import edgeR phyloseq
#' @return dispersion
#' @export
#' @examples
#' data_dir_test <- system.file("data", package = "PathoStat")
#' pstat_test <- loadPstat(indir=data_dir_test,
#' infileName="pstat_data_2_L1.rda")
#' phyloseq_to_edgeR(pstat_test, group="Sex")


phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){

    # Enforce orientation.
    if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
    x = as(otu_table(physeq), "matrix")
    # Add one to protect against overflow, log(0) issues.
    x = x + 1
    # Check `group` argument
    if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
        # Assume that group was a sample variable name (must be categorical)
        group = get_variable(physeq, group)
    }
    # Define gene annotations (`genes`) as tax_table
    taxonomy = tax_table(physeq, errorIfNULL=FALSE)
    if( !is.null(taxonomy) ){
        taxonomy = data.frame(as(taxonomy, "matrix"))
    }
    # Now turn into a DGEList
    y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
    # Calculate the normalization factors
    z = calcNormFactors(y, method=method)
    # Check for division by zero inside `calcNormFactors`
    if( !all(is.finite(z$samples$norm.factors)) ){
        stop("Something wrong with edgeR::calcNormFactors on this data,
    non-finite $norm.factors, consider changing `method` argument")
    }
    # Estimate dispersions
    return(estimateTagwiseDisp(estimateCommonDisp(z)))
}
################################################################################



#' Mann-whitney test for a dataframe
#'
#' @param df Input data object that contains the data to be tested. Required
#' @param label.vec.num The target binary condition. Required
#' @param pvalue.cutoff choose p-value cut-off
#' @return df.output object
#' @export
#' @examples
#' data('iris')
#' Wilcox_Test_df(t(iris[,1:4]),
#' c(rep(1,100), rep(0,50)))

Wilcox_Test_df <- function(df, label.vec.num, pvalue.cutoff = 0.05) {
    df.output <- NULL
    #save raw values
    label.vec.save <- unique(label.vec.num)

    # transform label into 1 and 0
    label.vec.num[label.vec.num == unique(label.vec.num)[1]] <- 1
    label.vec.num[label.vec.num != 1] <- 0
    for (i in 1:nrow(df)){
    # remove zero-variance rows
    if (sum(df[i,] == 1) == length(label.vec.num) |
    sum(df[i,] == 0) == length(label.vec.num)){
        next
    }
    tmp.result <- wilcox.test(df[i,which(label.vec.num == 1)],
    df[i,which(label.vec.num == 0)], correct=FALSE)
    if (tmp.result$p.value <= pvalue.cutoff & rownames(df)[i] != "others"){
        num.1 <- sum((df[i,which(label.vec.num == 1)] > 0))
        num.2 <- sum((df[i,which(label.vec.num == 0)] > 0))
        df.output <- rbind(df.output, c(rownames(df)[i],
        round(as.numeric(tmp.result$p.value), 4), num.1, num.2))
    }
    }
    if (is.null(df.output)){
        return(0)
    }
    df.output.prevalence <- percent(round((as.numeric(df.output[,3])+
    as.numeric(df.output[,4]))/ncol(df),4))
    df.output <- cbind(df.output, df.output.prevalence)
    colnames(df.output) <- c("Name", "FDR",
    label.vec.save[1],label.vec.save[2], "prevalence")
    # FDR adjustment
    df.output[,2] <- p.adjust(df.output[,2],
    method = "fdr",
    n = length(df.output[,2]))
    df.output <- df.output[order(df.output[,2]),]
    return(df.output)
}






#' transform cpm counts to presence-absence matrix
#'
#' @param df Input data object that contains the data to be tested.
#' Required
#' @return df.output object
#' @export
#' @examples
#' GET_PAM(data.frame(a = c(1,3,0), b = c(0,0.1,2)))

GET_PAM <- function(df) {
    for (i in 1:nrow(df)){
        df[i,] <- as.numeric(df[i,] > 0)
    }
    return(df)
}



#' Given PAM and disease/control annotation,
#' do Chi-square test for each row of PAM
#'
#' @param pam Input data object that contains the data to be tested. Required
#' @param label.vec.num The target binary condition. Required
#' @param pvalue.cutoff choose p-value cut-off
#' @return df.output object
#' @export
#' @examples
#' tmp <- matrix(rbinom(12,1,0.5), nrow = 3)
#' rownames(tmp) <- c("a", "b", "c")
#' Chisq_Test_Pam(tmp, c(1,1,0,0))


Chisq_Test_Pam <- function(pam, label.vec.num, pvalue.cutoff = 0.05) {
    df.output <- NULL

    #save raw values
    label.vec.save <- unique(label.vec.num)

    # transform label into 1 and 0
    label.vec.num[label.vec.num == unique(label.vec.num)[1]] <- 1
    label.vec.num[label.vec.num != 1] <- 0


    for (i in 1:nrow(pam)){
    # remove zero-variance rows
    if (sum(pam[i,] == 1) == length(label.vec.num) |
        sum(pam[i,] == 0) == length(label.vec.num)){
        next
    }
    tmp.result <- chisq.test(pam[i,], label.vec.num, correct=FALSE)
    if (tmp.result$p.value <= pvalue.cutoff  & rownames(pam)[i] != "others"){
        num.1 <- sum(pam[i,] == 1 & label.vec.num == 1)
        num.2 <- sum(pam[i,] == 1 & label.vec.num == 0)
        df.output <- rbind(df.output, c(rownames(pam)[i],
        round(as.numeric(tmp.result$p.value), 4), num.1, num.2))
    }
    }
    if (is.null(df.output)){
        return(0)
    }
    df.output.prevalence <- percent(round((as.numeric(df.output[,3])+
    as.numeric(df.output[,4]))/ncol(pam),4))
    df.output <- cbind(df.output, df.output.prevalence)
    colnames(df.output) <- c("Name", "FDR",
    label.vec.save[1], label.vec.save[2], "prevalence")
    # FDR adjustment
    df.output[,2] <- p.adjust(df.output[,2],
    method = "fdr", n = length(df.output[,2]))
    df.output <- df.output[order(df.output[,2]),]
    return(df.output)
}



#' Given PAM and disease/control annotation,
#' do Chi-square test for each row of PAM
#'
#' @param pam Input data object that contains the data to be tested. Required
#' @param label.vec.num The target binary condition. Required
#' @param pvalue.cutoff choose p-value cut-off
#' @return df.output object
#' @export
#' @examples
#' tmp <- matrix(rbinom(12,1,0.5), nrow = 3)
#' rownames(tmp) <- c("a", "b", "c")
#' Fisher_Test_Pam(tmp, c(1,1,0,0))

Fisher_Test_Pam <- function(pam, label.vec.num, pvalue.cutoff = 0.05) {
    df.output <- NULL

    #save raw values
    label.vec.save <- unique(label.vec.num)

    #  transform label into 1 and 0
    label.vec.num[label.vec.num == unique(label.vec.num)[1]] <- 1
    label.vec.num[label.vec.num != 1] <- 0

    for (i in 1:nrow(pam)){
        # remove zero-variance rows
        if (sum(pam[i,] == 1) == length(label.vec.num) |
        sum(pam[i,] == 0) == length(label.vec.num)){
        next
    }
    tmp.result <- fisher.test(pam[i,], label.vec.num)
    #print(tmp.result$p.value)
    if (tmp.result$p.value <= pvalue.cutoff  & rownames(pam)[i] != "others"){
        more.in.case <- sum(pam[i,] == 1 & label.vec.num == 1) >
        sum(pam[i,] == 1 & label.vec.num == 0)
        num.1 <- sum(pam[i,] == 1 & label.vec.num == 1)
        num.2 <- sum(pam[i,] == 1 & label.vec.num == 0)
        df.output <- rbind(df.output, c(rownames(pam)[i],
        round(as.numeric(tmp.result$p.value), 4), num.1, num.2))
    }
    }
    #return(df.output)
    if (is.null(df.output)){
        return(0)
    }
    df.output.prevalence <- percent(round((as.numeric(df.output[,3])+
    as.numeric(df.output[,4]))/ncol(pam),4))
    df.output <- cbind(df.output, df.output.prevalence)
    colnames(df.output) <- c("Name", "FDR", label.vec.save[1],
    label.vec.save[2], "prevalence")
    # FDR adjustment
    df.output[,2] <- p.adjust(df.output[,2], method = "fdr",
    n = length(df.output[,2]))
    df.output <- df.output[order(df.output[,2]),]
    return(df.output)
}
