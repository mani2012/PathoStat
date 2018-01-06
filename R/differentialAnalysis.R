#' Mann-whitney test for a dataframe
#'
#' @param df Input data object that contains the data to be tested. Required
#' @param label.vec.num The target binary condition. Required
#' @param pvalue.cutoff choose p-value cut-off
#' @return df.output object
#' @export
#' @examples
#' Wilcox_Test_df(df.list[[i]], label.vec.num, pvalue.cutoff)

Wilcox_Test_df <- function(df, label.vec.num, pvalue.cutoff = 0.05) {
  df.output <- NULL
  #save raw values
  label.vec.save <- unique(label.vec.num)
  
  # transform label into 1 and 0
  label.vec.num[label.vec.num == unique(label.vec.num)[1]] <- 1
  label.vec.num[label.vec.num != 1] <- 0
  
  for (i in 1:nrow(df)){
    # remove zero-variance rows
    if (sum(df[i,] == 1) == length(label.vec.num) | sum(df[i,] == 0) == length(label.vec.num)){
      next
    }
    tmp.result <- wilcox.test(df[i,which(label.vec.num == 1)], df[i,which(label.vec.num == 0)], correct=FALSE)
    if (tmp.result$p.value <= pvalue.cutoff){
      num.1 <- sum(df[i,which(label.vec.num == 1)])
      num.2 <- sum(df[i,which(label.vec.num == 0)])
      df.output <- rbind(df.output, c(rownames(df)[i], round(as.numeric(tmp.result$p.value), 4), num.1, num.2))
    }
  }
  if (is.null(df.output)){
    return(0)
  }
  colnames(df.output) <- c("Name", "P-value", label.vec.save[1],label.vec.save[2])
  return(df.output)
}






#' transform cpm counts to presence-absence matrix
#'
#' @param df Input data object that contains the data to be tested. Required
#' @return df.output object
#' @export
#' @examples
#' GET_PAM(df)

GET_PAM <- function(df) {
  for (i in 1:nrow(df)){
    df[i,] <- as.numeric(df[i,] > 0)
  }
  return(df)
}



#' Given PAM and disease/control annotation, do Chi-square test for each row of PAM
#'
#' @param pam Input data object that contains the data to be tested. Required
#' @param label.vec.num The target binary condition. Required
#' @param pvalue.cutoff choose p-value cut-off
#' @return df.output object
#' @export
#' @examples
#' Chisq_Test_Pam(pam, label.vec.num, pvalue.cutoff)

Chisq_Test_Pam <- function(pam, label.vec.num, pvalue.cutoff = 0.05) {
  df.output <- NULL
  
  #save raw values
  label.vec.save <- unique(label.vec.num)
  
  # transform label into 1 and 0
  label.vec.num[label.vec.num == unique(label.vec.num)[1]] <- 1
  label.vec.num[label.vec.num != 1] <- 0
  
  
  for (i in 1:nrow(pam)){
    # remove zero-variance rows
    if (sum(pam[i,] == 1) == length(label.vec.num) | sum(pam[i,] == 0) == length(label.vec.num)){
      next
    }
    tmp.result <- chisq.test(pam[i,], label.vec.num, correct=FALSE)
    if (tmp.result$p.value <= pvalue.cutoff){
      num.1 <- sum(pam[i,] == 1 & label.vec.num == 1)
      num.2 <- sum(pam[i,] == 1 & label.vec.num == 0)
      df.output <- rbind(df.output, c(rownames(pam)[i], round(as.numeric(tmp.result$p.value), 4), num.1, num.2))
    }
  }
  if (is.null(df.output)){
    return(0)
  }
  colnames(df.output) <- c("Name", "P-value", label.vec.save[1], label.vec.save[2])
  return(df.output)
}



#' Given PAM and disease/control annotation, do Chi-square test for each row of PAM
#'
#' @param pam Input data object that contains the data to be tested. Required
#' @param label.vec.num The target binary condition. Required
#' @param pvalue.cutoff choose p-value cut-off
#' @return df.output object
#' @export
#' @examples
#' Fisher_Test_Pam(pam, label.vec.num, pvalue.cutoff)

Fisher_Test_Pam <- function(pam, label.vec.num, pvalue.cutoff = 0.05) {
  df.output <- NULL
  
  #save raw values
  label.vec.save <- unique(label.vec.num)
  
  # transform label into 1 and 0
  label.vec.num[label.vec.num == unique(label.vec.num)[1]] <- 1
  label.vec.num[label.vec.num != 1] <- 0
  
  for (i in 1:nrow(pam)){
    # remove zero-variance rows
    if (sum(pam[i,] == 1) == length(label.vec.num) | sum(pam[i,] == 0) == length(label.vec.num)){
      next
    }
    tmp.result <- fisher.test(pam[i,], label.vec.num)
    #print(tmp.result$p.value)
    if (tmp.result$p.value <= pvalue.cutoff){
      more.in.case <- sum(pam[i,] == 1 & label.vec.num == 1) > sum(pam[i,] == 1 & label.vec.num == 0)
      num.1 <- sum(pam[i,] == 1 & label.vec.num == 1)
      num.2 <- sum(pam[i,] == 1 & label.vec.num == 0)
      df.output <- rbind(df.output, c(rownames(pam)[i], round(as.numeric(tmp.result$p.value), 4), num.1, num.2))
    }
  }
  #return(df.output)
  if (is.null(df.output)){
    return(0)
  }
  colnames(df.output) <- c("Name", "P-value", label.vec.save[1], label.vec.save[2])
  return(df.output)
}
