#' Simple function to convert binary string to decimal
#'
BinToDec <- function(x) 
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))

#' Compute log2(counts per mil reads) and library size for each sample
#'
#' @param qcounts quantile normalized counts
#' @param lib.size default is colsums(qcounts)
#' @return list containing log2(quantile counts per mil reads) and library sizes
#' @export
log2CPM <- function(qcounts, lib.size=NULL){
  if (is.null(lib.size)) 
    lib.size <- colSums(qcounts)
  minimum <- min(qcounts)
  if (minimum < 0)  {
    qcounts <- qcounts-minimum
  }
  avg <- mean(qcounts)
  qcounts <- apply(qcounts,1:2,FUN=function(x){ifelse(x<=0,avg,x)})
  y <- t(log2(t(qcounts + 0.5)/(lib.size + 1) * 1e+06))
  return(list(y=y, lib.size=lib.size))
}

readPathoscopeData <- function(input_dir=".")  {
  filenames <- list.files(input_dir, pattern="*.tsv", full.names=TRUE)
  #ltbl <- lapply(filenames, read.table, skip=1, header=TRUE, sep ='\t', nrows=10)
  #lhash <- lapply(ltbl, prop_hash)
  #genomes <- c()
  #for (i in 1:length(ltbl))  {
  #  genomes <- c(genomes, levels(ltbl[i][1]))
  #}
  #lprop <- lapply(lhash, proportion, genomes)
  genomes <- c()
  for (i in 1:length(filenames))  {
    filename <- filenames[i]
    tbl <- read.table(filename, skip=1, header=TRUE, sep ='\t', nrows=10)
    hasht <- prop_hash(tbl)
    genomes <- c(genomes, levels(tbl[,1]))
  }
  genomes <- unique(genomes)
  lprop <- vector('list', length(filenames)) 
  for (i in 1:length(filenames))  {
    filename <- filenames[i]
    tbl <- read.table(filename, skip=1, header=TRUE, sep ='\t', nrows=10)
    hasht <- prop_hash(tbl)
    lprop[[i]] <- proportion(hasht, genomes) # the column data
    #names(lprop)[i] <- paste('Col', i, sep='.')  # the column name 
  }
  do.call(cbind, lprop)
  dat <- data.frame(lprop)
  rownames(dat) <- genomes
  colnames(dat) <- 1:length(filenames)
  return(dat)
}

prop_hash <- function(tbl)  {
  prop_hash <- new.env()
  for(i in 1:nrow(tbl)) {
    prop_hash[[as.character(tbl[i, 1])]] <- tbl[i, 2]
  }
  return(prop_hash)
}

proportion <- function(hasht, genomes)  {
  prop <- c()
  for(genome in genomes)  {
    if (is.null(hasht[[as.character(genome)]]))  {
      prop <- c(prop, 0)
    } else  {
      prop <- c(prop, hasht[[as.character(genome)]])
      #prop <- c(prop, mget(as.character(genome), envir=hasht))
    }
  }
  return(prop)
}