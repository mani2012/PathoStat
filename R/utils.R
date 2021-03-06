#' Compute log2(counts per mil reads) and library size for each sample
#'
#' @param qcounts quantile normalized counts
#' @param lib.size default is colsums(qcounts)
#' @return list containing log2(quantile counts per mil reads) and library sizes
#' @export
#' @examples
#' log2CPM(matrix(1:12, nrow = 3))

log2CPM <- function(qcounts, lib.size = NULL) {
    if (is.null(lib.size))
        lib.size <- colSums(qcounts)
    qcounts <- apply(qcounts, seq_len(2), FUN = function(x) {
        ifelse(is.null(x) || is.na(x) || is.nan(x), 0, x)
    })
    minimum <- min(qcounts)
    if (minimum < 0) {
        qcounts <- qcounts - minimum
    }
    avg <- mean(colMeans(qcounts))
    qcounts <- apply(qcounts, seq_len(2), FUN = function(x) {
        ifelse(x < 0, avg, x)
    })
    y <- t(log2(t(qcounts + 0.5)/(lib.size + 1) * 1e+06))
    return(list(y = y, lib.size = lib.size))
}

#' Reads the data from PathoScope reports and returns a list of
#' final guess relative abundance and count data
#'
#' @param input_dir Directory where the tsv files from PathoScope are located
#' @param pathoreport_file_suffix PathoScope report files suffix
#' @param use.input.files whether input dir to pathoscope files
#' or directly pathoscope files
#' @param input.files.path.vec vector of pathoscope file paths
#' @param input.files.name.vec vector of pathoscope file names
#' @return List of final guess relative abundance and count data
#' @importFrom utils read.table
#' @export
#' @examples
#' example_data_dir <- system.file("example/data", package = "PathoStat")
#' pathoreport_file_suffix <- "-sam-report.tsv"
#' datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix,
#' input.files.name.vec = as.character(1:6))

readPathoscopeData <-
    function(input_dir=".", pathoreport_file_suffix="-sam-report.tsv",
    use.input.files = FALSE, input.files.path.vec = NULL,
    input.files.name.vec = NULL) {
    if (use.input.files == FALSE){
        if (input_dir == ".") {
            input_dir <- getwd()
        }
        pattern <- paste("*", pathoreport_file_suffix, sep="")
        filenames <- list.files(input_dir, pattern = pattern, full.names = TRUE)


    } else{
        filenames <- input.files.path.vec
    }
    ltbl <- lapply(filenames, read.table, skip=1, header=TRUE, sep="\t", stringsAsFactors=TRUE)

    lgenomes <- lapply(ltbl, function(tbl) {return(levels(tbl[,1]))})
    genomes <- unique(unlist(lgenomes))
    genomes <- c(genomes, "others")
    lfl <- lapply(filenames, readLines, n = 1)
    lnumReads <- unlist(lapply(lfl, function(fl)
    {return(as.numeric(strsplit(fl, "\t")[[1]][2]))}))
    samplenames <- unlist(lapply(input.files.name.vec, function(x) {
        return(strsplit(x, pathoreport_file_suffix)[[1]])
    }))
    
    dat <- matrix(0L, nrow = length(genomes), ncol = length(samplenames))
    countdat <- matrix(0L, nrow = length(genomes), ncol = length(samplenames))    
    
    for (i in seq(length(samplenames))){
      index.tmp <- match(ltbl[[i]][,1], genomes)
      dat[index.tmp,i] <- ltbl[[i]][,3]
      countdat[index.tmp,i] <- ltbl[[i]][,4]
      read.sum <- sum(ltbl[[i]][,4])
      read.num.other <- lnumReads[i] - read.sum
      dat[nrow(dat),i] <- 1 - sum(ltbl[[i]][,3])
      countdat[nrow(dat),i] <- read.num.other
    }
    
    dat <- data.frame(dat)
    rownames(dat) <- genomes
    colnames(dat) <- samplenames
    countdat <- data.frame(countdat)
    rownames(countdat) <- genomes
    colnames(countdat) <- samplenames
    return(list(data = dat, countdata = countdat))
}

#' Loads all data from a set of PathoID reports.
#' For each column in the PathoID report, construct a
#' matrix where the rows are genomes and the columns are samples.
#' Returns a list where each element is named according to
#' the PathoID column. For example, ret[["Final.Best.Hit.Read.Numbers"]]
#' on the result of this function will get you the final count matrix.
#' Also includes elements "total_reads" and "total_genomes" from the
#' first line of the PathoID report.
#'
#' @param reportfiles Paths to report files
#' @param nrows Option to read first N rows of PathoScope reports
#' @return Returns a list where each element is named according to
#' the PathoID column. For example, ret[["Final.Best.Hit.Read.Numbers"]]
#' on the result of this function will get you the final count matrix.
#' Also includes elements "total_reads" and "total_genomes" from the
#' first line of the PathoID report.
#' @export
#' @examples
#' input_dir <- system.file("example/data", package = "PathoStat")
#' reportfiles <- list.files(input_dir, pattern = "*-sam-report.tsv",
#' full.names = TRUE)

loadPathoscopeReports <- function(reportfiles, nrows=NULL) {
    # Report basenames
    report_base <- gsub('-sam-report', '', gsub('\\.tsv$', '',
    basename(reportfiles)))

        # Get total_reads and total_genomes from first line of report
    totals <- data.frame(row.names=report_base, do.call(rbind,
    lapply(reportfiles, function(rf){read.table(rf, sep='\t', nrows=1)})
    ))
    totals <- totals[,c('V2','V4')]
    colnames(totals) <- c('total_reads','total_genomes')

    # Read all reports into list
    mlist <- lapply(reportfiles, function(rf){
    read.table(rf, sep='\t', header=TRUE, stringsAsFactors=FALSE,
    row.names=1, skip=1, comment.char="")
    })

    if(!is.null(nrows)) mlist <- lapply(mlist, 
    function(tbl){tbl[seq_len(nrows), ]})

    # Get column names for variables in reports
    vals <- colnames(mlist[[1]])
    # Get unique genomes across all reports
    genomes <- unique(do.call(c,lapply(mlist,function(tbl){rownames(tbl)})))

    # Set unobserved genome to 0
    mlist.allgenomes <- lapply(mlist, function(tbl){
    z <- tbl[genomes,]
    rownames(z) <- genomes
    z[is.na(z)] <- 0
    z
    })
    # Sanity check: colnames and rownames should be the same for all tables
    stopifnot(all(vapply(mlist.allgenomes,
    function(tbl){all(rownames(tbl)==genomes)})))
    stopifnot(all(vapply(mlist.allgenomes,
    function(tbl){all(colnames(tbl)==vals)})))

    ret <- lapply(vals, function(v) {
    z <- data.frame(row.names=genomes, do.call(cbind, lapply(
    mlist.allgenomes,function(tbl){ tbl[,v]})))
    colnames(z) <- report_base
    z
    })
    names(ret) <- vals
    ret[['total_reads']] <- totals[report_base, 'total_reads']
    names(ret[['total_reads']]) <- report_base
    ret[['total_genomes']] <- totals[report_base, 'total_genomes']
    names(ret[['total_genomes']]) <- report_base

    ret
}

prop_hash <- function(tbl) {
    prop_hash <- new.env()
    sum <- 0
    for (i in seq_len(nrow(tbl))) {
        prop_hash[[as.character(tbl[i, 1])]] <- tbl[i, 2]
        sum <- sum + tbl[i, 2]
    }
    if (sum < 1) {
        prop_hash[["others"]] <- 1 - sum
    }
    return(prop_hash)
}

proportion <- function(hasht, genomes, countdata = FALSE, numReads = 1) {
    prop <- c()
    for (genome in genomes) {
        if (is.null(hasht[[as.character(genome)]])) {
            prop <- c(prop, 0)
        } else {
            if (countdata) {
                prop <- c(prop, round(numReads * hasht[[as.character(genome)]]))
            } else {
                prop <- c(prop, hasht[[as.character(genome)]])
            }
        }
    }
    return(prop)
}


proportionc <- function(lhasht, genomes, lnumReads, i) {
    propc <- proportion(lhasht[[i]], genomes, countdata = TRUE,
        numReads=lnumReads[i])
    return(propc)
}

#' Greps the tid from the given identifier string
#'
#' @param id Given identifier string
#' @return tid string
#' @export
#' @examples
#' grepTid("ti|700015|org|Coriobacterium_glomerans_PW2")

grepTid <- function(id) {
    tid <- strsplit(id, "\\|")
    tid <- sapply(tid, function(x) x[2])
    return(tid)
}


#' Save the pathostat object to R data(.rda) file
#'
#' @param pstat pathostat object
#' @param outdir Output Directory of the .rda file
#' @param outfileName File name of the .rda file
#' @return outfile .rda file
#' @export
#' @examples
#' data_dir_test <- system.file("data", package = "PathoStat")
#' pstat_test <- loadPstat(indir=data_dir_test,
#' infileName="pstat_data_2_L1.rda")
#' outfile <- savePstat(pstat_test)
savePstat <- function(pstat, outdir=".", outfileName="pstat_data.rda") {
    if (outdir == ".") {
        outdir <- getwd()
    }
    outfile <- file.path(outdir, outfileName)
    save(pstat, file=outfile)
    return(outfile)
}

#' Load the R data(.rda) file with pathostat object
#'
#' @param indir Input Directory of the .rda file
#' @param infileName File name of the .rda file
#' @return pstat pathostat object (NULL if it does not exist)
#' @export
#' @examples
#' data_dir <- system.file("data", package = "PathoStat")
#' infileName <- "pstat_data.rda"
#' pstat <- loadPstat(data_dir, infileName)

loadPstat <- function(indir=".", infileName="pstat_data.rda") {
    if (indir == ".") {
        indir <- getwd()
    }
    infile <- file.path(indir, infileName)
    pstat <- NULL
    load(file=infile)
    return(pstat)
}

#' Return the Relative Abundance (RA) data for the given count OTU table
#'
#' @param count_otu Count OTU table
#' @return ra_otu Relative Abundance (RA) OTU table
#' @export
#' @examples
#' data_dir <- system.file("data", package = "PathoStat")
#' infileName <- "pstat_data.rda"
#' pstat_test <- loadPstat(data_dir, infileName)
#' ra_otu <- findRAfromCount(phyloseq::otu_table(pstat_test))

findRAfromCount <- function(count_otu) {
    ra_otu <- otu_table(count_otu)
    numcol <- dim(count_otu)[2]
    for (i in seq_len(numcol))  {
        ra_otu[,i] <- ra_otu[,i]/sum(ra_otu[,i])
    }
    return(ra_otu)
}

#' Format taxonomy table for rendering
#'
#' @param ttable Taxonomy table
#' @import DT tidyr
#' @return Formatted table suitable for rendering with. DT::renderDataTable
formatTaxTable <- function(ttable) {
    ttable
}


#' Summarize sample
#'
#' Creates a table of summary metrics
#'
#' @param pstat Input pstat
#'
#' @return A data.frame object of summary metrics.
#' @export summarizeTable
#' @examples
#' data_dir_test <- system.file("data", package = "PathoStat")
#' pstat_test <- loadPstat(indir=data_dir_test,
#' infileName="pstat_data_2_L1.rda")
#' st.tmp <- summarizeTable(pstat_test)

summarizeTable <- function(pstat){
    return(data.frame("Metric" = c("Number of Samples",
    "Number of taxon elements",
    "Average number of reads per sample",
    "Average number of taxon elements per sample",
    "Samples with <5 taxon elements",
    "Taxon elements with reads across all samples"),
    "Value" = c(ncol(pstat@otu_table@.Data),
    nrow(pstat@otu_table@.Data),
    as.integer(mean(colSums(pstat@otu_table@.Data))),
    as.integer(mean(colSums(pstat@otu_table@.Data > 0))),
    sum(colSums(pstat@otu_table@.Data != 0) < 5),
    sum(rowSums(pstat@otu_table@.Data) == 0))))
}





#' Compute percentage
#'
#' @param x a number or a vector
#' @param digits how many digit of percentage
#' @param format numeric format, "f" for float
#' @return the percentage
#' @export
#' @examples
#' pecent.vec <- percent(c(0.9, 0.98))
percent <- function(x, digits = 2, format = "f") {
    paste0(formatC(100 * x, format = format, digits = digits), "%")
}
