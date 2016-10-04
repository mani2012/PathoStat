#' Compute log2(counts per mil reads) and library size for each sample
#'
#' @param qcounts quantile normalized counts
#' @param lib.size default is colsums(qcounts)
#' @return list containing log2(quantile counts per mil reads) and library sizes
#' @export
#' @examples
#' example_data_dir <- system.file("example/data", package = "PathoStat")
#' pathoreport_file_suffix <- "-sam-report.tsv"
#' datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix)
#' countdat <- datlist$countdata
#' lcpm <- log2CPM(countdat)
log2CPM <- function(qcounts, lib.size = NULL) {
    if (is.null(lib.size)) 
        lib.size <- colSums(qcounts)
    minimum <- min(qcounts)
    if (minimum < 0) {
        qcounts <- qcounts - minimum
    }
    avg <- mean(colMeans(qcounts))
    qcounts <- apply(qcounts, 1:2, FUN = function(x) {
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
#' @return List of final guess relative abundance and count data
#' @importFrom utils read.table
#' @export
#' @examples
#' example_data_dir <- system.file("example/data", package = "PathoStat")
#' readPathoscopeData(input_dir=example_data_dir)
readPathoscopeData <- 
    function(input_dir=".", pathoreport_file_suffix="-sam-report.tsv") {
    if (input_dir == ".") {
        input_dir <- getwd()
    }
    pattern <- paste("*", pathoreport_file_suffix, sep="")
    filenames <- list.files(input_dir, pattern = pattern, full.names = TRUE)

    ltbl <- lapply(filenames, read.table, skip=1, header=TRUE, sep="\t", 
        nrows=10)
    lgenomes <- lapply(ltbl, function(tbl) {return(levels(tbl[,1]))})
    genomes <- unique(unlist(lgenomes))
    genomes <- c(genomes, "others")
    
    lfl <- lapply(filenames, readLines, n = 1)
    lnumReads <- unlist(lapply(lfl, function(fl) 
        {return(as.numeric(strsplit(fl, "\t")[[1]][2]))}))
    lhasht <- lapply(ltbl, prop_hash)
    lprop <- lapply(lhasht, proportion, genomes)
    lcprop <- lapply(seq_along(lhasht), proportionc, lhasht=lhasht, 
        genomes=genomes, lnumReads=lnumReads)
    
    do.call(cbind, lprop)
    do.call(cbind, lcprop)
    samplenames <- unlist(lapply(filenames, function(x) {
        return(strsplit(basename(x), pathoreport_file_suffix)[[1]])
    }))
    dat <- data.frame(lprop)
    rownames(dat) <- genomes
    colnames(dat) <- samplenames
    countdat <- data.frame(lcprop)
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
#'     full.names = TRUE)
#' loadPathoscopeReports(reportfiles)
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
    
    if(!is.null(nrows)) mlist <- lapply(mlist, function(tbl){tbl[1:nrows, ]})

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
    stopifnot(all(sapply(mlist.allgenomes,
        function(tbl){all(rownames(tbl)==genomes)})))
    stopifnot(all(sapply(mlist.allgenomes,
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
#' tid <- grepTid("ti|367928|org|Bifidobacterium_adolescentis_ATCC_15703")
grepTid <- function(id) {
    tid <- unlist(strsplit(id, ".org"))[1]
    tid <- unlist(strsplit(tid, "ti."))[2]
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
#' data(pstat_data)
#' outfile <- savePstat(pstat)
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
#' pstat <- loadPstat(data_dir, infileName)
#' ra_otu <- findRAfromCount(phyloseq::otu_table(pstat))
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
#' @return Formatted table suitable for rendering with. DT::renderDataTable
formatTaxTable <- function(ttable) {
    ttable
}

