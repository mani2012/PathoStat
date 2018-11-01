
#' Find the taxonomy for the given taxon id name
#'
#' @param pstat pathostat object
#' @param input.id.vec names containing id
#' @param tax.level target taxon level
#' @return target taxon level names
#' @export
#' @examples
#' data_dir_test <- system.file("data", package = "PathoStat")
#' pstat_test <- loadPstat(indir=data_dir_test,
#' infileName="pstat_data_2_L1.rda")
#' names.new <- TranslateIdToTaxLevel(pstat_test,
#' c("ti|862962|org|Bacteroides_fragilis_638R",
#' "ti|697329|org|Ruminococcus_albus_7" ),
#' "genus")
TranslateIdToTaxLevel <- function(pstat, input.id.vec, tax.level){
    tax.df <- pstat@tax_table@.Data
    name.out <- as.character(tax.df[match(input.id.vec, rownames(tax.df)),
                                    which(colnames(tax.df) == tax.level)])
    return(name.out)
}


#' Find the taxonomy for maximum 300 tids
#'
#' @param tids Given taxonomy ids
#' @return taxondata Data with the taxonomy information
#' @import rentrez
#' @export
#' @examples
#' example_data_dir <- system.file("example/data", package = "PathoStat")
#' pathoreport_file_suffix <- "-sam-report.tsv"
#' datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix)
#' dat <- datlist$data
#' ids <- rownames(dat)
#' tids <- unlist(lapply(ids, FUN = grepTid))
#' taxonLevels <- findTaxonomy300(tids[1:5])

findTaxonomy300 <- function(tids) {
    if (is.null(tids)) {
        return(NULL)
    }

    na.vec <- c()
    for (i in seq_len(length(tids))){
        if(is.na(tids[i])){
            na.vec <- c(na.vec, i)
        }
    }

    r_fetch <- entrez_fetch(db = "taxonomy", id = tids, rettype = "xml")
    dat <- XML::xmlToList(r_fetch)
    taxonLevels <- lapply(dat, function(x) x$LineageEx)
    if(!is.null(na.vec)){
        for(i in seq_len(length(na.vec))){
        taxonLevels <- append(taxonLevels, list(NA), na.vec[i]-1)
        }
    }


    return(taxonLevels)
}



#' Find the taxonomy for unlimited tids
#'
#' @param tids Given taxonomy ids
#' @return taxondata Data with the taxonomy information
#' @import rentrez
#' @export
#' @examples
#' example_data_dir <- system.file("example/data", package = "PathoStat")
#' pathoreport_file_suffix <- "-sam-report.tsv"
#' datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix)
#' dat <- datlist$data
#' ids <- rownames(dat)
#' tids <- unlist(lapply(ids, FUN = grepTid))
#' taxonLevels <- findTaxonomy(tids[1:5])

findTaxonomy <- function(tids) {
    if (is.null(tids)) {
        return(NULL)
    }
    if (length(tids) <= 300){
        taxonLevels <- findTaxonomy300(tids)
    } else{
        taxonLevels <- list()
        batch.num <- ceiling(length(tids)/300)
        for (i in seq_len(batch.num)){
            if (i == batch.num){
                tids.batch <- tids[((i-1)*300 + 1):length(tids)]
            }else{
                tids.batch <- tids[((i-1)*300 + 1):(i*300)]
            }
        taxonLevels <- c(taxonLevels, findTaxonomy300(tids.batch))
        print(i) 
        }
    }


    return(taxonLevels)
}


findSelectedTaxonId <- function(tLineageEx, level) {
    id <- "others"
    if (is.null(tLineageEx)) {
        return(id)
    }
    if (is.na(tLineageEx[1])) {
        return(id)
    }
    for (i in seq_len(length(tLineageEx))) {
        rank <- tLineageEx[[i]]["Rank"]
        scientificName <- tLineageEx[[i]]["ScientificName"]
        taxid <- tLineageEx[[i]]["TaxId"]
        if (!is.null(rank)) {
            if (rank == level) {
                id <- paste0("ti|", taxid, "|", rank, "|", scientificName)
                break
            }
        }
    }
    return(id)
}

findTaxonLevels <- function(data) {
    ids <- rownames(data)
    tids <- unlist(lapply(ids, FUN = grepTid))
    taxonLevels <- findTaxonomy(tids)
    return(taxonLevels)
}

taxon_hash_update <- function(taxon_hash, taxon, length) {
    if (is.null(taxon_hash)) {
        taxon_hash <- new.env()
    }
    index <- length + 1
    newrow <- FALSE
    if (is.null(taxon_hash[[as.character(taxon)]])) {
        taxon_hash[[as.character(taxon)]] <- index
        length <- length + 1
        newrow <- TRUE
    } else {
        index <- taxon_hash[[as.character(taxon)]]
    }
    return(list(taxon_hash = taxon_hash, length = length, index = index,
        newrow = newrow))
}


#' Find the Taxonomy Information Matrix
#'
#' @param names Row names of the taxonomy matrix
#' @param taxonLevels Taxon Levels of all tids
#' @return taxmat Taxonomy Information Matrix
#' @export
#' @examples
#' example_data_dir <- system.file("example/data", package = "PathoStat")
#' pathoreport_file_suffix <- "-sam-report.tsv"
#' datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix)
#' dat <- datlist$data
#' ids <- rownames(dat)
#' tids <- unlist(lapply(ids, FUN = grepTid))
#' taxonLevels <- findTaxonomy(tids[1:5])
#' taxmat <- findTaxonMat(ids[1:5], taxonLevels)

findTaxonMat <- function(names, taxonLevels) {
    # tax.name <- c('superkingdom', 'kingdom', 'phylum', 'class', 'order',
    #     'suborder', 'family', 'subfamily', 'genus', 'subgenus', 'species',
    #     'subspecies', 'no rank')
    tax.name <- c('superkingdom', 'kingdom', 'phylum', 'class', 'order',
        'family', 'genus', 'species', 'no rank')
    tl <- c()
    for (i in seq_len(length(tax.name))) {
        tl[tax.name[i]] <- "others"
    }
    taxmat <- NULL
    for (i in seq_len(length(taxonLevels))) {
        taxrow <- tl
        tLineageEx <- taxonLevels[[i]]
        for (j in seq_len(length(tLineageEx))) {
            rank <- tLineageEx[[j]]["Rank"]
            #taxid <- tLineageEx[[j]]["TaxId"]
            scientificName <- tLineageEx[[j]]["ScientificName"]
            if (!is.null(rank) && !is.na(rank) && rank %in% tax.name) {
                #taxrow[as.character(rank)] <- as.character(taxid)
                taxrow[as.character(rank)] <- as.character(scientificName)
            }
        }
        taxmat <- rbind(taxmat, taxrow)
    }

    rownames(taxmat) <- names
    #print(taxmat)
    return(taxmat)
}
