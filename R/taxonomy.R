require(rentrez)

#' Find the taxonomy for the given taxon id
#'
#' @param tid Given taxon id
#' @return taxonomy LineageEx
#' @import rentrez
#' @export
findTaxonLevel <- function(tid) {
    if (is.na(tid)) 
        return(tid)
    # r_search <- entrez_search(db='taxonomy', term=tid, retmode = 'xml') 
    # r_summary <- entrez_summary(db='taxonomy', id=tid) 
    # r_fetch <- entrez_fetch(db='taxonomy', id=tid, rettype='xml',
    #   retmode='xml', parsed=TRUE)
    r_fetch <- entrez_fetch(db = "taxonomy", id = tid, rettype = "xml")
    dat <- XML::xmlToList(r_fetch)
    return(dat$Taxon$LineageEx)
}

#' Find the taxonomy for each taxon ids
#'
#' @param tids Given taxonomy ids
#' @return taxondata Data with the taxonomy information
#' @import rentrez
#' @export
findTaxonomy <- function(tids) {
    if (is.null(tids)) {
        return(NULL)
    }
    taxonLevels <- lapply(tids, FUN = findTaxonLevel)
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
    for (i in 1:length(tLineageEx)) {
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

#' Find the taxonomy level data for the given taxon level
#'
#' @param data Given Relative abundance data
#' @param taxonLevels Taxon Levels of all tids
#' @param level Selected Taxonomy Level
#' @return taxdata Taxonomy Level Data
#' @export
findTaxonLevelData <- function(data, taxonLevels, level) {
    if (level == "no rank") 
        return(data)
    taxdata <- data.frame(stringsAsFactors = FALSE)
    taxon_hash <- NULL
    length <- 0
    for (i in 1:dim(data)[1]) {
        id <- findSelectedTaxonId(taxonLevels[[i]], level)
        taxon_hash_list <- taxon_hash_update(taxon_hash, id, length)
        taxon_hash <- taxon_hash_list$taxon_hash
        index <- taxon_hash_list$index
        if (taxon_hash_list$newrow) {
            taxdata <- rbind(taxdata, data[i, ])
            rownames(taxdata)[length + 1] <- as.character(id)
        } else {
            for (j in 1:dim(data)[2]) {
                taxdata[index, j] <- taxdata[index, j] + data[i, j]
            }
        }
        length <- taxon_hash_list$length
    }
    return(taxdata)
} 
