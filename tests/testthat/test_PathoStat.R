library(PathoStat)
context("PathoStat functionality")

test_that("createPathoStat", {
    example_data_dir <- system.file("example/data", package = "PathoStat")
    pstat <- createPathoStat(input_dir=example_data_dir, 
        sample_data_file="sample_data.tsv")
    expect_equal(dim(otu_table(pstat)), c(31, 9))
    expect_equal(dim(sample_data(pstat)), c(9, 18))
    expect_equal(dim(tax_table(pstat)), c(31, 9))
})

test_that("runPathoStat", {
    outputfile <- runPathoStat(interactive = FALSE)
    expect_equal(basename(outputfile), "PathoStat_report.html")
})

test_that("log2CPM", {
    example_data_dir <- system.file("example/data", package = "PathoStat")
    pathoreport_file_suffix <- "-sam-report.tsv"
    datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix)
    countdat <- datlist$countdata
    lcpm <- log2CPM(countdat)
    expect_equal(dim(lcpm$y), dim(countdat))
    expect_equal(length(lcpm$lib.size), dim(countdat)[2])
})

test_that("readPathoscopeData", {
    example_data_dir <- system.file("example/data", package = "PathoStat")
    pdata <- readPathoscopeData(input_dir=example_data_dir)
    expect_equal(dim(pdata$data), c(31, 9))
    expect_equal(dim(pdata$countdata), c(31, 9))
})

test_that("loadPathoscopeReports", {
    input_dir <- system.file("example/data", package = "PathoStat")
    reportfiles <- list.files(input_dir, pattern = "*-sam-report.tsv", 
        full.names = TRUE)
    ret <- loadPathoscopeReports(reportfiles)
    expect_equal(length(ret$total_reads),9)
    expect_equal(length(ret$total_genomes),9)
})

test_that("grepTid", {
    tid <- grepTid("ti|367928|org|Bifidobacterium_adolescentis_ATCC_15703")
    expect_equal(tid, "367928")
})

test_that("savePstat", {
    data(pstat_data)
    outfile <- savePstat(pstat)
    expect_equal(basename(outfile), "pstat_data.rda")
})

test_that("loadPstat", {
    data_dir <- system.file("data", package = "PathoStat")
    infileName <- "pstat_data.rda"
    pstat <- loadPstat(data_dir, infileName)
    expect_equal(dim(otu_table(pstat)), c(41, 33))
    expect_equal(dim(sample_data(pstat)), c(33, 18))
    expect_equal(dim(tax_table(pstat)), c(41, 9))
})

test_that("findRAfromCount", {
    data_dir <- system.file("data", package = "PathoStat")
    infileName <- "pstat_data.rda"
    pstat <- loadPstat(data_dir, infileName)
    ra_otu <- findRAfromCount(phyloseq::otu_table(pstat))
    expect_equal(dim(ra_otu), c(41, 33))
})

test_that("findTaxonLevel", {
    example_data_dir <- system.file("example/data", package = "PathoStat")
    pathoreport_file_suffix <- "-sam-report.tsv"
    datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix)
    dat <- datlist$data
    ids <- rownames(dat)
    tids <- unlist(lapply(ids, FUN = grepTid))
    taxonLevel <- findTaxonomy(tids[1])
    expect_equal(taxonLevel[[1]]$Taxon$TaxId, "131567")
})

test_that("findTaxonomy", {
    example_data_dir <- system.file("example/data", package = "PathoStat")
    pathoreport_file_suffix <- "-sam-report.tsv"
    datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix)
    dat <- datlist$data
    ids <- rownames(dat)
    tids <- unlist(lapply(ids, FUN = grepTid))
    taxonLevels <- findTaxonomy(tids[1:5])
    expect_equal(length(taxonLevels), 5)
})

test_that("findTaxonMat", {
    example_data_dir <- system.file("example/data", package = "PathoStat")
    pathoreport_file_suffix <- "-sam-report.tsv"
    datlist <- readPathoscopeData(example_data_dir, pathoreport_file_suffix)
    dat <- datlist$data
    ids <- rownames(dat)
    tids <- unlist(lapply(ids, FUN = grepTid))
    taxonLevels <- findTaxonomy(tids[1:5])
    taxmat <- findTaxonMat(ids[1:5], taxonLevels)
    expect_equal(dim(taxmat), c(5,9))
})
