library(PathoStat)

### Example 1: Run PathoStat
runPathoStat()

### Example 2: L1 set
#example2_data_dir <- system.file("example/data2/L1", package = "PathoStat")
#pstat <- createPathoStat(input_dir=example2_data_dir, 
#    sample_data_file="sample_data.tsv")
# Save for future use
#savePstat(pstat, outdir=".", outfileName="pstat_data_2_L1.rda")
data_dir <- system.file("data", package = "PathoStat")
# Load and run PathoStat
pstat <- loadPstat(indir=data_dir, infileName="pstat_data_2_L1.rda")
runPathoStat(pstat)

### Example 2: L2 set
#example2_data_dir <- system.file("example/data2/L2", package = "PathoStat")
#pstat <- createPathoStat(input_dir=example2_data_dir, 
#    sample_data_file="sample_data.tsv")
# Save for future use
#savePstat(pstat, outdir=".", outfileName="pstat_data_2_L2.rda")
data_dir <- system.file("data", package = "PathoStat")
# Load and run PathoStat
pstat <- loadPstat(indir=data_dir, infileName="pstat_data_2_L2.rda")
runPathoStat(pstat)

### Example 3: Asthma dataset
data_dir <- system.file("data", package = "PathoStat")
pstat <- loadPstat(indir=data_dir, infileName="asthma_pstat_data.rda")
runPathoStat(pstat)
