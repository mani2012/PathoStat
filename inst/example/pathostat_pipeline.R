library(PathoStat)

### Example batch and condition
nbatch <- 11
ncond <- 3
npercond <- 11
#batch <- rep(1:nbatch, each=ncond*npercond/nbatch)
subject_id <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 
    10, 10, 10, 12, 12, 12, 13, 13, 13, 15, 15, 15)
batch <- unlist(lapply(subject_id, FUN=function(id) {paste("Person", id)}))
#condition <- rep(1:ncond, each=npercond)
diet <- c(1, 3, 2, 3, 1, 2, 2, 3, 1, 3, 2, 1, 3, 2, 1, 3, 1, 2, 1, 2, 3, 
    2, 1, 3, 3, 1, 2, 3, 1, 2, 2, 3, 1)
diet_key <- c("simple", "refined", "unrefined")
condition <- diet_key[diet]
#pdata <- data.frame(batch, condition)
#modmatrix = model.matrix(~as.factor(condition), data=pdata)

### Example PathoScope report files directory
example_data_dir <- system.file("example/data", package = "PathoStat")
### apply PathoStat
pathoStat(input_dir=example_data_dir, batch=batch, condition=condition, 
    report_file="pathostat_report.html", report_dir=".", report_option_binary=
    "111111111", view_report=FALSE, interactive=TRUE)

### Example 2 PathoScope report files directory
example2_data_dir <- system.file("example/data2", package = "PathoStat")
### apply PathoStat
pathoStat(input_dir=example2_data_dir, batch=batch, condition=condition, 
    report_file="pathostat_report.html", report_dir=".", report_option_binary=
    "111111111", view_report=FALSE, interactive=TRUE)

example2_data_dir <- system.file("example/data2/L1", package = "PathoStat")
### apply PathoStat
pathoStat(input_dir=example2_data_dir, batch=batch, condition=condition, 
    report_file="pathostat_report.html", report_dir=".", report_option_binary=
    "111111111", view_report=FALSE, interactive=TRUE)

example2_data_dir <- system.file("example/data2/L2", package = "PathoStat")
### apply PathoStat
pathoStat(input_dir=example2_data_dir, batch=batch, condition=condition, 
    report_file="pathostat_report.html", report_dir=".", report_option_binary=
    "111111111", view_report=FALSE, interactive=TRUE)
