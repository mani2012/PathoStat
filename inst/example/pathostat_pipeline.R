library(PathoStat)

nbatch <- 3
ncond <- 2
npercond <- 10

### apply PathoStat
batch <- rep(1:nbatch, each=ncond*npercond)
condition <- rep(rep(1:ncond, each=npercond), nbatch)
pdata <- data.frame(batch, condition)
modmatrix = model.matrix(~as.factor(condition), data=pdata)

example_data_dir <- system.file("example/data", package = "PathoStat")
pathoStat(input_dir=example_data_dir, batch=batch, condition=condition, 
        report_file="pathostat_report.html", report_dir=".", report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE)
