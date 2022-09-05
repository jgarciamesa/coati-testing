
if(!interactive()) {
    ARGS = commandArgs(trailingOnly = TRUE)
    data = read.csv(ARGS[1])
    summary = apply(data[,-1],2,sum)/nrow(data)
    write.table(file=ARGS[2], x = summary, quote = FALSE, sep = "\t\t", col.names = FALSE)
}