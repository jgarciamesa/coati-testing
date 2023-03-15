suppressMessages(suppressWarnings(library(dplyr)))


num_alns = function(filename) {
    results_summary = read.csv(file = filename, header = TRUE, stringsAsFactors = FALSE)

    # extract dseq values for all aligners per each reference alignment
    split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])
    split_results = split_tibble(results_summary, 'ref_name')
    dseq = t(sapply(split_results, function(x){ return(x$dseq)}))

    col_s = ncol(dseq)
    row_s = nrow(dseq)

    results = data.frame("perfect" = integer(col_s),
                         "best" = integer(col_s),
                         "imperfect" = integer(col_s))

    # perfect alignments (dseq = 0)
    if(col_s == 2) {
        results$perfect = length(which(dseq[,2] == 0))
    } else{
        results$perfect = apply(dseq[,],2,function(x){length(which(x == 0))})
    }

    # best alignments (lower d_seq)
    for(r in 1:(row_s)) {
        for(i in which(min(dseq[r, ])==dseq[r, ])){
            results$best[i] = results$best[i] + 1
        }
    }

    # imperfect alignments (dseq != 0 when another methods has dseq = 0)
    for(r in 1:(row_s)) {
        if(min(dseq[r, ]) == 0){
            for(i in which(dseq[r, ]>0)) {
            results$imperfect[i] = results$imperfect[i] + 1
            }
        }
    }

    rownames(results) = split_results[[1]]$aligner
    return(results)
}


