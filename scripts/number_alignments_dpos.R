suppressMessages(suppressWarnings(library(dplyr)))


num_alns = function(filename) {
    results_summary = read.csv(file = filename, header = TRUE, stringsAsFactors = FALSE)

    # extract dpos values for all aligners per each reference alignment
    split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])
    split_results = split_tibble(results_summary, 'ref_name')
    dpos = t(sapply(split_results, function(x){ return(x$dpos)}))

    col_s = ncol(dpos)
    row_s = nrow(dpos)

    results = data.frame("perfect" = integer(col_s),
                         "best" = integer(col_s),
                         "imperfect" = integer(col_s))

    # perfect alignments (dpos = 0)
    if(col_s == 2) {
        results$perfect = length(which(dpos[,2] == 0))
    } else{
        results$perfect = apply(dpos[,],2,function(x){length(which(x == 0))})
    }

    # best alignments (lower d_seq)
    for(r in 1:(row_s)) {
        for(i in which(min(dpos[r, ])==dpos[r, ])){
            results$best[i] = results$best[i] + 1
        }
    }

    # imperfect alignments (dpos != 0 when another methods has dpos = 0)
    for(r in 1:(row_s)) {
        if(min(dpos[r, ]) == 0){
            for(i in which(dpos[r, ]>0)) {
            results$imperfect[i] = results$imperfect[i] + 1
            }
        }
    }

    rownames(results) = split_results[[1]]$aligner
    return(results)
}


