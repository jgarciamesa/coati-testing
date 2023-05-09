

source("scripts/distance.R")
source("scripts/kaks.R")

main = function(aln, aligners) {
    # read info about reference alns creation
    reference_df = read.csv("data/ref_alignments.csv")
    len = length(aligners)

    ############################################################################
    # calculate dseq & dpos for inferred alignments
    dseqs = distance(aln, aligners, "dseq")
    dposs = distance(aln, aligners, "dpos")
    
    ############################################################################
    # omega for reference aln
    ref = process_aln(paste0("data/ref_alignments/"), aln)
    ref_omega = sapply(ref,omega)
    
    # omega for inferred alns
    models = sapply(aligners, function(x) {process_aln(paste0("aln/ref/",x), aln)})
    models_omega = sapply(models,omega)
    
    # score alignments
    scores = c()
    for(aligner in aligners) {
        s = suppressWarnings(system(paste0("bin/coati-alignpair -s aln/ref/", aligner, "/", aln), intern = TRUE, ignore.stderr = TRUE))
        scores = c(scores, as.double(s))
    }

    if(length(scores) != len) {
        scores = rep(NA, len)
    }

    # group distance (dseq) and selection (omega) results
    distance_selection_df = data.frame(ref_name = rep(aln, len),
                            aligner = aligners,
                            dseq = dseqs,
                            dpos = dposs,
                            ref_omega = rep(ref_omega, len),
                            aln_omega = models_omega,
                            score = scores)

    # merge results metrics with gap origin information
    results = merge(reference_df, distance_selection_df, by = "ref_name")
    write.table(results,
              file = "results/results_summary.csv",
              append = TRUE,
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE,
              sep = ",")
}

if(!interactive()) {
    ARGS = commandArgs(trailingOnly = TRUE)
    main(ARGS[1], ARGS[-1])
}