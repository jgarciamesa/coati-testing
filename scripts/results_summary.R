

source("scripts/pa_distance.R")
source("scripts/kaks.R")

main = function(aln, aligners) {
    # read info about reference alns creation
    reference_df = read.csv("data/ref_alignments.csv")
    len = length(aligners)

    ############################################################################
    # calculate dseq for inferred alignments
    dseqs = pa_distance_main(aln, aligners, 1)  # 1: dseq
    
    ############################################################################
    # calculate selection for reference and inferred alignments 
    selection = vector(mode = "double", length = len)

    # omega for reference aln
    ref = process_aln(paste0("data/ref_alignments/"), aln)
    ref_omega = sapply(ref,omega)
    
    # omega for inferred alns
    models = sapply(aligners, function(x) {process_aln(paste0("aln/ref/",x), aln)})
    models_omega = sapply(models,omega)
    
    # group distance (dseq) and selection (omega) results
    distance_selection_df = data.frame(ref_name = rep(aln, len),
                            aligner = aligners,
                            dseq = dseqs,
                            ref_omega = rep(ref_omega, len),
                            aln_omega = models_omega)

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