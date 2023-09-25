# Filter out CDS files where either sequence is longer than a maximum length,
#   or the reference (gorilla) contains an early stop or an incomplete codon.
# Return the CDS files that pass the filtering with the gorilla sequence first.

library(seqinr)

filter_seqs <- function(N, max_len, outfile) {
    CDS <- list.files("raw_data/", pattern = "*.fasta", full.names = TRUE)[1:N]
    filtered = c()

    for(file in CDS) {
        fasta <- read.alignment(file = file, format = "fasta", forceToLower = FALSE)
        fasta$seq <- toupper(fasta$seq) # 'forceToLower' on read.alignment is broken

        gor_len <- nchar(fasta$seq[2])
        hum_len <- nchar(fasta$seq[1])

        # check for length
        if(gor_len + hum_len > max_len) {
            next
        }

        # check for incomplete codons in the gorilla sequence
        if(gor_len %% 3 != 0) {
            next
        }

        # check for early stop codons in the gorilla sequence
        codons = unlist(strsplit(fasta$seq[2], "(?<=.{3})", perl = TRUE))
        if(any(c("TAA", "TAG", "TGA") %in% codons[-length(codons)])) {
            next
        }

        # check for ambiguous nucs on gorilla sequence
        gor_nucs = sort(unique(unlist(strsplit(fasta$seq[2], split = ""))))
        if(length(gor_nucs) != 4) {
            next
        } else if(!all.equal(gor_nucs, c("A", "C", "G", "T"))) {
            next
        }

        # save name of CDS
        filtered = c(filtered, basename(file))
        
        # rewrite CDS with gorilla as first sequence
        fasta$seq = c(fasta$seq[2], fasta$seq[1])
        fasta$nam = c(fasta$nam[2], fasta$nam[1])
        write.fasta(sequences = as.list(fasta$seq),
                    as.string = TRUE, 
                    names = fasta$nam,
                    file.out = paste0("data/cds/", basename(file)))
    }
    # save list of filtered CDSs
    write.table(x = filtered,
                file = outfile,
                col.names = FALSE,
                row.names = FALSE,
                quote = FALSE,
                sep = ",")
}

if(!interactive()) {
    library(seqinr)

    ARGS = commandArgs(trailingOnly = TRUE)
    stopifnot(length(ARGS) == 3)
    filter_seqs(ARGS[1], ARGS[2], ARGS[3])
}