library(seqinr)
library(stringr)
library(purrr)

ARGS <- commandArgs(trailingOnly = TRUE)

coati_bin <- ARGS[[1]]
coati_model <- ARGS[[2]]
input_path <- ARGS[[3]]
output_path <- ARGS[[4]]

# Flip the order of the sequences
input_data <- seqinr::read.fasta(input_path, seqtype = "DNA",
    forceDNAtolower = FALSE, as.string = TRUE)
input_data <- rev(input_data)

rev_file <- tempfile(fileext=".fasta")
seqinr::write.fasta(input_data, names(input_data), file.out = rev_file, as.string = TRUE)

results <- system2(coati_bin,
    args = c("-m", coati_model,
        "-o", "fasta:-",
        rev_file),
    stdout = TRUE
)

if(is.integer(results)) {
    quit(save = "no", status = results)
} else if(is.integer(attr(results, "status"))) {
    quit(save = "no", status = attr(results, "status"))
}

rev_aligned <- seqinr::read.fasta(textConnection(results), seqtype = "DNA",
    forceDNAtolower = FALSE, as.string = TRUE)
rev_aligned <- rev(rev_aligned)

seqinr::write.fasta(rev_aligned , names(rev_aligned), file.out = output_path, as.string = TRUE)
