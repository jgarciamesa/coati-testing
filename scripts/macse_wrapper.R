library(seqinr)
library(stringr)
library(purrr)

# ENSG00000010278

ARGS <- commandArgs(trailingOnly = TRUE)

macse_jar <- ARGS[[1]]
input_path <- ARGS[[2]]
output_path <- ARGS[[3]]

# Split input file into reliable and less reliable input files
input_data <- seqinr::read.fasta(input_path, seqtype = "DNA", forceDNAtolower = FALSE)
seq_data <- input_data[1]
seq_lr_data <- input_data[-1]

seq_file <- tempfile(fileext=".fasta")
seq_lr_file <- tempfile(fileext=".fasta")
out_aa_file <- tempfile(fileext=".fasta")
out_nt_file <- tempfile(fileext=".fasta")

seqinr::write.fasta(seq_data, names(seq_data), file.out = seq_file)
seqinr::write.fasta(seq_lr_data, names(seq_lr_data), file.out = seq_lr_file)

results <- system2("java",
    args = c(
        "-jar", macse_jar,
        "-prog", "alignSequences",
        "-seq", seq_file,
        "-seq_lr", seq_lr_file,
        "-out_NT", out_nt_file,
        "-out_AA", out_aa_file
        ),
    stdout = FALSE
)

if(results) {
    quit(save = "no", status = results)
}

# MACSE uses "!" to mark deletions of one or two nucleotides that create frameshifts
# Due to downstream analyses, we will replace these with "-"

dna_aligned <- seqinr::read.fasta(out_nt_file, seqtype = "DNA", forceDNAtolower = FALSE)

dna_aligned <- dna_aligned |> map(\(x) {x[x == "!"] <- "-"; x})

seqinr::write.fasta(dna_aligned, names(dna_aligned), file.out = output_path)
