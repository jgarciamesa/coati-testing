library(fs)
library(seqinr)
library(readr)
library(purrr)
library(stringr)

table_path <- "raw_fasta/hs-gg_gene_pairs.csv.gz"
hsapiens_path <- "raw_fasta/Homo_sapiens.GRCh38.cds.all.fa.gz"
ggorilla_path <- "raw_fasta/Gorilla_gorilla.gorGor4.cds.all.fa.gz"
output_path <- "raw_fasta"

write_raw_fasta_main <- function(table_path, hsapiens_path, ggorilla_path, output_path) {

    dir_create(output_path)

    ortho_data <- read_csv(table_path, show_col_types = FALSE)
    hsapiens_data <- read.fasta(hsapiens_path)
    ggorilla_data <- read.fasta(ggorilla_path)

    ortho_data <- transpose(ortho_data)

    write_file <- function(x) {
        fasta <- path(output_path, x$hsapiens_ensembl_gene_id, ext = "fasta")

        h <- hsapiens_data[x$hsapiens_ensembl_transcript_id_version]
        g <- ggorilla_data[x$ggorilla_ensembl_transcript_id_version]

        hg <- map(c(h, g), \(x) str_to_upper(str_flatten(x)))

        write.fasta(hg, names(hg), file.out = fasta, as.string = TRUE)
    }

    walk(ortho_data, write_file, .progress = TRUE)
}


if(!interactive()) {
    ARGS <- commandArgs(trailingOnly = TRUE)
    write_raw_fasta_main(ARGS[1], ARGS[2], ARGS[3], ARGS[4])
}