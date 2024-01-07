library(seqinr)
library(stringr)
library(purrr)

ARGS <- commandArgs(trailingOnly = TRUE)

clustalo_bin <- ARGS[[1]]
input_path <- ARGS[[2]]
output_path <- ARGS[[3]]

universal_genetic_code <- function(remove_stops = FALSE) {
    # genetic code in TCGA order
    aa    <- "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    base1 <- "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
    base2 <- "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    base3 <- "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

    aa <- str_split_1(aa, "")
    base1 <- str_split_1(base1, "")
    base2 <- str_split_1(base2, "")
    base3 <- str_split_1(base3, "")

    names(aa) <- str_c(base1, base2, base3)
    
    if(isTRUE(remove_stops)) {
        aa <- aa[aa != "*"]
    }

    # return code in ACGT order
    aa[order(names(aa))]
}

array_to_codons <- function(x) {
    s <- seq(1, length(x), 3)
    str_c(x[s], x[s+1], x[s+2])
}

uni64 <- universal_genetic_code()

translate_seq <- function(x) {
    aa <- uni64[array_to_codons(x)]
    aa[is.na(aa)] <- "X"
    aa
}

input_data <- seqinr::read.fasta(input_path, seqtype = "DNA", forceDNAtolower = FALSE)

aa_data <- map(input_data, translate_seq)

aa_fasta_str <- aa_data |> imap(\(x, y) {
        s <- str_flatten(x)
        str_glue(">{y}\n{s}\n\n\n")
    }) |> str_flatten()


results <- system2(clustalo_bin,
    args = c( "--seqtype=Protein",
              "--infile=-",
              "--output-order=input-order"
              ),
    input = aa_fasta_str,
    stdout = TRUE )

if(is.integer(results)) {
    quit(save = "no", status = results)
} else if(is.integer(attr(results, "status"))) {
    quit(save = "no", status = attr(results, "status"))
}

aa_aligned <- seqinr::read.fasta(textConnection(results), seqtype = "AA")

untranslate_aln <- function(x, y) {
    xx <- rep(x, each = 3)
    # if y is too short, pad it with gaps
    yy <- c(y, rep("-", sum(xx != "-") - length(y)))
    w <- which(xx != "-")

    xx[ xx != "-" ] <- yy
    xx
}

dna_aligned <- map2(aa_aligned, input_data, untranslate_aln)

# remove any resiude-free columns
no_nucs <- dna_aligned |>
    map(\(x) which(x == "-")) |>
    reduce(intersect)

if(length(no_nucs)) {
    dna_aligned <- dna_aligned |> map(\(x) x[-no_nucs])
}

seqinr::write.fasta(dna_aligned, names(dna_aligned), file.out = output_path)
