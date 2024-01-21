options(tidyverse.quiet = TRUE,
    readr.show_col_types = FALSE,
    conflicts.policy = list(warn = FALSE ))

library(tidyverse)
library(stringr)
library(seqinr)
library(ape)
library(fs)

calculate_omega <- function(aln) { #dN/dS
    k <- kaks(aln)
    if(length(k) == 1) {
        return(NA)
    }
    ka <- max(k$ka[1], 0)
    ks <- max(k$ks[1], 0)

    if(ka == 0 && ks == 0) {
        return(0)
    } else if(ks == 0) {
        return(10)
    } else {
        return(ka/ks)
    }
}

create_cigar <- function(anc, dec, use_match = FALSE) {
    x <- case_when(
        anc == "-" ~ "I",
        dec == "-" ~ "D",
        use_match ~ "M",
        anc == dec ~ "=",
        anc == "N" ~ "=",
        dec == "N" ~ "=",
        TRUE ~ "X"
    )
    x <- rle(x)
    str_c(x$lengths, x$values, collapse="")
}

# path <- "raw_fasta_aligned/coati-tri-mg/ENSG00000000419.coati-tri-mg.fasta"

process_file <- function(path, method = NULL) {
    file <- path_file(path) |> path_ext_remove()
    gene <- path_ext_remove(file)
    method <- method %||% path_ext(file)

    x <- try(seqinr::read.fasta(path, forceDNAtolower = FALSE), silent = TRUE)

    if(inherits(x, "try-error")) {
        return(list(
            gene = gene,
            method = method,
            cigar = NA_character_,
            omega = NA_real_,
            k2p = NA_real_,
            checksum = NA_character_,
            error = "blank_file"
        ))
    }

    anc <- x[[1]]
    dec <- x[[2]]

    # Create an alignment using only columns with ancestor
    a <- map(x, \(y) y[anc != "-"]) |> ape::as.alignment()

    allowed_chars <- c("A", "C", "G", "T", "N", "-")

    # Sanity checks
    error <- NA_character_
    if(length(anc) != length(dec)) {
        error <- "different_lengths"
    } else if(!all(anc %in% allowed_chars & dec %in% allowed_chars)) {
        error <- "unallowed_chars"
    } else if(any(anc == "-" & dec == "-")) {
        error <- "empty_columns"
    }
    if(is.na(error)) {
        cigar <- create_cigar(anc, dec)
        omega <- calculate_omega(a)
        k2p <- dist.dna(as.DNAbin(a), model = "K80")
    } else {
        cigar <- NA_character_
        omega <- NA_real_
        k2p <- NA_real_
    }

    # Create a checksum of the unaligned sequences
    anc_s <- str_flatten(anc[anc != "-"])
    dec_s <- str_flatten(dec[dec != "-"])
    s <- str_glue(">{getName(anc)}\n{anc_s}\n>{getName(dec)}\n{dec_s}\n")
    checksum <- digest::digest(s, "xxhash32", serialize = FALSE)

    list(
        gene = gene,
        method = method,
        cigar = cigar,
        omega = omega,
        k2p = c(k2p),
        checksum = checksum,
        error = error
    )

}

aln_data_main <- function(input_dir, method = NULL) {
    files <- dir_ls(input_dir, type = "file", regexp = "\\b(ENS|TEST)*.fasta$")

    dat <- files |> map(\(x) process_file(x, method))
    dat <- dat |> map(as_tibble_row) |>
        list_rbind()

    dat
}

if(!rlang::is_interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    if(length(args) >= 2) {
        method <- args[[2]]
    } else {
        method <- NULL
    }
    dat <- aln_data_main(args[[1]], method)
    cat(format_csv(dat))
}
