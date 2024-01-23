options(tidyverse.quiet = TRUE,
    readr.show_col_types = FALSE,
    conflicts.policy = list(warn = FALSE ))

library(tidyverse)
library(stringr)
library(seqinr)
library(ape)
library(fs)
library(here)

COATI_BIN <- here("bin", "coati-alignpair")

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

calculate_score <- function(file) {
    results <- system2(COATI_BIN,
        args = c("-s", "-m", "mar-mg",
            file), stdout = TRUE)
    if(is.integer(results)) {
        return(NA_real_)
    } else if(is.integer(attr(results, "status"))) {
        return(NA_real_)
    } else if(length(results) == 0) {
        return(NA_real_)
    }
    as.double(results)
}

read_fasta <- purrr::possibly(seqinr::read.fasta, quiet = TRUE)

process_file <- function(path, method) {
    file <- path_file(path) |> path_ext_remove()
    gene <- path_ext_remove(file)
    method <- method %||% path_ext(file)

    x <- read_fasta(path, forceDNAtolower = FALSE)

    if(is.null(x)) {
        return(list(
            gene = gene,
            method = method,
            cigar = NA_character_,
            score = NA_real_,
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
        score <- calculate_score(path)
    } else {
        cigar <- NA_character_
        omega <- NA_real_
        k2p <- NA_real_
        score <- NA_real_
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
        score = score,
        omega = omega,
        k2p = c(k2p),
        checksum = checksum,
        error = error
    )

}

aln_data_main <- function(input_dir, method) {
    files <- dir_ls(input_dir, type = "file", regexp = "\\b(ENS|TEST)*.fasta$")

    dat <- files |> map(\(x) process_file(x, method), .progress = TRUE)
    dat <- dat |> map(as_tibble_row) |>
        list_rbind()

    dat
}

if(!rlang::is_interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    dat <- aln_data_main(args[[1]], args[[2]])
    write_csv(dat, args[[3]])
}
