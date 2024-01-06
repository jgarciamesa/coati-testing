options(tidyverse.quiet = TRUE, readr.show_col_types = FALSE)
library(seqinr)
library(tidyverse)
library(fs)

input_dir <- "raw_fasta"

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

create_metrics <- function(input_dir) {
    files <- dir_ls(input_dir, type = "file", glob = "*\\bENS*.fasta")

    uni64 <- universal_genetic_code()

    process_seq <- function(x) {
        aa <- uni64[array_to_codons(x)]
        good_length <- (length(x) %% 3) == 0
        num_early_stops <- sum(head(aa, -1) == "*", na.rm = TRUE)
        num_amb_chars <- sum(!(x %in% c("A", "C", "G", "T")))

        list(good_length = good_length,
            num_early_stops = num_early_stops,
            num_amb_chars = num_amb_chars,
            length = length(x))
    }

    process_file <- function(x) {
        y <- read.fasta(x, forceDNAtolower = FALSE)

        # human data
        a <- list(
            file = path_file(x) |> path_ext_remove(),
            species = "human",
            transcript = names(y)[[1]]
        )
        a <- c(a, process_seq(y[[1]]))

        #gorilla data
        b <- list(
            file = path_file(x) |> path_ext_remove(),
            species = "gorilla",
            transcript = names(y)[[2]]
        )
        b <- c(b, process_seq(y[[2]]))

        a <- tibble_row(!!!a)
        b <- tibble_row(!!!b)

        bind_rows(a, b)
    }

    results <- map(files, process_file, .progress = TRUE)
    results <- list_rbind(results)

    results
}

if(!rlang::is_interactive()) {
    ARGS <- commandArgs(trailingOnly = TRUE)
    dat <- create_metrics(ARGS[1])
    readr::write_csv(dat, ARGS[2])
}
