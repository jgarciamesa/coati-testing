# insert gap patterns into gapless genes
options(tidyverse.quiet = TRUE,
    readr.show_col_types = FALSE,
    conflicts.policy = list(warn = FALSE))

library(seqinr)
library(tidyverse)
library(stringr)
library(here)
library(fs)
library(cli)

NUCS <- c("A" = 0.308, "C" = 0.185, "G" = 0.199, "T" = 0.308)

# possible spacings between gaps before glue starts
SPACINGS <- c(96, 48, 24, 12, 6, 3)

# construct stationary codon frequency
codons <- expand_grid(X = names(NUCS), Y = names(NUCS), Z = names(NUCS))
codons <- codons |> mutate(
    codon = str_c(X, Y, Z),
    weight = NUCS[X]*NUCS[Y]*NUCS[Z])
codons <- codons |> filter(!(codon %in% c("TAG", "TAA", "TGA")))

codons <- set_names(codons$weight, codons$codon)
codons <- codons / sum(codons)

simulate_alignment <- function(sequences, gap_patterns) {
    A <- sequences[[1]]
    B <- sequences[[2]]

    # sanity checks
    stopifnot(length(A) == length(B), length(A) %% 3 == 0)

    # randomize order of cigar strings
    pat <- sample(gap_patterns$cigar)

    for(p in pat) {
        op <- str_extract_all(p, "\\D+")[[1]]
        len <- as.integer(str_extract_all(p, "\\d+")[[1]])

        dlen <- sum(len[op == "D"])

        # identify size of static space, reducing criteria if needed
        spacing  <- detect(SPACINGS, \(sp) {
            any(op == "M" & len >= sp) &&
            sum(pmin(len[op == "M"], sp)) + dlen <= length(A)
        })
        if(is.null(spacing)) {
            # This pattern is too big, go to next pattern.
            next
        }

        is_flex <- op == "M" & len >= spacing
        if(any(is_flex)) {
            # shorten flexible matches, preserving phase
            len[is_flex] <- spacing + (len[is_flex] %% 3)

            # static space
            static_len <- sum(len[op %in% c("M", "D")])
            # flexible space
            flex_len <- length(A) - static_len
            stopifnot(flex_len %% 3 == 0)

            # divide flexible space among each flexible match
            flex_len <- flex_len %/% 3
            w <- sample(0:flex_len, size = sum(is_flex)-1, replace = TRUE)
            w <- sort(w)
            w <- diff(c(0,w,flex_len))
            len[is_flex] <- len[is_flex] + 3L*w
        }

        # sanity checks
        stopifnot(sum(len[op %in% c("M", "D")]) == length(A))
        stopifnot(all(op %in% c("M", "D", "I")))

        # write new alignments
        aln <- rep.int(op, len)
        # copy the ancestor's nucleotides for every position
        # that is a match or deletion
        AA <- rep("-", length(aln))
        AA[aln %in% c("M", "D")] <- A

        # initialize descent sequence with random codons
        # this ensures that insertions obey phase
        B0 <- sample(names(codons), size = ceiling(sum(aln %in% c("M", "I"))/3),
            replace = TRUE, prob = codons)
        B0 <- str_split(B0, "") |> list_c()
        B0 <- B0[seq_len(sum(aln %in% c("M", "I")))]
        BB <- rep("N", length(aln))
        BB[aln %in% c("M", "I")] <- B0

        # Copy over descendant nucleotides
        BB[aln %in% c("M", "D")] <- B
        # mark deletions
        BB[aln == "D"] <- "-"

        new_cigar <- str_c(len, op, collapse="")
        ret <- list(
            cigar = p,
            anc = AA,
            dec = BB,
            aln = new_cigar
        )
        return(ret)
    }
    stop("No compatible gap pattern found.")
}

anc_length <- function(x) {
    len <- str_extract_all(x, "\\d+")
    val <- str_extract_all(x, "\\D+")

    result <- map2_int(len, val, function(n, op) {
        sum(as.integer(n)[op %in% c("M", "D")])
    })
    result
}

sim_main <- function() {
    # Make results reproducible by setting seed
    seed <- digest::digest2int("simulate_benchmarks.R", seed = 1L)
    withr::local_seed(seed, .rng_kind = "default")

    # Load simulation "parameters"
    gapless_genes <- read_lines(here("benchmark_fasta",
        "gapless_genes.txt"))
    gap_data <- read_csv(here("benchmark_fasta",
        "gap_patterns.csv.gz"))
    gap_data <- gap_data |> mutate(anc_len = anc_length(cigar))

    # split cigars into groups by method
    gap_pat <- split(select(gap_data, cigar, anc_len), gap_data$method)

    # assign each gene to use one of the groups
    # we will use each group evenly and then randomize the association
    g <- sample(rep(seq_along(gap_pat), length.out = length(gapless_genes)))

    # clear destination folders (removes any leftovers in case the number of
    # output files changes)
    f <- dir_ls(here("benchmark_fasta"), regexp="\\bTEST\\d+[.]fasta$")
    file_delete(f)
    f <- dir_ls(here("benchmark_fasta_nogaps"), regexp="\\bTEST\\d+[.]fasta$")
    file_delete(f)

    # process each gene
    tab <- NULL
    cli_progress_bar("Simulating Benchmark Alignments", total = length(g))
    for(i in seq_along(g)) {
        gene <- gapless_genes[[i]]
        model <- g[[i]]
        fasta <- read.fasta(here("raw_fasta", str_glue("{gene}.fasta")),
            forceDNAtolower = FALSE, set.attributes = FALSE)
        result <- simulate_alignment(fasta, gap_pat[[model]])

        stem <- sprintf("TEST%06d", i)

        # Aligned Sequences
        seqs <- list(result$anc, result$dec)
        write.fasta(seqs, str_glue("{stem}_{n}", n = c("A", "D")),
            file.out = here("benchmark_fasta", str_glue("{stem}.fasta")))
        
        # Unaligned Sequences
        seqs <- map(seqs, \(x) x[x != "-"])
        write.fasta(seqs, str_glue("{stem}_{n}", n = c("A", "D")),
            file.out = here("benchmark_fasta_nogaps",
                str_glue("{stem}.nogaps.fasta")))

        r <- tibble_row(gene = gene,
            model = names(gap_pat)[[model]],
            pattern = result$cigar,
            aln = result$aln)

        tab <- bind_rows(tab, r)
        cli_progress_update()
    }
    cli_progress_done()

    # save meta data
    write_csv(tab, here("benchmark_fasta", "metadata.csv.gz"))
}

if(!rlang::is_interactive()) {
    sim_main()
}
