options(tidyverse.quiet = TRUE,
    readr.show_col_types = FALSE,
    conflicts.policy = list(warn = FALSE ))

library(tidyverse)
library(stringr)
library(here)

# convert cigar characters that are "=" and "X" to "M" and reencode
simplify_cigar <- function(x) {
    len <- str_extract_all(x, "\\d+")
    val <- str_extract_all(x, "\\D+")

    result <- map2_chr(len, val, function(n, v) {
        if(length(n) == 1L && is.na(n)) {
            return(NA_character_)
        }
        v[v == "=" | v == "X"] <- "M"
        m <- rep.int(v, as.integer(n))

        m <- rle(m)
        str_c(m$lengths, m$values, collapse="")
    })
    result
}

# collect cigar patterns from the five primary aligners
gather_main <- function() {
    aln_data <- read_csv(here("results", "raw_fasta_aligned_stats.csv.gz"))

    aln_data <- aln_data |> 
        filter(method %in%
            c("clustalo", "coati-tri-mg", "macse", "mafft", "prank")) |>
        mutate(cigar = simplify_cigar(cigar))

    gap_data <- aln_data |>
        filter(str_detect(cigar, "[ID]")) |>
        select(gene, method, cigar)

    gapped_genes <- sort(unique(gap_data$gene))
    gapless_genes <- setdiff(sort(unique(aln_data$gene)), gapped_genes)

    list(cigar = gap_data,
        gapped = gapped_genes,
        gapless = gapless_genes)
}

if(!rlang::is_interactive()) {
    results <- gather_main()

    write_csv(results$cigar, here("benchmark_fasta", "gap_patterns.csv.gz"))

    write_lines(results$gapped, here("benchmark_fasta", "gapped_genes.txt"))
    write_lines(results$gapless, here("benchmark_fasta", "gapless_genes.txt"))
}
