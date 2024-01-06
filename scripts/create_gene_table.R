# Get gene ID list from Ensembl

library(biomaRt)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)

TOTAL_LEN_LIMIT <- 6000

# As explained in bioconductor.org, if R version is 3.6 or above:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("biomaRt")
#
# For documentation of this version:
# browseVignettes("biomaRT")

get_gene_ids <- function() {
    # Query to Ensembl and select homo sapiens gene dataset
    ensembl_h <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    # Filters to use
    filters <- list(
        with_ccds = TRUE,
        biotype = "protein_coding",
        transcript_biotype = "protein_coding",
        chromosome_name = 1:22
    )
    # Attributes to retrieve
    attributes <- c("ensembl_gene_id", "ensembl_gene_id_version",
        "ggorilla_homolog_ensembl_gene", "ggorilla_homolog_orthology_type")

    results <- getBM(attributes, filters,
        mart = ensembl_h, uniqueRows = TRUE)

    # Filter one one2one (121) orthologs
    # Sort by ensembl_gene_id
    results <- results |> 
        filter(ggorilla_homolog_orthology_type == "ortholog_one2one") |>
        arrange(ensembl_gene_id) |>
        select(-ggorilla_homolog_orthology_type)
    
    as_tibble(results)
}

get_human_canonical <- function() {
    ensembl_h <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    # Filters to use
    filters <- list(
        with_ccds = TRUE,
        biotype = "protein_coding",
        transcript_biotype = "protein_coding",
        transcript_is_canonical = TRUE,
        chromosome_name = 1:22
    )
    # Attributes to retrieve
    attributes <- c("ensembl_gene_id", "ensembl_gene_id_version",
        "ensembl_transcript_id", "ensembl_transcript_id_version",
        "cds_length")

    results <- getBM(attributes, filters,
        mart = ensembl_h, uniqueRows = TRUE)

    results <- results |>
        arrange(ensembl_gene_id)

    as_tibble(results)
}

get_gorilla_canonical <- function() {
    ensembl_h <- useMart("ensembl", dataset = "ggorilla_gene_ensembl")
    # Filters to use
    filters <- list(
        biotype = "protein_coding",
        transcript_biotype = "protein_coding",
        transcript_is_canonical = TRUE,
        chromosome_name = c(1, "2A", "2B", 3:22)
    )
    # Attributes to retrieve
    attributes <- c("ensembl_gene_id", "ensembl_gene_id_version",
        "ensembl_transcript_id", "ensembl_transcript_id_version",
        "cds_length")

    results <- getBM(attributes, filters,
        mart = ensembl_h, uniqueRows = TRUE)

    results <- results |>
        arrange(ensembl_gene_id)

    as_tibble(results)
}

get_gene_table_main <- function() {
    dat_gene_ids <- get_gene_ids()
    dat_hum_trans <- get_human_canonical()
    dat_gor_trans <- get_gorilla_canonical()

    results <- left_join(dat_gene_ids, dat_hum_trans) |>
        rename(hsapiens_ensembl_gene_id = ensembl_gene_id,
            hsapiens_ensembl_gene_id_version = ensembl_gene_id_version,
            hsapiens_ensembl_transcript_id = ensembl_transcript_id,
            hsapiens_ensembl_transcript_id_version = ensembl_transcript_id_version,
            hsapiens_cds_length = cds_length,
            ggorilla_ensembl_gene_id = ggorilla_homolog_ensembl_gene)

    dat_gor_trans <- dat_gor_trans |> rename(
            ggorilla_ensembl_gene_id = ensembl_gene_id,
            ggorilla_ensembl_gene_id_version = ensembl_gene_id_version,
            ggorilla_ensembl_transcript_id = ensembl_transcript_id,
            ggorilla_ensembl_transcript_id_version = ensembl_transcript_id_version,
            ggorilla_cds_length = cds_length
        )

    results <- left_join(results, dat_gor_trans)

    results <- results |> relocate(ggorilla_ensembl_gene_id,
        .after = hsapiens_cds_length)

    drop_na(results) |>
        filter(hsapiens_cds_length + ggorilla_cds_length <= TOTAL_LEN_LIMIT)
}

if(!interactive()) {
    ARGS <- commandArgs(trailingOnly = TRUE)
    dat <- get_gene_table_main()
    readr::write_csv(dat, ARGS[1])
}
