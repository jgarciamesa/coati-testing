# Get gene ID list from Ensembl

library(biomaRt)
library(dplyr, warn.conflicts = FALSE)

# As explained in bioconductor.org, if R version is 3.6 or above:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("biomaRt")
#
# For documentation of this version:
# browseVignettes("biomaRT")

get_geneId = function(species, output) {
	dir.create(file.path("raw_data",species), showWarnings = FALSE,
	           recursive = TRUE)

	#output = paste0("raw_data/",species,"_geneId.tsv")

	# Query to Ensembl and select homo sapiens gene dataset
	ensembl_h = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	# Filters to use
	filters = c("with_ccds")
	# Attributes to retrieve
	if(species == "gorilla") {
		attributes = c("ensembl_gene_id", "ggorilla_homolog_ensembl_gene", "ggorilla_homolog_orthology_type")
	} else if(species == "mouse") {
		attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type")
	} else if(species == "dmelanogaster") {
		attributes = c("ensembl_gene_id", "dmelanogaster_homolog_ensembl_gene", "dmelanogaster_homolog_orthology_type")
	} else {
		stop(paste0("species ", species," not supported."))
	}
	values = list(TRUE)#, TRUE)

	human_otherspecies = getBM(attributes, filters, values,
		mart = ensembl_h, uniqueRows = TRUE)

	# Filter one one2one (121) orthologs
	if(species == "gorilla") {
		human_otherspecies_121_orth = filter(human_otherspecies, ggorilla_homolog_orthology_type == "ortholog_one2one")
	} else if(species == "mouse") {
		human_otherspecies_121_orth = filter(human_otherspecies, mmusculus_homolog_orthology_type == "ortholog_one2one")
	} else if(species == "dmelanogaster") {
		human_otherspecies_121_orth = filter(human_otherspecies, dmelanogaster_homolog_orthology_type == "ortholog_one2one")
	} else {
	    stop(paste0("species ", species, " not supported."))
	}

	# sort
	genes = human_otherspecies_121_orth[order(human_otherspecies_121_orth$ensembl_gene_id),]
	genes = genes[, -3]

	write.table(genes, file = output, quote = FALSE, sep = "\t",
	            row.names = FALSE, col.names = FALSE)
}

if(!interactive()) {
	ARGS = commandArgs(trailingOnly = TRUE)
	get_geneId(species = ARGS[1], output = ARGS[2])
}

