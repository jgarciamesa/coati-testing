################################################################################
# Global variables                                                             #
################################################################################
RSCRIPT = Rscript --vanilla
SHELL = /bin/bash -e -o pipefail

# GENE_IDS variable
include raw_fasta/gene_id_list.mk

.DELETE_ON_ERROR:

################################################################################
# STEP 1: Download empirical CDS sequences and create unaligned fasta files    #
################################################################################

step/1_download_raw_fasta:

.PHONY: step/1_download_raw_fasta

# Download CDS sequences for every human-gorilla orthologous gene pair where
#   - Both genes are on autosomes and encode proteins
#   - The human sequence has a CCDS annotation
#   - There are one-to-one orthologs
#   - Use the canonical isoforms for both species
#   - Total nucleotide length of both isoforms is <= 6000

raw_fasta/hs-gg_gene_pairs.csv.gz: scripts/create_gene_table.R
	@echo "Downloading gene ids from ENSEMBL..."
	$(RSCRIPT) scripts/create_gene_table.R $@

raw_fasta/Homo_sapiens.GRCh38.cds.all.fa.gz:
	curl https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz > $@

raw_fasta/Gorilla_gorilla.gorGor4.cds.all.fa.gz:
	curl https://ftp.ensembl.org/pub/release-110/fasta/gorilla_gorilla/cds/Gorilla_gorilla.gorGor4.cds.all.fa.gz > $@

raw_fasta/.script_done: scripts/write_raw_fasta.R raw_fasta/hs-gg_gene_pairs.csv.gz \
	raw_fasta/Homo_sapiens.GRCh38.cds.all.fa.gz raw_fasta/Gorilla_gorilla.gorGor4.cds.all.fa.gz
	@echo "Creating unaligned FASTA files for every gene...."	
	$(RSCRIPT) scripts/write_raw_fasta.R raw_fasta/hs-gg_gene_pairs.csv.gz \
		raw_fasta/Homo_sapiens.GRCh38.cds.all.fa.gz raw_fasta/Gorilla_gorilla.gorGor4.cds.all.fa.gz \
		raw_fasta
	touch $@

raw_fasta/gene_id_list.mk: | raw_fasta/hs-gg_gene_pairs.csv.gz
	echo "GENE_IDS = \\" > $@
	zcat raw_fasta/hs-gg_gene_pairs.csv.gz | awk -F, 'NR > 1 { print "    " $$1 " \\" }' >> $@
	echo "" >> $@

raw_fasta/gene_id_list.txt: | raw_fasta/hs-gg_gene_pairs.csv.gz
	zcat raw_fasta/hs-gg_gene_pairs.csv.gz | awk -F, 'NR > 1 { print $$1 }' | sort > $@

results/raw_fasta_metrics.csv: scripts/create_raw_fasta_metrics.R raw_fasta/.script_done
	$(RSCRIPT) scripts/create_raw_fasta_metrics.R raw_fasta $@

raw_fasta: raw_fasta/.script_done raw_fasta/gene_id_list.mk results/raw_fasta_metrics.csv

.PHONY: raw_fasta

step/1_download_raw_fasta: raw_fasta

################################################################################
# STEP 2: Align empirical CDS sequences                                        #
################################################################################

METHODS = clustalo coati-dna coati-mar-ecm coati-mar-mg coati-tri-ecm \
    coati-tri-mg macse mafft prank rev-coati-tri-mg

# aligner commands are stored in a different file so they can be used
# without make having to process the full dependency list
include align.mk

step/2_empirical_alignments: $(addprefix step/2_empirical_alignments_,$(METHODS))

.PHONY: step/2_empirical_alignments

# NOTE: The empirical alignment rules are at the end of this file because they
# require secondary expansion.

################################################################################
# STEP 3: Generate alignment statistics                                        #
################################################################################

raw_fasta_aligned/%/stats.csv.gz: scripts/aln_data.R
	$(RSCRIPT) scripts/aln_data.R raw_fasta_aligned/$* $* $@

results/raw_fasta_aligned_stats.csv.gz: $(addprefix raw_fasta_aligned/,$(addsuffix /stats.csv.gz,$(METHODS)))
	$(RSCRIPT) -e "commandArgs(TRUE) |> \
	purrr::map(\(x) readr::read_csv(x, col_types = readr::cols(.default = \"c\"))) |> \
	purrr::list_rbind() |> readr::write_csv(\"$@\")" $^

step/3_generate_stats: results/raw_fasta_aligned_stats.csv.gz

.PHONY: step/3_generate_stats

################################################################################
# STEP 4: Create benchmark alignments                                          #
################################################################################

benchmark_fasta/gap_patterns.csv.gz benchmark_fasta/gapped_genes.txt \
benchmark_fasta/gapless_genes.txt: results/raw_fasta_aligned_stats.csv.gz \
scripts/gather_benchmark_data.R
	$(RSCRIPT) scripts/gather_benchmark_data.R

benchmark_fasta/.script_done: benchmark_fasta/gap_patterns.csv.gz benchmark_fasta/gapless_genes.txt scripts/simulate_benchmarks.R
	$(RSCRIPT) scripts/simulate_benchmarks.R
	touch $@

benchmark_fasta/test_id_list.mk: | benchmark_fasta/gapless_genes.txt
	echo "TEST_IDS = \\" > $@
	cat $| | awk '{ printf "    TEST%06d \\\n", NR }' >> $@
	echo "" >> $@

benchmark_fasta/test_id_list.txt: | benchmark_fasta/gapless_genes.txt
	cat $| | awk '{ printf "TEST%06d\n", NR }' > $@

step/4_simulate_benchmarks: benchmark_fasta/.script_done

.PHONY: step/4_simulate_benchmarks

################################################################################
# STEP 5: Align benchmark alignments                                           #
################################################################################

include benchmark_fasta/test_id_list.mk

step/5_benchmark_alignments: $(addprefix step/5_benchmark_alignments_,$(METHODS))

.PHONY: step/5_benchmark_alignments

################################################################################
# STEP 6: Generate benchmark statistics                                        #
################################################################################

benchmark_fasta/stats.csv.gz: scripts/aln_data.R
	$(RSCRIPT) scripts/aln_data.R benchmark_fasta benchmark $@

benchmark_fasta_aligned/%/stats.csv.gz: scripts/aln_data.R
	$(RSCRIPT) scripts/aln_data.R raw_fasta_aligned/$* $* $@

results/benchmark_fasta_aligned_stats.csv.gz: benchmark_fasta/stats.csv.gz \
$(addprefix benchmark_fasta_aligned/,$(addsuffix /stats.csv.gz,$(METHODS)))
	$(RSCRIPT) -e "commandArgs(TRUE) |> \
	purrr::map(\(x) readr::read_csv(x, col_types = readr::cols(.default = \"c\"))) |> \
	purrr::list_rbind() |> readr::write_csv(\"$@\")" $^

step/6_generate_stats: results/benchmark_fasta_aligned_stats.csv.gz

.PHONY: step/6_generate_stats


################################################################################

################################################################################
# TARGETS NEEDING SECONDEXPANSION                                              #
################################################################################

.SECONDEXPANSION:

## EMPIRICAL ALIGNMENT #########################################################

$(addprefix step/2_empirical_alignments_,$(METHODS)): step/2_empirical_alignments_%: \
	$$(addprefix raw_fasta_aligned/$$*/,$$(addsuffix .$$*.fasta,$$(GENE_IDS)))

.PHONY: $(addprefix step/2_empirical_alignments_,$(METHODS))

## BENCHMARK ALIGNMENT #########################################################

$(addprefix step/5_benchmark_alignments_,$(METHODS)): step/5_benchmark_alignments_%: \
	$$(addprefix benchmark_fasta_aligned/$$*/,$$(addsuffix .$$*.fasta,$$(TEST_IDS)))

.PHONY: $(addprefix step/5_benchmark_alignments_,$(METHODS))

# find raw_fasta_aligned -type f -name '*.fasta' -exec touch {} +
