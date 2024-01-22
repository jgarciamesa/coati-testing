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

results/raw_fasta_metrics.csv: scripts/create_raw_fasta_metrics.R raw_fasta/.script_done
	$(RSCRIPT) scripts/create_raw_fasta_metrics.R raw_fasta $@

raw_fasta: raw_fasta/.script_done raw_fasta/gene_id_list.mk results/raw_fasta_metrics.csv

.PHONY: raw_fasta

step/1_download_raw_fasta: raw_fasta

################################################################################
# STEP 2: Align empirical CDS sequences                                        #
################################################################################

# aligner commands are stored in a different file so they can be used
# without make having to process the full dependency list
include align.mk

step/2_empirical_alignments:

.PHONY: step/2_empirical_alignments

## COATI TRI-MG ################################################################

RAWALN_FILES_COATI_TRIMG=$(addprefix raw_fasta_aligned/coati-tri-mg/,\
    $(addsuffix .coati-tri-mg.fasta, $(GENE_IDS)))

raw_fasta_aligned/coati-tri-mg: $(RAWALN_FILES_COATI_TRIMG)

.PHONY: raw_fasta_aligned/coati-tri-mg

step/2_empirical_alignments: raw_fasta_aligned/coati-tri-mg

raw_fasta_aligned/coati-tri-mg/coati-tri-mg.archive.tar.gz: raw_fasta/hs-gg_gene_pairs.csv.gz $(RAWALN_FILES_COATI_TRIMG)
	zcat $< | awk -F, 'NR > 1 { print "raw_fasta_aligned/coati-tri-mg/" $$1 ".coati-tri-mg.fasta" }' | tar cvzf $@ -T -

step/2_empirical_alignments: raw_fasta_aligned/coati-tri-mg/coati-tri-mg.archive.tar.gz

## COATI TRI-ECM ###############################################################

RAWALN_FILES_COATI_TRIECM=$(addprefix raw_fasta_aligned/coati-tri-ecm/,\
    $(addsuffix .coati-tri-ecm.fasta, $(GENE_IDS)))

raw_fasta_aligned/coati-tri-ecm: $(RAWALN_FILES_COATI_TRIECM)

.PHONY: raw_fasta_aligned/coati-tri-ecm

step/2_empirical_alignments: raw_fasta_aligned/coati-tri-ecm

raw_fasta_aligned/coati-tri-ecm/coati-tri-ecm.archive.tar.gz: raw_fasta/hs-gg_gene_pairs.csv.gz $(RAWALN_FILES_COATI_TRIECM)
	zcat $< | awk -F, 'NR > 1 { print "raw_fasta_aligned/coati-tri-ecm/" $$1 ".coati-tri-ecm.fasta" }' | tar cvzf $@ -T -

step/2_empirical_alignments: raw_fasta_aligned/coati-tri-ecm/coati-tri-ecm.archive.tar.gz

## COATI MAR-MG ################################################################

RAWALN_FILES_COATI_MARMG=$(addprefix raw_fasta_aligned/coati-mar-mg/,\
    $(addsuffix .coati-mar-mg.fasta, $(GENE_IDS)))

raw_fasta_aligned/coati-mar-mg: $(RAWALN_FILES_COATI_MARMG)

.PHONY: raw_fasta_aligned/coati-mar-mg

step/2_empirical_alignments: raw_fasta_aligned/coati-mar-mg

raw_fasta_aligned/coati-mar-mg/coati-mar-mg.archive.tar.gz: raw_fasta/hs-gg_gene_pairs.csv.gz $(RAWALN_FILES_COATI_MARMG)
	zcat $< | awk -F, 'NR > 1 { print "raw_fasta_aligned/coati-mar-mg/" $$1 ".coati-mar-mg.fasta" }' | tar cvzf $@ -T -

step/2_empirical_alignments: raw_fasta_aligned/coati-mar-mg/coati-mar-mg.archive.tar.gz

## COATI MAR-ECM ###############################################################

RAWALN_FILES_COATI_MARECM=$(addprefix raw_fasta_aligned/coati-mar-ecm/,\
    $(addsuffix .coati-mar-ecm.fasta, $(GENE_IDS)))

raw_fasta_aligned/coati-mar-ecm: $(RAWALN_FILES_COATI_MARECM)

.PHONY: raw_fasta_aligned/coati-mar-ecm

step/2_empirical_alignments: raw_fasta_aligned/coati-mar-ecm

raw_fasta_aligned/coati-mar-ecm/coati-mar-ecm.archive.tar.gz: raw_fasta/hs-gg_gene_pairs.csv.gz $(RAWALN_FILES_COATI_MARECM)
	zcat $< | awk -F, 'NR > 1 { print "raw_fasta_aligned/coati-mar-ecm/" $$1 ".coati-mar-ecm.fasta" }' | tar cvzf $@ -T -

step/2_empirical_alignments: raw_fasta_aligned/coati-mar-ecm/coati-mar-ecm.archive.tar.gz

## COATI DNA ###################################################################

RAWALN_FILES_COATI_DNA=$(addprefix raw_fasta_aligned/coati-dna/,\
    $(addsuffix .coati-dna.fasta, $(GENE_IDS)))

raw_fasta_aligned/coati-dna: $(RAWALN_FILES_COATI_DNA)

.PHONY: raw_fasta_aligned/coati-dna

step/2_empirical_alignments: raw_fasta_aligned/coati-dna

raw_fasta_aligned/coati-dna/coati-dna.archive.tar.gz: raw_fasta/hs-gg_gene_pairs.csv.gz $(RAWALN_FILES_COATI_DNA)
	zcat $< | awk -F, 'NR > 1 { print "raw_fasta_aligned/coati-dna/" $$1 ".coati-dna.fasta" }' | tar cvzf $@ -T -

step/2_empirical_alignments: raw_fasta_aligned/coati-dna/coati-dna.archive.tar.gz

## REV COATI TRI-MG ############################################################

RAWALN_FILES_REV_COATI=$(addprefix raw_fasta_aligned/rev-coati-tri-mg/,\
    $(addsuffix .rev-coati-tri-mg.fasta, $(GENE_IDS)))

raw_fasta_aligned/rev-coati-tri-mg: $(RAWALN_FILES_REV_COATI)

.PHONY: raw_fasta_aligned/rev-coati-tri-mg

step/2_empirical_alignments: raw_fasta_aligned/rev-coati-tri-mg

raw_fasta_aligned/rev-coati-tri-mg/rev-coati-tri-mg.archive.tar.gz: raw_fasta/hs-gg_gene_pairs.csv.gz $(RAWALN_FILES_REV_COATI)
	zcat $< | awk -F, 'NR > 1 { print "raw_fasta_aligned/rev-coati-tri-mg/" $$1 ".rev-coati-tri-mg.fasta" }' | tar cvzf $@ -T -

step/2_empirical_alignments: raw_fasta_aligned/rev-coati-tri-mg/rev-coati-tri-mg.archive.tar.gz

## PRANK #######################################################################

RAWALN_FILES_PRANK=$(addprefix raw_fasta_aligned/prank/,\
    $(addsuffix .prank.fasta, $(GENE_IDS)))

raw_fasta_aligned/prank: $(RAWALN_FILES_PRANK)

.PHONY: raw_fasta_aligned/prank

step/2_empirical_alignments: raw_fasta_aligned/prank

raw_fasta_aligned/prank/prank.archive.tar.gz: raw_fasta/hs-gg_gene_pairs.csv.gz $(RAWALN_FILES_PRANK)
	zcat $< | awk -F, 'NR > 1 { print "raw_fasta_aligned/prank/" $$1 ".prank.fasta" }' | tar cvzf $@ -T -

step/2_empirical_alignments: raw_fasta_aligned/prank/prank.archive.tar.gz

## MAFFT #######################################################################

RAWALN_FILES_MAFFT=$(addprefix raw_fasta_aligned/mafft/,\
    $(addsuffix .mafft.fasta, $(GENE_IDS)))

raw_fasta_aligned/mafft: $(RAWALN_FILES_MAFFT)

.PHONY: raw_fasta_aligned/mafft

step/2_empirical_alignments: raw_fasta_aligned/mafft

raw_fasta_aligned/mafft/mafft.archive.tar.gz: raw_fasta/hs-gg_gene_pairs.csv.gz $(RAWALN_FILES_MAFFT)
	zcat $< | awk -F, 'NR > 1 { print "raw_fasta_aligned/mafft/" $$1 ".mafft.fasta" }' | tar cvzf $@ -T -

step/2_empirical_alignments: raw_fasta_aligned/mafft/mafft.archive.tar.gz

## CLUSTAL OMEGA ################################################################

RAWALN_FILES_CLUSTALO=$(addprefix raw_fasta_aligned/clustalo/,\
    $(addsuffix .clustalo.fasta, $(GENE_IDS)))

raw_fasta_aligned/clustalo: $(RAWALN_FILES_CLUSTALO)

.PHONY: raw_fasta_aligned/clustalo

step/2_empirical_alignments: raw_fasta_aligned/clustalo

raw_fasta_aligned/clustalo/clustalo.archive.tar.gz: raw_fasta/hs-gg_gene_pairs.csv.gz $(RAWALN_FILES_CLUSTALO)
	zcat $< | awk -F, 'NR > 1 { print "raw_fasta_aligned/clustalo/" $$1 ".clustalo.fasta" }' | tar cvzf $@ -T -

step/2_empirical_alignments: raw_fasta_aligned/clustalo/clustalo.archive.tar.gz

## MACSE #######################################################################

RAWALN_FILES_MACSE=$(addprefix raw_fasta_aligned/macse/,\
    $(addsuffix .macse.fasta, $(GENE_IDS)))

raw_fasta_aligned/macse: $(RAWALN_FILES_MACSE)

.PHONY: raw_fasta_aligned/macse

step/2_empirical_alignments: raw_fasta_aligned/macse

raw_fasta_aligned/macse/macse.archive.tar.gz: raw_fasta/hs-gg_gene_pairs.csv.gz $(RAWALN_FILES_MACSE)
	zcat $< | awk -F, 'NR > 1 { print "raw_fasta_aligned/macse/" $$1 ".macse.fasta" }' | tar cvzf $@ -T -

step/2_empirical_alignments: raw_fasta_aligned/macse/macse.archive.tar.gz

################################################################################
# STEP 2 (Alt): Extract empirical alignments from archives.                    #
################################################################################

step/2a_extract_empirical_alignments:
	tar -xvmf raw_fasta_aligned/clustalo/clustalo.archive.tar.gz
	tar -xvmf raw_fasta_aligned/coati-dna/coati-dna.archive.tar.gz
	tar -xvmf raw_fasta_aligned/coati-mar-ecm/coati-mar-ecm.archive.tar.gz
	tar -xvmf raw_fasta_aligned/coati-mar-mg/coati-mar-mg.archive.tar.gz
	tar -xvmf raw_fasta_aligned/coati-tri-ecm/coati-tri-ecm.archive.tar.gz
	tar -xvmf raw_fasta_aligned/coati-tri-mg/coati-tri-mg.archive.tar.gz
	tar -xvmf raw_fasta_aligned/macse/macse.archive.tar.gz
	tar -xvmf raw_fasta_aligned/mafft/mafft.archive.tar.gz
	tar -xvmf raw_fasta_aligned/prank/prank.archive.tar.gz
	tar -xvmf raw_fasta_aligned/rev-coati-tri-mg/rev-coati-tri-mg.archive.tar.gz

.PHONY: step/2a_extract_empirical_alignments

################################################################################
# STEP 3: Generate alignment statistics                                        #
################################################################################

METHODS = clustalo coati-dna coati-mar-ecm coati-mar-mg coati-tri-ecm \
    coati-tri-mg macse mafft prank rev-coati-tri-mg

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

## Benchmark Data ##############################################################

benchmark_fasta/gap_patterns.csv.gz benchmark_fasta/gapped_genes.txt benchmark_fasta/gapless_genes.txt: results/raw_fasta_aligned_stats.csv.gz scripts/gather_benchmark_data.R
	$(RSCRIPT) scripts/gather_benchmark_data.R

benchmark_fasta/.script_done: benchmark_fasta/gap_patterns.csv.gz benchmark_fasta/gapless_genes.txt scripts/simulate_benchmarks.R
	$(RSCRIPT) scripts/simulate_benchmarks.R
	touch $@

benchmark_fasta/stats.csv: scripts/aln_data.R
	$(RSCRIPT) scripts/aln_data.R benchmark_fasta benchmark > $@

benchmark_fasta/test_id_list.mk: | benchmark_fasta/gapless_genes.txt
	echo "TEST_IDS = \\" > $@
	cat $| | awk '{ printf "    TEST%06d \\\n", NR }' >> $@
	echo "" >> $@

benchmark_fasta/test_id_list.txt: | benchmark_fasta/gapless_genes.txt
	cat $| | awk '{ printf "TEST%06d\n", NR }' > $@

step/4_simulate_benchmarks: benchmark_fasta/.script_done

.PHONY: step/4_simulate_benchmarks

BENCHMARK_FILES=$(addprefix benchmark_fasta/,$(addsuffix .fasta, $(TEST_IDS)))

benchmark_fasta/benchmark.archive.tar.gz: benchmark_fasta/gapless_genes.txt $(BENCHMARK_FILES)
	cat $< | awk '{ printf "benchmark_fasta/TEST%06d.fasta\n", NR }' | tar cvzf $@ -T -

BENCHMARK_NOGAPS_FILES=$(addprefix benchmark_fasta_nogaps/,$(addsuffix .nogaps.fasta, $(TEST_IDS)))

benchmark_fasta_nogaps/benchmark_nogaps.archive.tar.gz: benchmark_fasta/gapless_genes.txt $(BENCHMARK_NOGAPS_FILES)
	cat $< | awk '{ printf "benchmark_fasta_nogaps/TEST%06d.nogaps.fasta\n", NR }' | tar cvzf $@ -T -

################################################################################
# STEP 4 (Alt): Extract benchmarks from archives.                              #
################################################################################

step/4a_extract_benchmarks:
	tar -xvmf benchmark_fasta/benchmark.archive.tar.gz
	tar -xvmf benchmark_fasta_nogaps/benchmark_nogaps.archive.tar.gz

.PHONY: step/4a_extract_benchmarks

################################################################################
# STEP 5: Align benchmark alignments                                           #
################################################################################

include benchmark_fasta/test_id_list.mk

step/5_benchmark_alignments:

.PHONY: step/5_benchmark_alignments

## COATI TRI-MG ################################################################

TSTALN_FILES_COATI_TRIMG=$(addprefix benchmark_fasta_aligned/coati-tri-mg/, \
    $(addsuffix .coati-tri-mg.fasta, $(TEST_IDS)))

benchmark_fasta_aligned/coati-tri-mg: $(TSTALN_FILES_COATI_TRIMG)

.PHONY: benchmark_fasta_aligned/coati-tri-mg

step/5_benchmark_alignments: benchmark_fasta_aligned/coati-tri-mg

benchmark_fasta_aligned/coati-tri-mg/coati-tri-mg.archive.tar.gz: benchmark_fasta/test_id_list.txt $(TSTALN_FILES_COATI_TRIMG)
	cat $< | awk '{ print "benchmark_fasta_aligned/coati-tri-mg/" $$1 ".coati-tri-mg.fasta" }' | tar cvzf $@ -T -

step/5_benchmark_alignments: benchmark_fasta_aligned/coati-tri-mg/coati-tri-mg.archive.tar.gz

## COATI TRI-ECM ###############################################################

TSTALN_FILES_COATI_TRIECM=$(addprefix benchmark_fasta_aligned/coati-tri-ecm/, \
    $(addsuffix .coati-tri-ecm.fasta, $(TEST_IDS)))

benchmark_fasta_aligned/coati-tri-ecm: $(TSTALN_FILES_COATI_TRIECM)

.PHONY: benchmark_fasta_aligned/coati-tri-ecm

step/5_benchmark_alignments: benchmark_fasta_aligned/coati-tri-ecm

benchmark_fasta_aligned/coati-tri-ecm/coati-tri-ecm.archive.tar.gz: benchmark_fasta/test_id_list.txt $(TSTALN_FILES_COATI_TRIECM)
	cat $< | awk '{ print "benchmark_fasta_aligned/coati-tri-ecm/" $$1 ".coati-tri-ecm.fasta" }' | tar cvzf $@ -T -

step/5_benchmark_alignments: benchmark_fasta_aligned/coati-tri-ecm/coati-tri-ecm.archive.tar.gz

## COATI MAR-MG ################################################################

TSTALN_FILES_COATI_MARMG=$(addprefix benchmark_fasta_aligned/coati-mar-mg/, \
    $(addsuffix .coati-mar-mg.fasta, $(TEST_IDS)))

benchmark_fasta_aligned/coati-mar-mg: $(TSTALN_FILES_COATI_MARMG)

.PHONY: benchmark_fasta_aligned/coati-mar-mg

step/5_benchmark_alignments: benchmark_fasta_aligned/coati-mar-mg

benchmark_fasta_aligned/coati-mar-mg/coati-mar-mg.archive.tar.gz: benchmark_fasta/test_id_list.txt $(TSTALN_FILES_COATI_MARMG)
	cat $< | awk '{ print "benchmark_fasta_aligned/coati-mar-mg/" $$1 ".coati-mar-mg.fasta" }' | tar cvzf $@ -T -

step/5_benchmark_alignments: benchmark_fasta_aligned/coati-mar-mg/coati-mar-mg.archive.tar.gz

## COATI MAR-ECM ###############################################################

TSTALN_FILES_COATI_MARECM=$(addprefix benchmark_fasta_aligned/coati-mar-ecm/, \
    $(addsuffix .coati-mar-ecm.fasta, $(TEST_IDS)))

benchmark_fasta_aligned/coati-mar-ecm: $(TSTALN_FILES_COATI_MARECM)

.PHONY: benchmark_fasta_aligned/coati-mar-ecm

step/5_benchmark_alignments: benchmark_fasta_aligned/coati-mar-ecm

benchmark_fasta_aligned/coati-mar-ecm/coati-mar-ecm.archive.tar.gz: benchmark_fasta/test_id_list.txt $(TSTALN_FILES_COATI_MARECM)
	cat $< | awk '{ print "benchmark_fasta_aligned/coati-mar-ecm/" $$1 ".coati-mar-ecm.fasta" }' | tar cvzf $@ -T -

step/5_benchmark_alignments: benchmark_fasta_aligned/coati-mar-ecm/coati-mar-ecm.archive.tar.gz

## COATI DNA ###################################################################

TSTALN_FILES_COATI_DNA=$(addprefix benchmark_fasta_aligned/coati-dna/, \
    $(addsuffix .coati-dna.fasta, $(TEST_IDS)))

benchmark_fasta_aligned/coati-dna: $(TSTALN_FILES_COATI_DNA)

.PHONY: benchmark_fasta_aligned/coati-dna

step/5_benchmark_alignments: benchmark_fasta_aligned/coati-dna

benchmark_fasta_aligned/coati-dna/coati-dna.archive.tar.gz: benchmark_fasta/test_id_list.txt $(TSTALN_FILES_COATI_DNA)
	cat $< | awk '{ print "benchmark_fasta_aligned/coati-dna/" $$1 ".coati-dna.fasta" }' | tar cvzf $@ -T -

step/5_benchmark_alignments: benchmark_fasta_aligned/coati-dna/coati-dna.archive.tar.gz

## REV COATI TRI-MG ############################################################

TSTALN_FILES_REV_COATI=$(addprefix benchmark_fasta_aligned/rev-coati-tri-mg/, \
    $(addsuffix .rev-coati-tri-mg.fasta, $(TEST_IDS)))

benchmark_fasta_aligned/rev-coati-tri-mg: $(TSTALN_FILES_REV_COATI)

.PHONY: benchmark_fasta_aligned/rev-coati-tri-mg

step/5_benchmark_alignments: benchmark_fasta_aligned/rev-coati-tri-mg

benchmark_fasta_aligned/rev-coati-tri-mg/rev-coati-tri-mg.archive.tar.gz: benchmark_fasta/test_id_list.txt $(TSTALN_FILES_REV_COATI)
	cat $< | awk '{ print "benchmark_fasta_aligned/rev-coati-tri-mg/" $$1 ".rev-coati-tri-mg.fasta" }' | tar cvzf $@ -T -

step/5_benchmark_alignments: benchmark_fasta_aligned/rev-coati-tri-mg/rev-coati-tri-mg.archive.tar.gz

## PRANK #######################################################################

TSTALN_FILES_PRANK=$(addprefix benchmark_fasta_aligned/prank/, \
    $(addsuffix .prank.fasta, $(TEST_IDS)))

benchmark_fasta_aligned/prank: $(TSTALN_FILES_PRANK)

.PHONY: benchmark_fasta_aligned/prank

step/5_benchmark_alignments: benchmark_fasta_aligned/prank

benchmark_fasta_aligned/prank/prank.archive.tar.gz: benchmark_fasta/test_id_list.txt $(TSTALN_FILES_PRANK)
	cat $< | awk '{ print "benchmark_fasta_aligned/prank/" $$1 ".prank.fasta" }' | tar cvzf $@ -T -

step/5_benchmark_alignments: benchmark_fasta_aligned/prank/prank.archive.tar.gz

## MAFFT #######################################################################

TSTALN_FILES_MAFFT=$(addprefix benchmark_fasta_aligned/mafft/, \
    $(addsuffix .mafft.fasta, $(TEST_IDS)))

benchmark_fasta_aligned/mafft: $(TSTALN_FILES_MAFFT)

.PHONY: benchmark_fasta_aligned/mafft

step/5_benchmark_alignments: benchmark_fasta_aligned/mafft

benchmark_fasta_aligned/mafft/mafft.archive.tar.gz: benchmark_fasta/test_id_list.txt $(TSTALN_FILES_MAFFT)
	cat $< | awk '{ print "benchmark_fasta_aligned/mafft/" $$1 ".mafft.fasta" }' | tar cvzf $@ -T -

step/5_benchmark_alignments: benchmark_fasta_aligned/mafft/mafft.archive.tar.gz

## CLUSTAL OMEGA ################################################################

TSTALN_FILES_CLUSTALO=$(addprefix benchmark_fasta_aligned/clustalo/, \
    $(addsuffix .clustalo.fasta, $(TEST_IDS)))

benchmark_fasta_aligned/clustalo: $(TSTALN_FILES_CLUSTALO)

.PHONY: benchmark_fasta_aligned/clustalo

step/5_benchmark_alignments: benchmark_fasta_aligned/clustalo

benchmark_fasta_aligned/clustalo/clustalo.archive.tar.gz: benchmark_fasta/test_id_list.txt $(TSTALN_FILES_CLUSTALO)
	cat $< | awk '{ print "benchmark_fasta_aligned/clustalo/" $$1 ".clustalo.fasta" }' | tar cvzf $@ -T -

step/5_benchmark_alignments: benchmark_fasta_aligned/clustalo/clustalo.archive.tar.gz

## MACSE #######################################################################

TSTALN_FILES_MACSE=$(addprefix benchmark_fasta_aligned/macse/, \
    $(addsuffix .macse.fasta, $(TEST_IDS)))

benchmark_fasta_aligned/macse: $(TSTALN_FILES_MACSE)

.PHONY: benchmark_fasta_aligned/macse

step/5_benchmark_alignments: benchmark_fasta_aligned/macse

benchmark_fasta_aligned/macse/macse.archive.tar.gz: benchmark_fasta/test_id_list.txt $(TSTALN_FILES_MACSE)
	cat $< | awk '{ print "benchmark_fasta_aligned/macse/" $$1 ".macse.fasta" }' | tar cvzf $@ -T -

step/5_benchmark_alignments: benchmark_fasta_aligned/macse/macse.archive.tar.gz

################################################################################
# STEP 5a (Alt): Extract benchmark alignments from archives.                   #
################################################################################

step/5a_extract_benchmark_alignments:
	tar -xvmf benchmark_fasta_aligned/clustalo/clustalo.archive.tar.gz
	tar -xvmf benchmark_fasta_aligned/coati-dna/coati-dna.archive.tar.gz
	tar -xvmf benchmark_fasta_aligned/coati-mar-ecm/coati-mar-ecm.archive.tar.gz
	tar -xvmf benchmark_fasta_aligned/coati-mar-mg/coati-mar-mg.archive.tar.gz
	tar -xvmf benchmark_fasta_aligned/coati-tri-ecm/coati-tri-ecm.archive.tar.gz
	tar -xvmf benchmark_fasta_aligned/coati-tri-mg/coati-tri-mg.archive.tar.gz
	tar -xvmf benchmark_fasta_aligned/macse/macse.archive.tar.gz
	tar -xvmf benchmark_fasta_aligned/mafft/mafft.archive.tar.gz
	tar -xvmf benchmark_fasta_aligned/prank/prank.archive.tar.gz
	tar -xvmf benchmark_fasta_aligned/rev-coati-tri-mg/rev-coati-tri-mg.archive.tar.gz

.PHONY: step/5a_extract_benchmark_alignments

################################################################################
# Clean pipeline results except gene id list and raw fasta downloads		   #
################################################################################

.PHONY: clean_pipeline

clean_pipeline:
	rm -f data/{filtered,gaps,gaps_cigar,nogaps}.csv
	rm -f data/no_gaps_ref/*.fasta
	rm -f data/ref_alignments/*.fasta
	rm -f data/ref_alignments.csv
	rm -f data/cigar.rda
	rm -f aln/{tri-mg,mar-mg,dna,tri-ecm,mar-ecm,prank,mafft,clustalo,macse}/*.fasta
	rm -f aln/ref/{tri-mg,mar-mg,dna,tri-ecm,mar-ecm,prank,mafft,clustalo,macse}/*.fasta
	rm -f results/results_summary.csv

# find raw_fasta_aligned -type f -name '*.fasta' -exec touch {} +