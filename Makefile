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

raw_fasta/gene_id_list.mk: raw_fasta/hs-gg_gene_pairs.csv.gz
	echo "GENE_IDS = \\" > $@
	zcat $< | awk -F, 'NR > 1 { print "   " $$1 " \\" }' >> $@
	echo "" >> $@

results/raw_fasta_metrics.csv: scripts/create_raw_fasta_metrics.R raw_fasta/.script_done
	$(RSCRIPT) scripts/create_raw_fasta_metrics.R raw_fasta $@

raw_fasta: raw_fasta/.script_done raw_fasta/gene_id_list.mk results/raw_fasta_metrics.csv

.PHONY: raw_fasta

step/1_download_raw_fasta: raw_fasta

################################################################################
# STEP 2: Align empirical CDS sequences                                        #
################################################################################

COATI_BIN = ./bin/coati-alignpair
PRANK_BIN = ./bin/prank
MAFFT_BIN = ./bin/mafft
CLUSTALO_BIN = ./bin/clustalo
MACSE_JAR = ./bin/macse_v2.06.jar

step/2_empirical_alignments:

.PHONY: step/2_empirical_alignments

## COATI TRI-MG ################################################################

raw_fasta_aligned/coati-tri-mg/%.coati-tri-mg.fasta: raw_fasta/%.fasta
	$(COATI_BIN) $< -m tri-mg -o $@

raw_fasta_aligned/coati-tri-mg: $(addprefix raw_fasta_aligned/coati-tri-mg/, $(addsuffix .coati-tri-mg.fasta, $(GENE_IDS)))

.PHONY: raw_fasta_aligned/coati-tri-mg

step/2_empirical_alignments: raw_fasta_aligned/coati-tri-mg

## COATI TRI-ECM ###############################################################

raw_fasta_aligned/coati-tri-ecm/%.coati-tri-ecm.fasta: raw_fasta/%.fasta
	$(COATI_BIN) $< -m tri-ecm -o $@

raw_fasta_aligned/coati-tri-ecm: $(addprefix raw_fasta_aligned/coati-tri-ecm/, $(addsuffix .coati-tri-ecm.fasta, $(GENE_IDS)))

.PHONY: raw_fasta_aligned/coati-tri-ecm

step/2_empirical_alignments: raw_fasta_aligned/coati-tri-ecm

## COATI MAR-MG ################################################################

raw_fasta_aligned/coati-mar-mg/%.coati-mar-mg.fasta: raw_fasta/%.fasta
	$(COATI_BIN) $< -m mar-mg -o $@

raw_fasta_aligned/coati-mar-mg: $(addprefix raw_fasta_aligned/coati-mar-mg/, $(addsuffix .coati-mar-mg.fasta, $(GENE_IDS)))

.PHONY: raw_fasta_aligned/coati-mar-mg

step/2_empirical_alignments: raw_fasta_aligned/coati-mar-mg

## COATI MAR-ECM ###############################################################

raw_fasta_aligned/coati-mar-ecm/%.coati-mar-ecm.fasta: raw_fasta/%.fasta
	$(COATI_BIN) $< -m mar-ecm -o $@

raw_fasta_aligned/coati-mar-ecm: $(addprefix raw_fasta_aligned/coati-mar-ecm/, $(addsuffix .coati-mar-ecm.fasta, $(GENE_IDS)))

.PHONY: raw_fasta_aligned/coati-mar-ecm

step/2_empirical_alignments: raw_fasta_aligned/coati-mar-ecm

## COATI DNA ###################################################################

raw_fasta_aligned/coati-dna/%.coati-dna.fasta: raw_fasta/%.fasta
	$(COATI_BIN) $< -m dna -o $@

raw_fasta_aligned/coati-dna: $(addprefix raw_fasta_aligned/coati-dna/, $(addsuffix .coati-dna.fasta, $(GENE_IDS)))

.PHONY: raw_fasta_aligned/coati-dna

step/2_empirical_alignments: raw_fasta_aligned/coati-dna

## REV COATI TRI-MG ############################################################

raw_fasta_aligned/rev-coati-tri-mg/%.rev-coati-tri-mg.fasta: raw_fasta/%.fasta
	$(RSCRIPT) scripts/rev_coati_wrapper.R "$(COATI_BIN)" tri-mg $< $@

raw_fasta_aligned/rev-coati-tri-mg: $(addprefix raw_fasta_aligned/rev-coati-tri-mg/, $(addsuffix .rev-coati-tri-mg.fasta, $(GENE_IDS)))

.PHONY: raw_fasta_aligned/rev-coati-tri-mg

step/2_empirical_alignments: raw_fasta_aligned/rev-coati-tri-mg

## PRANK #######################################################################

raw_fasta_aligned/prank/%.prank.fasta: raw_fasta/%.fasta
	$(PRANK_BIN) -codon -d="$<" -o=$@ -quiet &>/dev/null
	if [ -f $@.best.fas ]; then mv $@.best.fas $@; else touch $@; fi

raw_fasta_aligned/prank: $(addprefix raw_fasta_aligned/prank/, $(addsuffix .prank.fasta, $(GENE_IDS)))

.PHONY: raw_fasta_aligned/prank

step/2_empirical_alignments: raw_fasta_aligned/prank

## MAFFT #######################################################################

raw_fasta_aligned/mafft/%.mafft.fasta: raw_fasta/%.fasta
	$(MAFFT_BIN) --nuc --globalpair --maxiterate 1000 --preservecase --quiet $< > $@
	
raw_fasta_aligned/mafft: $(addprefix raw_fasta_aligned/mafft/, $(addsuffix .mafft.fasta, $(GENE_IDS)))

.PHONY: raw_fasta_aligned/mafft

step/2_empirical_alignments: raw_fasta_aligned/mafft

## CLUSTAL OMEGA ################################################################

# Wrap the following command:
#    clustalo --seqtype=Protein --infile=- --output-order=input-order

raw_fasta_aligned/clustalo/%.clustalo.fasta: raw_fasta/%.fasta
	$(RSCRIPT) scripts/clustalo_wrapper.R "$(CLUSTALO_BIN)" $< $@
	
raw_fasta_aligned/clustalo: $(addprefix raw_fasta_aligned/clustalo/, $(addsuffix .clustalo.fasta, $(GENE_IDS)))

.PHONY: raw_fasta_aligned/clustalo

step/2_empirical_alignments: raw_fasta_aligned/clustalo

## MACSE #######################################################################

# Wrap the following command:
#    java -jar macse.jar -prog alignSequences -seq temp1 -seq_lr temp2

raw_fasta_aligned/macse/%.macse.fasta: raw_fasta/%.fasta
	$(RSCRIPT) scripts/macse_wrapper.R "$(MACSE_JAR)" $< $@
	
raw_fasta_aligned/macse: $(addprefix raw_fasta_aligned/macse/, $(addsuffix .macse.fasta, $(GENE_IDS)))

.PHONY: raw_fasta_aligned/macse

step/2_empirical_alignments: raw_fasta_aligned/macse

################################################################################
# Identify which initial alignments have gaps                                  #
################################################################################
MODELS = $(COATI_MODEL) prank mafft clustalo macse
GAPS_FILE = data/gaps.csv

$(GAPS_FILE): scripts/gaps.sh
	@echo "Find alignments with gaps               "
	@bash $< ${MODELS}

################################################################################
# Create cigar strings with indel types for simulation                         #
################################################################################
GAPS = $(addsuffix .gap,$(shell cat $(GAPS_FILE) 2> /dev/null))

data/gaps_cigar.csv: scripts/gaps2cigar.R $(GAPS_FILE) $(GAPS)
	@echo "Encoded gap patterns using CIGAR strings              "
	@sed -i '1 i\raw_name,cigar,origin' $@

# extract gap information and encode it using CIGAR strings
%.gap: scripts/gaps2cigar.R
	@echo -ne "Encoding gaps from " $* "     \r"
	@${RSCRIPT} $< $* | cut -d '"' -f 2 >> data/gaps_cigar.csv

################################################################################
# Create file with un-gaped initial alignments
################################################################################

data/nogaps.csv: scripts/nogaps.sh $(GAPS_FILE)
	@echo "Find alignments without gaps                       "
	@bash $<

################################################################################
# Generate alignments with biologically inferred indel information (cigar str) #
################################################################################
REF_ALIG = $(addprefix data/ref_alignments/,$(shell cat data/nogaps.csv 2> /dev/null))

.PHONY: reference
reference: $(REF_ALIG) | data/nogaps.csv # data/gaps_cigar.csv
	@echo "Done creating reference alignments                   "
	@sed -i '1 i\raw_name,cigar,origin,ref_name' data/ref_alignments.csv

data/ref_alignments/%: scripts/simulate2.R data/cigar.rda scripts/write_fasta.R
	@echo -ne "Creating reference alignment $*\r"
	@timeout 20s ${RSCRIPT} $< ${RAW_PATH}/$* $@ | cut -d '"' -f 2 >> data/ref_alignments.csv

data/cigar.rda: scripts/create_cigar_list.R
	@${RSCRIPT} $<

################################################################################
# Create reference alignments with no gaps for testing                         #
################################################################################

.PHONY: no_gaps_reference
no_gaps_reference: data/nogaps.csv# data/gaps_cigar.csv
	@bash scripts/rem_gaps_ref.sh

################################################################################
# Align simulated alignments with all methods (aln recipes in Makefile_aln.mak)#
################################################################################
REF_PATH = data/no_gaps_ref
ALN = $(shell ls $(REF_PATH)/ )

ALIGN_REFERENCE = no_gaps_reference \
				$(addprefix aln/ref/macse/,$(ALN)) \
				$(addprefix aln/ref/mafft/,$(ALN)) \
				$(addprefix aln/ref/$(COATI_MODEL)/,$(ALN)) \
				$(addprefix aln/ref/clustalo/,$(ALN))\
				$(addprefix aln/ref/prank/,$(ALN))

.PHONY: align_ref
align_ref: $(ALIGN_REFERENCE)

################################################################################
# Summary statistics                                                           #
################################################################################
ALL_ALNS = $(shell ls data/ref_alignments/)
RES = $(addprefix results/, $(addsuffix .res, $(ALL_ALNS)))

results/results_summary.csv: $(RES)
	@sed -i '1 i\ref_name,raw_name,cigar,origin,aligner,dseq,dpos,ref_omega,aln_omega,ref_score,score' $@
	@echo -ne "Built $@   \n"

results/%.res: scripts/results_summary.R
	@echo -ne "summary stats $*\r"
	@$(RSCRIPT) $< $* ${MODELS}

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

include Makefile_aln.mak
