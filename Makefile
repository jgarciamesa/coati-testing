################################################################################
# Global variables                                                             #
################################################################################
RSCRIPT = Rscript --vanilla
SHELL = /bin/bash

RAW_PATH = raw_data/
N ?= 16000
LEN ?= 6000

################################################################################
# Download gene ID table from Ensembl (use to update files)                    #
################################################################################

.PHONY: download_geneid
download_geneid: | raw_data/gorilla_geneId.tsv

raw_data/gorilla_geneId.tsv: scripts/get_geneId.R
	@echo "Downloading" $* "gene IDs"
	@$(RSCRIPT) $< $* $@

################################################################################
# Download first N sequences from ENSEMBL                                      #
################################################################################
GENES = $(addsuffix .fasta,$(shell head -n${N} ${RAW_PATH}/gorilla_geneId.tsv | cut -f1))
DOWNLOAD_GENES = $(addprefix $(RAW_PATH)/,$(GENES))

.PHONY: download_genes
download_genes: $(DOWNLOAD_GENES)

$(RAW_PATH)/%.fasta: | scripts/get_sequences.R
	@$(shell $(RSCRIPT) scripts/get_sequences.R raw_data/gorilla_geneId.tsv $* $(@D))
	@echo -ne "Downloaded " $* "("$(shell ls ${RAW_PATH}/*.fasta | wc -l)" out of ${N}).\r"

################################################################################
# Filter sequences by length												   #
################################################################################

.PHONY: filter
filter: data/filtered.csv

data/filtered.csv: $(DOWNLOAD_GENES) scripts/filter_seqs.sh
	@bash scripts/filter_seqs.sh $(N) $(LEN) $@

################################################################################
# Initial alignments with all methods (aln recipes in Makefile_aln.mak)        #
################################################################################
.PHONY: initial_alignment

FILTERED = $(shell cat data/filtered.csv 2> /dev/null)
INITIAL_ALIGNMENTS = data/filtered.csv \
					$(addprefix aln/macse/,$(FILTERED)) \
					$(addprefix aln/mafft/,$(FILTERED)) \
					$(addprefix aln/ecm/,$(FILTERED)) \
					$(addprefix aln/clustalo/,$(FILTERED)) \
					$(addprefix aln/prank/,$(FILTERED))

initial_alignment: $(INITIAL_ALIGNMENTS)

################################################################################
# Identify which initial alignments have gaps                                  #
################################################################################
MODELS = ecm prank mafft clustalo macse
GAPS_FILE = data/gaps.csv

$(GAPS_FILE): scripts/gaps.sh
	@echo "Find alignments with gaps               "
	@$(shell bash $< ${MODELS})

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
	@$(shell ${RSCRIPT} $< $* | cut -d '"' -f 2 >> data/gaps_cigar.csv)

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

data/ref_alignments/%: scripts/simulate2.R scripts/write_fasta.R
	@echo -ne "Creating reference alignment $*\r"
	@$(shell timeout 60s ${RSCRIPT} $< ${RAW_PATH}/$* $@ | cut -d '"' -f 2 >> data/ref_alignments.csv)

################################################################################
# Create reference alignments with no gaps for testing                         #
################################################################################

.PHONY: no_gaps_reference
no_gaps_reference: $(REF_ALIG) data/nogaps.csv# data/gaps_cigar.csv
	@bash scripts/rem_gaps_ref.sh

################################################################################
# Align simulated alignments with all methods (aln recipes in Makefile_aln.mak)#
################################################################################
REF_PATH = data/no_gaps_ref
ALN = $(shell ls $(REF_PATH)/ )

ALIGN_REFERENCE = no_gaps_reference \
				$(addprefix aln/ref/macse/,$(ALN)) \
				$(addprefix aln/ref/mafft/,$(ALN)) \
				$(addprefix aln/ref/ecm/,$(ALN)) \
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
	@sed -i '1 i\ref_name,raw_name,cigar,origin,aligner,dseq,ref_omega,aln_omega' $@
	@echo -ne "Built $@   \n"

results/%.res: scripts/results_summary.R
	@echo -ne "summary stats $*\r"
	$(shell $(RSCRIPT) $< $* ${MODELS})

supplementary_materials.pdf: supplementary_materials.Rmd
	@Rscript -e "rmarkdown::render('supplementary_materials.Rmd')"

################################################################################
# Clean pipeline results except gene id list and raw fasta downloads		   #
################################################################################

.PHONY: clean_pipeline

clean_pipeline:
	rm -f data/{filtered,gaps,gaps_cigar,nogaps}.csv
	rm -f data/no_gaps_ref/*.fasta
	rm -f data/ref_alignments/*.fasta
	rm -f data/ref_alignments.csv
	rm -f aln/{coati,mcoati,dna,ecm,mecm,prank,mafft,clustalo,macse}/*.fasta
	rm -f aln/ref/{coati,mcoati,dna,ecm,mecm,prank,mafft,clustalo,macse}/*.fasta
	rm -f results/results_summary.csv

include Makefile_aln.mak
