################################################################################
# Global variables                                                             #
################################################################################
RSCRIPT = Rscript --vanilla
SHELL = /bin/bash

SPECIES ?= gorilla
RAW_PATH = raw_data/$(SPECIES)
N ?= 4000

################################################################################
# Download gene ID table from Ensembl (use to update files)                    #
################################################################################

raw_data/%_geneId.tsv: scripts/get_geneId.R
	@echo "Downloading" $* "gene IDs"
	@$(RSCRIPT) $< $* $@

################################################################################
# Download first N sequences from ENSEMBL                                      #
################################################################################
GENES = $(addsuffix .fasta,$(shell head -n${N} ${RAW_PATH}_geneId.tsv | cut -f1))
DOWNLOAD_GENES = $(addsuffix .fasta,$(addprefix $(RAW_PATH)/,$(shell head -n${N} ${RAW_PATH}_geneId.tsv | cut -f1)))

.PHONY: download_genes
download_genes: raw_data/$(SPECIES)_geneId.tsv $(DOWNLOAD_GENES)

$(RAW_PATH)/%.fasta: scripts/get_sequences.R | raw_data/$(SPECIES)_geneId.tsv
	@$(shell $(RSCRIPT) $^ raw_data/${SPECIES}_geneId.tsv $* $(@D))
	@echo -ne "Downloaded "$(shell ls ${RAW_PATH}/*.fasta | wc -l)" out of ${N}.\r"

################################################################################
# Initial alignments with all methods (aln recipes in Makefile_aln.mak)        #
################################################################################
INITIAL_ALIGNMENTS = $(DOWNLOAD_GENES) \
					$(addprefix aln/mecm/,$(GENES)) \
					$(addprefix aln/macse/,$(GENES)) \
					$(addprefix aln/mafft/,$(GENES)) \
 					$(addprefix aln/mcoati/,$(GENES)) \
					$(addprefix aln/clustalo/,$(GENES))\
 					$(addprefix aln/prank/,$(GENES))\
					$(addprefix aln/coati/,$(GENES)) \
					$(addprefix aln/dna/,$(GENES)) \
					$(addprefix aln/ecm/,$(GENES)) \

################################################################################
# Identify which initial alignments have gaps                                  #
################################################################################
MODELS = coati mcoati dna ecm mecm prank mafft clustalo macse
GAPS_FILE = data/$(SPECIES)/gaps.csv
$(GAPS_FILE): scripts/gaps.sh | $(INITIAL_ALIGNMENTS)
	@echo "Find alignments with gaps"
	@$(shell bash $< ${SPECIES} ${MODELS})
	@$(eval GAPS=$(addsuffix .gap,$(shell cat data/$(SPECIES)/gaps.csv)))


################################################################################
# Create cigar strings with indel types for simulation                         #
################################################################################
# GAPS = $(addsuffix .gap,$(shell cat $(GAPS_FILE)))

data/$(SPECIES)/gaps_cigar.csv: scripts/gaps2cigar.R $(GAPS_FILE)
	@$(eval GAPS=$(addsuffix .gap,$(shell cat data/$(SPECIES)/gaps.csv))); \
	$(MAKE) $(GAPS)
	@echo "Encoded gap patterns using CIGAR strings"

# extract gap information and encode it using CIGAR strings
%.gap: scripts/gaps2cigar.R
	@$(shell ${RSCRIPT} $< $* | cut -d '"' -f 2 >> $(DATA)/gaps_cigar.csv)

################################################################################
# Create file with un-gaped initial alignments
################################################################################

data/$(SPECIES)/nogaps.csv: scripts/nogaps.sh $(GAPS_FILE)
	@echo "Find alignments without gaps"
	@bash $< ${SPECIES} ${N}
	@$(eval REF_ALIG=$(addprefix ${DATA}/ref_alignments/,$(addsuffix .fasta,$(shell cat ${DATA}/nogaps.csv))))

################################################################################
# Generate alignments with biologically inferred indel information (cigar str) #
################################################################################
DATA = data/$(SPECIES)
# REF_ALIG = $(addprefix ${DATA}/ref_alignments/,$(addsuffix .fasta,$(shell cat ${DATA}/nogaps.csv)))

.PHONY: reference
reference: $(DATA)/gaps_cigar.csv $(DATA)/nogaps.csv
	@$(eval REF_ALIG=$(addprefix ${DATA}/ref_alignments/,$(addsuffix .fasta,$(shell cat ${DATA}/nogaps.csv)))); \
	$(MAKE) $(REF_ALIG)
	@echo "Created reference alignments"

data/$(SPECIES)/ref_alignments/%: scripts/simulate2.R scripts/write_fasta.R $(DATA)/gaps_cigar.csv
# 	@echo "Creating reference alignment $*"
	@echo -ne "Creating reference alignment $*\r"
	@$(shell ${RSCRIPT} $< ${SPECIES} ${RAW_PATH}/$* $@)

################################################################################
# Create reference alignments with no gaps for testing                         #
################################################################################

.PHONY: no_gaps_reference
no_gaps_reference: reference
	$(shell bash scripts/rem_gaps_ref.sh $(SPECIES))
	@$(eval ALN=$(shell ls ${DATA}/no_gaps_ref))

################################################################################
# Align simulated alignments with all methods (aln recipes in Makefile_aln.mak)#
################################################################################
REF_PATH = $(DATA)/no_gaps_ref
ALN = $(shell ls $(REF_PATH)/)

ALIGN_REFERENCE = no_gaps_reference \
				$(addprefix aln/ref/mecm/,$(ALN)) \
				$(addprefix aln/ref/macse/,$(ALN)) \
				$(addprefix aln/ref/mafft/,$(ALN)) \
 				$(addprefix aln/ref/mcoati/,$(ALN)) \
				$(addprefix aln/ref/clustalo/,$(ALN))\
 				$(addprefix aln/ref/prank/,$(ALN))\
				$(addprefix aln/ref/coati/,$(ALN)) \
				$(addprefix aln/ref/dna/,$(ALN)) \
				$(addprefix aln/ref/ecm/,$(ALN)) 

.PHONY: align_ref
align_ref: $(ALIGN_REFERENCE)

################################################################################
# Summary statistics                                                           #
################################################################################
ALN_DSEQ = $(addprefix $(DATA)/, $(addsuffix .dseq, $(ALN)))
REF = aln/ref

$(DATA)/dseq_summary.csv: scripts/dseq_summary.R $(DATA)/dseq.csv
	$(shell Rscript --vanilla $^ $@)
	@cat $@

$(DATA)/dseq.csv: $(ALN_DSEQ) $(ALIGN_REFERENCE)
	$(shell sed -i '1 i\filename,mecm,macse,mafft,mcoati,clustalo,prank,coati,dna,ecm' $@)

$(DATA)/%.dseq: scripts/pa_distance.R 
	@echo "dseq $*"
	$(shell $(RSCRIPT) $< $(DATA)/ref_alignments/$* mecm/$* macse/$* mafft/$* \
						mcoati/$* clustalo/$* prank/$* coati/$* dna/$* ecm/$* \
						1 >> $(DATA)/dseq.csv)
	$(shell printf "\n" >> $(DATA)/dseq.csv)

################################################################################
# Clean pipeline results except gene id list and raw fasta downloads		   #
################################################################################

.PHONY: clean_pipeline

clean_pipeline:
	rm -f data/$(SPECIES)/*.csv
	rm -f data/$(SPECIES)/no_gaps_ref/*
	rm -f data/$(SPECIES)/ref_alignments/*
	rm -f aln/{coati,mcoati,dna,ecm,mecm,prank,mafft,clustalo,macse,baliphy}/*
	rm -f aln/ref/{coati,mcoati,dna,ecm,mecm,prank,mafft,clustalo,macse,baliphy}/*
	rm -f data/dseq.csv
	rm -f data/dseq_summary.csv

include Makefile_aln.mak
# include Makefile_stats.mak
