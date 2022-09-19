################################################################################
# Global variables                                                             #
################################################################################
RSCRIPT = Rscript --vanilla
SHELL = /bin/bash

SPECIES ?= gorilla
RAW_PATH = raw_data/$(SPECIES)
N ?= 4000
LEN ?= 7000

################################################################################
# Download gene ID table from Ensembl (use to update files)                    #
################################################################################

.PHONY: download_geneid
download_geneid: | raw_data/$(SPECIES)_geneId.tsv

raw_data/%_geneId.tsv: scripts/get_geneId.R
	@echo "Downloading" $* "gene IDs"
	@$(RSCRIPT) $< $* $@

################################################################################
# Download first N sequences from ENSEMBL                                      #
################################################################################
GENES = $(addsuffix .fasta,$(shell head -n${N} ${RAW_PATH}_geneId.tsv | cut -f1))
DOWNLOAD_GENES = $(addprefix $(RAW_PATH)/,$(GENES))

.PHONY: download_genes
download_genes: $(DOWNLOAD_GENES)

$(RAW_PATH)/%.fasta: | scripts/get_sequences.R raw_data/$(SPECIES)_geneId.tsv
	@$(shell $(RSCRIPT) scripts/get_sequences.R raw_data/${SPECIES}_geneId.tsv $* $(@D))
	@echo -ne "Downloaded " $* "("$(shell ls ${RAW_PATH}/*.fasta | wc -l)" out of ${N}).\r"

################################################################################
# Filter sequences by length												   #
################################################################################

.PHONY: filter
filter: data/$(SPECIES)/filtered.csv

data/$(SPECIES)/filtered.csv: $(DOWNLOAD_GENES) scripts/filter_seqs.sh
	@bash scripts/filter_seqs.sh $(SPECIES) $(N) $(LEN) $@

################################################################################
# Initial alignments with all methods (aln recipes in Makefile_aln.mak)        #
################################################################################
.PHONY: initial_alignment

FILTERED = $(shell cat data/${SPECIES}/filtered.csv 2> /dev/null)
INITIAL_ALIGNMENTS = $(DOWNLOAD_GENES) \
					$(addprefix aln/mecm/,$(FILTERED)) \
					$(addprefix aln/macse/,$(FILTERED)) \
					$(addprefix aln/mafft/,$(FILTERED)) \
					$(addprefix aln/mcoati/,$(FILTERED)) \
					$(addprefix aln/clustalo/,$(FILTERED)) \
					$(addprefix aln/prank/,$(FILTERED)) \
					$(addprefix aln/coati/,$(FILTERED)) \
					$(addprefix aln/dna/,$(FILTERED)) \
					$(addprefix aln/ecm/,$(FILTERED)) \

initial_alignment: $(INITIAL_ALIGNMENTS)

################################################################################
# Identify which initial alignments have gaps                                  #
################################################################################
MODELS = coati mcoati dna ecm mecm prank mafft clustalo macse
GAPS_FILE = data/$(SPECIES)/gaps.csv

$(GAPS_FILE): scripts/gaps.sh | $(INITIAL_ALIGNMENTS)
	@echo "Find alignments with gaps               "
	@$(shell bash $< ${SPECIES} ${MODELS})

################################################################################
# Create cigar strings with indel types for simulation                         #
################################################################################
GAPS = $(addsuffix .gap,$(shell cat $(GAPS_FILE) 2> /dev/null))

data/$(SPECIES)/gaps_cigar.csv: scripts/gaps2cigar.R | $(GAPS_FILE) $(GAPS)
	@echo "Encoded gap patterns using CIGAR strings"

# extract gap information and encode it using CIGAR strings
%.gap: scripts/gaps2cigar.R
	@echo -ne "Encoding gaps from " $* "     \r"
	@$(shell ${RSCRIPT} $< $* | cut -d '"' -f 2 >> $(DATA)/gaps_cigar.csv)

################################################################################
# Create file with un-gaped initial alignments
################################################################################

data/$(SPECIES)/nogaps.csv: scripts/nogaps.sh $(GAPS_FILE)
	@echo "Find alignments without gaps"
	@bash $< ${SPECIES}
	@mkdir -p data/$(SPECIES)/ref_alignments

################################################################################
# Generate alignments with biologically inferred indel information (cigar str) #
################################################################################
DATA = data/$(SPECIES)
REF_ALIG = $(addprefix ${DATA}/ref_alignments/,$(shell cat ${DATA}/nogaps.csv 2> /dev/null))
$(mkdir -p data/$(SPECIES)/ref_alignments)

.PHONY: reference
reference: $(REF_ALIG) | $(DATA)/nogaps.csv
	@echo "Done creating reference alignments"

$(DATA)/ref_alignments/%: scripts/simulate2.R scripts/write_fasta.R
	@echo -ne "Creating reference alignment $*\r"
	@$(shell ${RSCRIPT} $< ${SPECIES} ${RAW_PATH}/$* $@)

################################################################################
# Create reference alignments with no gaps for testing                         #
################################################################################
$(mkdir -p data/$(SPECIES)/no_gaps_ref)

.PHONY: no_gaps_reference
no_gaps_reference: $(REF_ALIG)
	$(shell bash scripts/rem_gaps_ref.sh $(SPECIES))

################################################################################
# Align simulated alignments with all methods (aln recipes in Makefile_aln.mak)#
################################################################################
REF_PATH = $(DATA)/no_gaps_ref
ALN = $(shell mkdir -p $(REF_PATH); ls $(REF_PATH)/ )

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

$(DATA)/dseq.csv: $(ALN_DSEQ)
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
	rm -f aln/{coati,mcoati,dna,ecm,mecm,prank,mafft,clustalo,macse}/*
	rm -f aln/ref/{coati,mcoati,dna,ecm,mecm,prank,mafft,clustalo,macse}/*
	rm -f data/dseq.csv
	rm -f data/dseq_summary.csv

include Makefile_aln.mak
# include Makefile_stats.mak
