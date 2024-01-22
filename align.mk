# This file contains the rules for aligning sequences. This allows alignment to
# be parallelized without having to process the full dependency list.

RSCRIPT ?= Rscript --vanilla
SHELL ?= /bin/bash -e -o pipefail

COATI_BIN = ./bin/coati-alignpair
PRANK_BIN = ./bin/prank
MAFFT_BIN = ./bin/mafft
CLUSTALO_BIN = ./bin/clustalo
MACSE_JAR = ./bin/macse_v2.06.jar

## COATI TRI-MG ################################################################

raw_fasta_aligned/coati-tri-mg/%.coati-tri-mg.fasta: raw_fasta/%.fasta
	$(COATI_BIN) $< -m tri-mg -o $@ || [ $$? -eq 1 ] && touch $@

benchmark_fasta_aligned/coati-tri-mg/%.coati-tri-mg.fasta: benchmark_fasta_nogaps/%.nogaps.fasta
	$(COATI_BIN) $< -m tri-mg -o $@ || [ $$? -eq 1 ] && touch $@

## COATI TRI-ECM ###############################################################

raw_fasta_aligned/coati-tri-ecm/%.coati-tri-ecm.fasta: raw_fasta/%.fasta
	$(COATI_BIN) $< -m tri-ecm -o $@ || [ $$? -eq 1 ] && touch $@

benchmark_fasta_aligned/coati-tri-ecm/%.coati-tri-ecm.fasta: benchmark_fasta_nogaps/%.nogaps.fasta
	$(COATI_BIN) $< -m tri-ecm -o $@ || [ $$? -eq 1 ] && touch $@

## COATI MAR-MG ################################################################

raw_fasta_aligned/coati-mar-mg/%.coati-mar-mg.fasta: raw_fasta/%.fasta
	$(COATI_BIN) $< -m mar-mg -o $@ || [ $$? -eq 1 ] && touch $@

benchmark_fasta_aligned/coati-mar-mg/%.coati-mar-mg.fasta: benchmark_fasta_nogaps/%.nogaps.fasta
	$(COATI_BIN) $< -m mar-mg -o $@ || [ $$? -eq 1 ] && touch $@

## COATI MAR-ECM ###############################################################

raw_fasta_aligned/coati-mar-ecm/%.coati-mar-ecm.fasta: raw_fasta/%.fasta
	$(COATI_BIN) $< -m mar-ecm -o $@ || [ $$? -eq 1 ] && touch $@

benchmark_fasta_aligned/coati-mar-ecm/%.coati-mar-ecm.fasta: benchmark_fasta_nogaps/%.nogaps.fasta
	$(COATI_BIN) $< -m mar-ecm -o $@ || [ $$? -eq 1 ] && touch $@

## COATI DNA ###################################################################

raw_fasta_aligned/coati-dna/%.coati-dna.fasta: raw_fasta/%.fasta
	$(COATI_BIN) $< -m dna -o $@ || [ $$? -eq 1 ] && touch $@

benchmark_fasta_aligned/coati-dna/%.coati-dna.fasta: benchmark_fasta_nogaps/%.nogaps.fasta
	$(COATI_BIN) $< -m dna -o $@ || [ $$? -eq 1 ] && touch $@

## REV COATI TRI-MG ############################################################

raw_fasta_aligned/rev-coati-tri-mg/%.rev-coati-tri-mg.fasta: raw_fasta/%.fasta
	$(RSCRIPT) scripts/rev_coati_wrapper.R "$(COATI_BIN)" tri-mg $< $@ || [ $$? -eq 1 ] && touch $@

benchmark_fasta_aligned/rev-coati-tri-mg/%.rev-coati-tri-mg.fasta: benchmark_fasta_nogaps/%.nogaps.fasta
	$(RSCRIPT) scripts/rev_coati_wrapper.R "$(COATI_BIN)" tri-mg $< $@ || [ $$? -eq 1 ] && touch $@

## PRANK #######################################################################

raw_fasta_aligned/prank/%.prank.fasta: raw_fasta/%.fasta
	$(PRANK_BIN) -codon -d="$<" -o=$@ -quiet &>/dev/null
	if [ -f $@.best.fas ]; then mv $@.best.fas $@; else touch $@; fi

benchmark_fasta_aligned/prank/%.prank.fasta: benchmark_fasta_nogaps/%.nogaps.fasta
	$(PRANK_BIN) -codon -d="$<" -o=$@ -quiet &>/dev/null
	if [ -f $@.best.fas ]; then mv $@.best.fas $@; else touch $@; fi

## MAFFT #######################################################################

raw_fasta_aligned/mafft/%.mafft.fasta: raw_fasta/%.fasta
	$(MAFFT_BIN) --nuc --globalpair --maxiterate 1000 --preservecase --quiet $< > $@

benchmark_fasta_aligned/mafft/%.mafft.fasta: benchmark_fasta_nogaps/%.nogaps.fasta
	$(MAFFT_BIN) --nuc --globalpair --maxiterate 1000 --preservecase --quiet $< > $@

## CLUSTAL OMEGA ################################################################

# Wrap the following command:
#    clustalo --seqtype=Protein --infile=- --output-order=input-order

raw_fasta_aligned/clustalo/%.clustalo.fasta: raw_fasta/%.fasta
	$(RSCRIPT) scripts/clustalo_wrapper.R "$(CLUSTALO_BIN)" $< $@

benchmark_fasta_aligned/clustalo/%.clustalo.fasta: benchmark_fasta_nogaps/%.nogaps.fasta
	$(RSCRIPT) scripts/clustalo_wrapper.R "$(CLUSTALO_BIN)" $< $@

## MACSE #######################################################################

# Wrap the following command:
#    java -jar macse.jar -prog alignSequences -seq temp1 -seq_lr temp2

raw_fasta_aligned/macse/%.macse.fasta: raw_fasta/%.fasta
	$(RSCRIPT) scripts/macse_wrapper.R "$(MACSE_JAR)" $< $@

benchmark_fasta_aligned/macse/%.macse.fasta: benchmark_fasta_nogaps/%.nogaps.fasta
	$(RSCRIPT) scripts/macse_wrapper.R "$(MACSE_JAR)" $< $@
