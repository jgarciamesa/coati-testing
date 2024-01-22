
default:

.PHONY: default

METHODS = clustalo coati-dna coati-mar-ecm coati-mar-mg coati-tri-ecm \
    coati-tri-mg macse mafft prank rev-coati-tri-mg

################################################################################
# STEP 2 (Alt): Empirical alignments archives.                                 #
################################################################################

raw_fasta_aligned/%/archive.tar.gz: raw_fasta/gene_id_list.txt
	awk '{ print "$(@D)/" $$1 ".$*.fasta" }' $< | tar cvzf $@ -T -

$(addprefix archive/2_empirical_alignments_,$(METHODS)): \
archive/2_empirical_alignments_%: raw_fasta_aligned/%/archive.tar.gz

$(addprefix extract/2_empirical_alignments_,$(METHODS)): \
extract/2_empirical_alignments_%: 
	tar -xvmf raw_fasta_aligned/$*/archive.tar.gz

$(addprefix clean/2_empirical_alignments_,$(METHODS)): \
clean/2_empirical_alignments_%: 
	rm -f raw_fasta_aligned/$*/archive.tar.gz

archive/2_empirical_alignments: $(addprefix archive/2_empirical_alignments_,$(METHODS))
extract/2_empirical_alignments: $(addprefix extract/2_empirical_alignments_,$(METHODS))
clean/2_empirical_alignments: $(addprefix clean/2_empirical_alignments_,$(METHODS))

.PHONY: archive/2_empirical_alignments
.PHONY: extract/2_empirical_alignments
.PHONY: clean/2_empirical_alignments
.PHONY: $(addprefix archive/2_empirical_alignments_,$(METHODS))
.PHONY: $(addprefix extract/2_empirical_alignments_,$(METHODS))
.PHONY: $(addprefix clean/2_empirical_alignments_,$(METHODS))

################################################################################
# STEP 4 (Alt): Benchmarks archives.                                           #
################################################################################

benchmark_fasta/archive.tar.gz: benchmark_fasta/test_id_list.txt
	awk '{ print "$(@D)/" $$1 ".fasta" }' $< | tar cvzf $@ -T -

benchmark_fasta_nogaps/archive.tar.gz: benchmark_fasta/test_id_list.txt
	awk '{ print "$(@D)/" $$1 ".nogaps.fasta" }' $< | tar cvzf $@ -T -

archive/4_benchmarks_fasta: benchmark_fasta/archive.tar.gz
archive/4_benchmarks_fasta_nogaps: benchmark_fasta_nogaps/archive.tar.gz

extract/4_benchmarks_fasta:
	tar -xvmf benchmark_fasta/archive.tar.gz

extract/4_benchmarks_fasta_nogaps:
	tar -xvmf benchmark_fasta_nogaps/archive.tar.gz

clean/4_benchmarks_fasta:
	rm -f benchmark_fasta/archive.tar.gz

clean/4_benchmarks_fasta_nogaps:
	rm -f benchmark_fasta_nogaps/archive.tar.gz

archive/4_benchmarks: archive/4_benchmarks_fasta archive/4_benchmarks_fasta_nogaps
extract/4_benchmarks: extract/4_benchmarks_fasta extract/4_benchmarks_fasta_nogaps
clean/4_benchmarks: clean/4_benchmarks_fasta clean/4_benchmarks_fasta_nogaps

.PHONY: archive/4_benchmarks
.PHONY: archive/4_benchmarks_fasta
.PHONY: archive/4_benchmarks_fasta_nogaps
.PHONY: extract/4_benchmarks
.PHONY: extract/4_benchmarks_fasta
.PHONY: extract/4_benchmarks_fasta_nogaps
.PHONY: clean/4_benchmarks
.PHONY: clean/4_benchmarks_fasta
.PHONY: clean/4_benchmarks_fasta_nogaps

################################################################################
# STEP 5 (Alt): Benchmark alignments archives.                                 #
################################################################################

benchmark_fasta_aligned/%/archive.tar.gz: benchmark_fasta/test_id_list.txt
	awk '{ print "$(@D)/" $$1 ".$*.fasta" }' $< | tar cvzf $@ -T -

$(addprefix archive/5_benchmark_alignments_,$(METHODS)): \
archive/5_benchmark_alignments_%: benchmark_fasta_aligned/%/archive.tar.gz

$(addprefix extract/5_benchmark_alignments_,$(METHODS)): \
extract/5_benchmark_alignments_%: 
	tar -xvmf benchmark_fasta_aligned/$*/archive.tar.gz

$(addprefix clean/5_benchmark_alignments_,$(METHODS)): \
clean/5_benchmark_alignments_%: 
	rm -f benchmark_fasta_aligned/$*/archive.tar.gz

archive/5_benchmark_alignments: $(addprefix archive/5_benchmark_alignments_,$(METHODS))
extract/5_benchmark_alignments: $(addprefix extract/5_benchmark_alignments_,$(METHODS))
clean/5_benchmark_alignments: $(addprefix clean/5_benchmark_alignments_,$(METHODS))

.PHONY: archive/5_benchmark_alignments
.PHONY: extract/5_benchmark_alignments
.PHONY: clean/5_benchmark_alignments
.PHONY: $(addprefix archive/5_benchmark_alignments_,$(METHODS))
.PHONY: $(addprefix extract/5_benchmark_alignments_,$(METHODS))
.PHONY: $(addprefix clean/5_benchmark_alignments_,$(METHODS))
