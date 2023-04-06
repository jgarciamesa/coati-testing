################################################################################
# aligner recipes: coati, prank, mafft, clustal omega, macse                   #
################################################################################

aln/tri-mg/%: $(RAW_PATH)/%
	@echo -ne "coati tri-mg align $*\t\t\r"
	@./bin/coati-alignpair $< -m tri-mg -o $@

aln/mar-mg/%: $(RAW_PATH)/%
	@echo -ne "mcoati ar-mg align $*\t\t\r"
	@./bin/coati-alignpair $< -m mar-mg -o $@ 2> /dev/null

aln/dna/%: $(RAW_PATH)/%
	@echo -ne "coati dna align $*\t\t\r"
	@./bin/coati-alignpair $< -m dna -o $@

aln/tri-ecm/%: $(RAW_PATH)/%
	@echo -ne "coati tri-ecm align $*\t\t\r"
	@./bin/coati-alignpair $< -m tri-ecm -o $@

aln/mar-ecm/%: $(RAW_PATH)/%
	@echo -ne "coati mar-ecm align $*\t\t\r"
	@./bin/coati-alignpair $< -m mar-ecm -o $@

aln/prank/%: $(RAW_PATH)/%
	@echo -ne "prank align $*\t\t\r"
	@./bin/prank -codon -d="$<" -o=$@ -quiet &>/dev/null
	@if [ -f $@.best.fas ]; then mv $@.best.fas $@; else touch $@; fi

aln/mafft/%: $(RAW_PATH)/%
	@echo -ne "mafft align $*\t\t\r"
	@./bin/mafft --quiet --preservecase --globalpair --maxiterate 1000 $< > $@

aln/clustalo/%: $(RAW_PATH)/%
	@echo -ne "clustalo align $*\t\t\r"
	@./bin/clustalo --force -i $< -o $@

aln/macse/%: $(RAW_PATH)/%
	@echo -ne "macse align $*\t\t\r"
	@java -jar ./bin/macse_v2.06.jar -fs_lr 10 -stop_lr 10 -prog alignSequences -seq $< -out_NT $@ -out_AA output_AA_$* &>/dev/null
	@rm output_AA_$*
	@sed -i 's/!/-/g' $@

clean_initial_%:
	rm aln/$*/*

clean_initial_all:
	rm -f aln/tri-mg/*.fasta
	rm -f aln/mar-mg/*.fasta
	rm -f aln/dna/*.fasta
	rm -f aln/ecm/*.fasta
	rm -f aln/macse/*.fasta
	rm -f aln/mar-ecm/*.fasta
	rm -f aln/prank/*.fasta
	rm -f aln/mafft/*.fasta
	rm -f aln/clustalo/*.fasta

################################################################################
# coati models vs prank vs mafft vs clustalO using REFERENCE aligments         #
################################################################################

aln/ref/tri-mg/%: $(REF_PATH)/%
	@echo -ne "coati tri-mg align $*\r"
	@./bin/coati-alignpair $< -m tri-mg -o $@

aln/ref/mar-mg/%: $(REF_PATH)/%
	@echo -ne "coati mar-mg align $*\r"
	@./bin/coati-alignpair $< -m mar-mg -o $@

aln/ref/dna/%: $(REF_PATH)/%
	@echo -ne "coati dna align $*\r"
	@./bin/coati-alignpair $< -m dna -o $@

aln/ref/tri-ecm/%: $(REF_PATH)/%
	@echo -ne "coati tri-ecm align $*\r"
	@./bin/coati-alignpair $< -m tri-ecm -o $@

aln/ref/mar-ecm/%: $(REF_PATH)/%
	@echo -ne "coati mar-ecm align $*\r"
	@./bin/coati-alignpair $< -m mar-ecm -o $@

aln/ref/prank/%: $(REF_PATH)/%
	@echo -ne "prank align $*\r"
	@./bin/prank -codon -d="$<" -o=$@ -quiet &>/dev/null
	@if [ -f $@.best.fas ]; then mv $@.best.fas $@; else touch $@; fi

aln/ref/mafft/%: $(REF_PATH)/%
	@echo -ne "mafft align $*\r"
	@./bin/mafft --quiet --preservecase --globalpair --maxiterate 1000 $< > $@

aln/ref/clustalo/%: $(REF_PATH)/%
	@echo -ne "clustalo align $*\r"
	@./bin/clustalo -i $< -o $@

aln/ref/macse/%: $(REF_PATH)/%
	@echo -ne "macse align $*\r"
	@java -jar bin/macse_v2.06.jar -fs_lr 10 -stop_lr 10 -prog alignSequences -seq $< -out_NT $@ -out_AA output_AA_$* &>/dev/null
	@rm output_AA_$*
	@sed -i 's/!/-/g' $@
	

clean_reference_%:
	rm aln/ref/$*/*

clean_reference_all:
	rm -f aln/ref/tri-mg/*
	rm -f aln/ref/mar-mg/*
	rm -f aln/ref/dna/*
	rm -f aln/ref/ecm/*
	rm -f aln/ref/mecm/*
	rm -f aln/ref/prank/*
	rm -f aln/ref/mafft/*
	rm -f aln/ref/clustalo/*

