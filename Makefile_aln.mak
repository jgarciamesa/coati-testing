################################################################################
# coati models vs prank vs mafft alignments                                    #
################################################################################

aln/coati/%: $(RAW_PATH)/%
	@echo -ne "coati align $*\t\t\r"
	@./bin/coati-alignpair $< -m coati -o $@

aln/mcoati/%: $(RAW_PATH)/%
	@echo -ne "mcoati align $*\t\t\r"
	@./bin/coati-alignpair $< -m m-coati -o $@

aln/dna/%: $(RAW_PATH)/%
	@echo -ne "dna align $*\t\t\r"
	@./bin/coati-alignpair $< -m dna -o $@

aln/ecm/%: $(RAW_PATH)/%
	@echo -ne "ecm align $*\t\t\r"
	@./bin/coati-alignpair $< -m ecm -o $@

aln/mecm/%: $(RAW_PATH)/%
	@echo -ne "mecm align $*\t\t\r"
	@./bin/coati-alignpair $< -m m-ecm -o $@

aln/prank/%: $(RAW_PATH)/%
	@echo -ne "prank align $*\t\t\r"
	@./bin/prank -codon -d="$<" -o=$@ -quiet &>/dev/null
	@if [ -f $@.best.fas ]; then mv $@.best.fas $@; else touch $@; fi

aln/mafft/%: $(RAW_PATH)/%
	@echo -ne "mafft align $*\t\t\r"
	@mafft --quiet --preservecase --globalpair --maxiterate 1000 $< > $@

aln/clustalo/%: $(RAW_PATH)/%
	@echo -ne "clustalo align $*\t\t\r"
	@clustalo --force -i $< -o $@

aln/macse/%: $(RAW_PATH)/%
	@echo -ne "macse align $*\t\t\r"
	@java -jar ./bin/macse_v2.06.jar -fs_lr 10 -stop_lr 10 -prog alignSequences -seq $< -out_NT $@ -out_AA output_AA_$* &>/dev/null
	@rm output_AA_$*
	@sed -i 's/!/-/g' $@

clean_initial_%:
	rm aln/$*/*

clean_initial_all:
	rm -f aln/coati/*.fasta
	rm -f aln/mcoati/*.fasta
	rm -f aln/dna/*.fasta
	rm -f aln/ecm/*.fasta
	rm -f aln/macse/*.fasta
	rm -f aln/mecm/*.fasta
	rm -f aln/prank/*.fasta
	rm -f aln/mafft/*.fasta
	rm -f aln/clustalo/*.fasta

################################################################################
# coati models vs prank vs mafft vs clustalO using REFERENCE aligments         #
################################################################################

aln/ref/coati/%: $(REF_PATH)/%
	@echo -ne "coati align $*\r"
	@./bin/coati-alignpair $< -m coati -o $@

aln/ref/mcoati/%: $(REF_PATH)/%
	@echo -ne "mcoati align $*\r"
	@./bin/coati-alignpair $< -m m-coati -o $@

aln/ref/dpmcoati/%: $(REF_PATH)/%
	@echo -ne "dpmcoati align $*\r"
	@./bin/coati-alignpair $< -m dp-mcoati -o $@

aln/ref/dna/%: $(REF_PATH)/%
	@echo -ne "dna align $*\r"
	@./bin/coati-alignpair $< -m dna -o $@

aln/ref/ecm/%: $(REF_PATH)/%
	@echo -ne "ecm align $*\r"
	@./bin/coati-alignpair $< -m ecm -o $@

aln/ref/mecm/%: $(REF_PATH)/%
	@echo -ne "mecm align $*\r"
	@./bin/coati-alignpair $< -m m-ecm -o $@

aln/ref/prank/%: $(REF_PATH)/%
	@echo -ne "prank align $*\r"
	@./bin/prank -codon -d="$<" -o=$@ -quiet &>/dev/null
	@if [ -f $@.best.fas ]; then mv $@.best.fas $@; else touch $@; fi

aln/ref/mafft/%: $(REF_PATH)/%
	@echo -ne "mafft align $*\r"
	@mafft --quiet --preservecase --globalpair --maxiterate 1000 $< > $@

aln/ref/clustalo/%: $(REF_PATH)/%
	@echo -ne "clustalo align $*\r"
	@clustalo -i $< -o $@

aln/ref/macse/%: $(REF_PATH)/%
	@echo -ne "macse align $*\r"
	@java -jar bin/macse_v2.06.jar -fs_lr 10 -stop_lr 10 -prog alignSequences -seq $< -out_NT $@ -out_AA output_AA_$* &>/dev/null
	@rm output_AA_$*
	@sed -i 's/!/-/g' $@
	

#aln/ref/baliphy/%.fasta:
#	./bin/bali-phy --smodel mg94 --iterations 100 data/$(SPECIES)/trim_ref_ungap/$*.fasta
#	tac $*-1/C1.P1.fastas | sed '/iterations/Q' | tac > $@
#	rm -r $*-1

clean_reference_%:
	rm aln/ref/$*/*

clean_reference_all:
	rm -f aln/ref/coati/*
	rm -f aln/ref/mcoati/*
	rm -f aln/ref/dpmcoati/*
	rm -f aln/ref/dna/*
	rm -f aln/ref/ecm/*
	rm -f aln/ref/mecm/*
	rm -f aln/ref/prank/*
	rm -f aln/ref/mafft/*
	rm -f aln/ref/clustalo/*
	rm -f aln/ref/baliphy/*

