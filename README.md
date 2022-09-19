# coati-testing

Pipeline to test COATi by simulating alignments with biologically-like indels.

Modify default values for species, number of sequences, max length of sequences,
and number of parallel jobs in `pipeline.sh`. Run for full pipeline with:
```bash
bash pipeline.sh
```

A detailed summary of this pipeline:

- Download pairwise homologous CDS genes between human and either gorilla,
mouse, or *Drosophila melanogaster* from ENSEMBL.
- Filter our sequences that are longer than threshold.
- Perform initial alignment.
- Split alignments that had gaps and those without.
- Extract gap patterns from initial alignments.
- Randomly insert gaps to sequences initially without, simulating data set of
"true alignments".
- Compare performance of aligners Clustal Omega, MACSE, MAFFT, PRANK, and
different COATi models in retrieving the data set of true alignments.

Results comparing the distance from retrieved and simulated alignments are saved
in `data/${SPECIES}/dseq_summary.csv`. Lower values represent a better
performance. Other informative metrics are displayed on the terminal upon end of
the pipeline.