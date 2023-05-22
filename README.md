# coati-testing

Pipeline to test COATi by simulating alignments with biologically-like indels.

Modify default values for number of sequences, max length of sequences,
and number of parallel jobs in `pipeline.sh`. Run for full pipeline with:
```bash
bash pipeline.sh
```

A detailed summary of this pipeline:

- Download pairwise homologous CDS genes between human and gorilla from ENSEMBL.
- Filter our sequences that are longer than threshold.
- Perform initial alignment.
- Split alignments that had gaps and those without.
- Extract gap patterns from initial alignments.
- Randomly insert gaps to sequences initially without, simulating data set of
"true alignments".
- Compare performance of aligners Clustal Omega, MACSE, MAFFT, PRANK, and COATi
in retrieving the data set of true alignments.

A detailed log of the results is saved in `results/results_summary.csv`.

This pipeline is used in "COATi: statistical pairwise alignment of protein coding sequences" and the supplementary data can be found [here](https://figshare.com/articles/dataset/Supplementary_data_for_COATi_statistical_pairwise_alignment_of_protein_coding_sequences_/23064227).
