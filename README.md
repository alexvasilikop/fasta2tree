# fasta2tree
A simple computational workflow for inferring phylogenetic trees from a set of unaligned fasta files

## Usage
``` nextflow run fasta2tree.nf -params-file params.json```

## Requirements
- Nextflow

Alignment trimming:
- trimal

Phylogenetic reconstruction software:
- IQ-TREE
- FASTTREE

Multiple sequence alignment
- MAFFT
- MUSCLE
- FSA
