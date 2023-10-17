# fasta2tree
A simple computational workflow for inferring phylogenetic trees from a set of unaligned fasta files (input directory is provided in the params.json file) and 
visualizing these trees after rooting at midpoint.

## Usage
``` nextflow run fasta2tree.nf -params-file params.json```

## Requirements
The following should be on path:
- Nextflow

Multiple sequence alignment
- MAFFT
- MUSCLE
- FSA

Alignment trimming:
- trimAl

Phylogenetic reconstruction software:
- IQ-TREE
- FASTTREE

In addition installation of ETE3 (http://etetoolkit.org/) and python3 is required.
