# fasta2tree
A simple computational workflow written in nextflow for inferring phylogenetic trees from a set of unaligned fasta files of protein sequences and visualizing these trees after rooting at midpoint. Input directory and options should be provided in the params.json file. The tool performs the following tasks given a directory of unaligned fasta files:
1. Multiple sequence alignment (MUSCLE, MAFFT or FSA)
2. Alignment trimming (trimAl)
3. Phylogenetic reconstruction (IQ-TREE v. 1.6.X or fasttree v. 2.1.X). If IQ-TREE is chosen, then model selection is performed automatically. If fasttree is selected the LG model is used (Le-Gascuel 2008 model). In IQ-TREE mode branch support is based on 2000 SH-aLRT replicates
5. Tree rerooting at midpoint and extraction of rooted tree figures in .svg format using ETE3

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
