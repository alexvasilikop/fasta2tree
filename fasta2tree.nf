#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Using a params file specifying 1) path to fasta files 2) multiple seq alignment software and 3) tree reconstruction software, perform:
// 1) multiple sequence alignments (options: mafft, muscle, fsa)
// 2) trimming of individual alignment files (trimal)
// 3) phylogenetic reconstruction (options: iqtree, fasttree)

process ALIGNMENT {

    publishDir params.alignments_outdir, mode: 'copy', pattern: "*aligned.fasta"
    debug true

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}.aligned.fasta", emit: fasta

    script:

    if (params.al_software=="mafft") {
        
        //Use auto otions
        println("Selected mafft for MSA...")
        """
        mafft --auto --thread ${params.cpus} ${fasta_file} > "${fasta_file.baseName}.aligned.fasta"
        """
    }

    else if (params.al_software=="muscle") {

        //Use default options
        println("Selected muscle for MSA...")
        """
        muscle -in ${fasta_file} -out "${fasta_file.baseName}.aligned.fasta"
        """
    }

    else if (params.al_software=="fsa") {

        //Use default options
        println("Selected fsa for MSA...")
        """
        fsa --fast ${fasta_file} > "${fasta_file.baseName}.aligned.fasta"
        """
    }

    else {
        println("### Incorrect alignment software specified in parameters file: \"${params.al_software}\" ###\n")
        """
        ${params.al_software} -h
        """
    }

}

process TRIM_ALIGNMENTS {

    publishDir params.alignments_outdir, mode: 'copy'
    debug true

    input:
    path alignment

    output:
    path "${alignment.baseName}.aln.trimmed.fa"

    script:
    """
    trimal -in ${alignment} -out "${alignment.baseName}.aln.trimmed.fa" -gappyout -fasta
    """
}

process PHYLOGENY {

    publishDir params.tree_outdir, mode: 'copy', pattern: "*treefile"
    publishDir params.tree_log_outdir, mode: 'copy', pattern: "*log"
    debug true

    input:
    path aligned_fasta

    output:
    path "${aligned_fasta.baseName}.treefile", emit: tree
    path "${aligned_fasta.baseName}*.log", emit: log

    script:
    if (params.tree_software=="fasttree") {
        """
        fasttree -log "${aligned_fasta.baseName}.fasttree.log" ${aligned_fasta} > "${aligned_fasta.baseName}.treefile"
        """
    }

    else if (params.tree_software=="iqtree") {
        """
        iqtree -s ${aligned_fasta} -m MFP -nt ${params.cpus} -pre ${aligned_fasta.baseName} -safe -alrt 2000
        """
    }

    else {
        """
        echo -e "### Incorrect alignment software specified in parameters file: \"${params.al_software}\" ###\n"
        ${params.al_software} -h
        """
    }
}

workflow {
    fasta_ch = Channel.fromPath(params.fasta, checkIfExists: true)

    //Align
    ALIGNMENT(fasta_ch)
    //Trim alignments
    TRIM_ALIGNMENTS(ALIGNMENT.out.fasta)
    //Phylogenetic tree reconstruction
    PHYLOGENY(TRIM_ALIGNMENTS.out)
    PHYLOGENY.out.log.view()
    PHYLOGENY.out.tree.view()
}