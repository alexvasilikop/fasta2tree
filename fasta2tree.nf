#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Using a params file specifying 1) path to fasta files 2) multiple seq alignment software and 3) tree reconstruction software, perform:
// 1) multiple sequence alignments (options: mafft, muscle, fsa)
// 2) trimming of individual alignment files (trimal)
// 3) phylogenetic reconstruction (options: iqtree, fasttree)
// 4) tree midpoint rerooting, visualization and export in .svg format using ETE3 (Huerta-Cepas et al. 2016, MBE)

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

    publishDir params.alignments_outdir, mode: 'copy', pattern: "*.trimmed.fa"
    debug true

    input:
    path alignment

    output:
    path "${alignment.baseName}.trimmed.fa", emit: fasta

    script:
    """
    trimal -in ${alignment} -out "${alignment.baseName}.trimmed.fa" -gappyout -fasta
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
    path "${aligned_fasta.baseName}.log", emit: log

    script:
    if (params.tree_software=="fasttree") {
        """
        fasttree -lg -log "${aligned_fasta.baseName}.log" ${aligned_fasta} > "${aligned_fasta.baseName}.treefile"
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

process TREE_REROOTING_VISUALIZATION {

    publishDir params.tree_figs_outdir, mode: 'copy', pattern: "*rerooted.svg"
    debug true
    
    input:
    path treefiles

    output:
    path "${treefiles.baseName}.rerooted.svg", emit: svg

    script:
    """
    #!/usr/bin/env python3
    from ete3 import Tree, TreeStyle

    t=Tree("${treefiles}", format=0)

    midpoint_ancestor=t.get_midpoint_outgroup()
    t.set_outgroup(midpoint_ancestor)

    ts = TreeStyle()
    ts.show_branch_support = True
    t.render("${treefiles.baseName}.rerooted.svg", tree_style=ts)
    """
}

workflow {
    fasta_ch = Channel.fromPath(params.fasta, checkIfExists: true)

    //Align
    ALIGNMENT(fasta_ch)
    //Trim alignments
    TRIM_ALIGNMENTS(ALIGNMENT.out.fasta)
    //Phylogenetic tree reconstruction
    PHYLOGENY(TRIM_ALIGNMENTS.fasta.out)
    PHYLOGENY.out.log.view()
    PHYLOGENY.out.tree.view()

    //visualize the trees
    TREE_REROOTING_VISUALIZATION(PHYLOGENY.out.tree)
    TREE_REROOTING_VISUALIZATION.out.svg.view()
}
