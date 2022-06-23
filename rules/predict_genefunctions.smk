#!/usr/bin/env python
# coding: utf-8

# In[ ]:


FILES = glob_wildcards(config["dir_input"]+"/annotation_cogid_{individual}.tsv")
INDIVIDUALS = sorted(FILES.individual)

rule all:
    input:
        expand(config["dir_input"]+"/tpm_join_{individual}.txt", individual=INDIVIDUALS),
        config["dir_output"] + "/gene_expression/spatial_covariation.csv",

rule concat_fasta:
    input:
        expand(config["dir_input"] +  "/unknown_aminoacid_{individual}.fasta", individual=INDIVIDUALS),
    output:
        config["dir_output"] + "/unknown_aminoacid.fasta"
    shell:
        "cat {input} > {output}"

rule make_gene_seq_database:
    input:
        rules.concat_fasta.output
    output:
        config["dir_output"] + "/unknowngene_cluster/mmseqs_database/mmseqs.DB"
    shell:
        "mmseqs createdb {input} {output}"
        
rule gene_clustering:
    input:
        rules.make_gene_seq_database.output
    output:
        cluster = config["dir_output"] + "/unknowngene_cluster/cluster.list",
        dir_tmp = config["dir_output"] + "/unknowngene_cluster/tmp"
    params:
        cluster_mode = config["gene_clustering"]["cluster_mode"],
        cov_mode = config["gene_clustering"]["cov_mode"],
        c = config["gene_clustering"]["c"],
        s = config["gene_clustering"]["s"],
        kmer = config["gene_clustering"]["kmer"]
    threads:
        config["threads"]
    shell:
        """
        mmseqs cluster --threads {threads} --cluster-mode {params.cluster_mode} --cov-mode {params.cov_mode} \
        -c {params.c} -s {params.s} --kmer-per-seq {params.kmer} {input} {output.cluster} {output.dir_tmp}
        """
        
rule gene_clustering_tsv:
    input:
        mmseqs_db = rules.make_gene_seq_database.output,
        cluster = rules.gene_clustering.output.cluster
    output:
        config["dir_output"] + "/unknowngene_cluster/cluster.tsv"
    shell:
        "mmseqs createtsv {input.mmseq_db} {input.mmseq_db} {input.cluster} {output}"
        
rule convert_cluster:
    input: 
        rules.gene_clustering_tsv.output
    output:
        config["dir_output"] + "/unknowngene_cluster/cluster_ugc.tsv"
    shell:
        """
        cat {input} | awk 'BEGIN {l="X"} {FS=OFS="\t"} {if($1!=l){i++} {print "UGC"i,$2} {l=$1}}' > {output}
        """

rule merge_genecluster:
    input:
        unknowngene = rules.convert_cluster.output,
        tpm = config["dir_input"]+"/tpm_join_{individual}.txt",
        annotatedgene = config["dir_input"]+"/annotation_cogid_{individual}.tsv"
        individual_name = "{individual}"
    output:
        config["dir_output"] + "/gene_expression/tpm_{individual}.tsv"
    shell:
        """
        python ../scripts/merge_geneclusters.py {input.unknowngene} {input.tpm} \
        {input.annotatedgene} {input.individual_name} {output}
        """
        
rule merge_individuals:
    input:
        expand(config["dir_input"]+"/tpm_join_{individual}.txt", individual=INDIVIDUALS),
    output:
        final = config["dir_output"] + "/gene_expression/tpm_merged.tsv",
        tmp = config["dir_output"] + "/gene_expression/tpm_merged_tmp.tsv",
    shell:
        "../scripts/multijoin_tpm.sh {output.final} {output.tmp} {input}"

rule removal:
    input:
        rules.merge_individuals.output.tmp
    shell:
        "rm {input}"
        
rule spatial_covariation:
    input:
        tpm = rules.merge_individuals.output.final,
        pseudo_expression = config["spatial_covariation"]["pseudo_expression"],
        variance = config["spatial_covariation"]["variance"],
        b = config["spatial_covariation"]["b"],
        distance = config["spatial_covariation"]["distance"],
        L_threshold = config["spatial_covariation"]["L_threshold"]
    output:
        linked_genepair = config["dir_output"] + "/gene_expression/spatial_covariation.csv"
    shell:
        """
        python ../scripts/spatial_covariation.py {input.tpm} {output} \
        {input.pseudo_expression} {input.variance} {input.b} {input.distance} {input.L_threshold}
        """

