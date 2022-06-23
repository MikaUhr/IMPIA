#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Merging contigs

rule all:
    input:
        config["merge_dir"] + "/merged_out.fasta"

rule merge1:
    input:
        contig1 = config["contig1"],
        contig2 = config["contig2"]
    output:
        config["dir_out"]
    params:
        hco = config["merge"]["hco"],
        c = config["merge"]["c"],
        l = config["merge"]["l"],
        mlength = config["merge"]["mlength"] 
    shell:
        """
        cd {output}
        merge_wrapper.py -hco ${params.hco} -c ${params.c} -l ${params.l} -ml ${params.mlength} {input.contig1} {input.contig2}
        """

