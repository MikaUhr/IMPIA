#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Metagenome
# - preprocessing
# - assembly
# - scaffolding

FILES = glob_wildcards(config["dir_raw_metagenome"]+"/{sample}_1.fastq")
SAMPLES = sorted(FILES.sample)

rule all:
    input:
        expand(config["dir_out"] + "/scaffold/{sample}", sample=SAMPLES)

rule trim:
    input:
        read1 = config["dir_raw_metagenome"] + "/{sample}_1.fastq",
        read2 = config["dir_raw_metagenome"] + "/{sample}_2.fastq",
        adapter = config["adapter"]
    output:
        paired1 = config["dir_out"] + "/trim/{sample}_metagenome_paired_1.fastq",
        unpaired1 = config["dir_out"] + "/trim/{sample}_metagenome_unpaired_1.fastq",
        paired2 = config["dir_out"] + "/trim/{sample}_metagenome_paired_2.fastq",
        unpaired2 = config["dir_out"] + "/trim/{sample}_metagenome_unpaired_2.fastq"
    log:
        config["dir_out"] + "/trim/trimmomatic.log"        
    threads:
        config["threads"]
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 {input.read1} {input.read2} \
        {output.paired1} {output.unpaired1} {output.paired2} {output.unpaired2} \
        ILLUMINACLIP:{input.adapter}:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 2> {log}
        """

rule filt:
    input:
        read1 = rules.trim.output.paired1,
        read2 = rules.trim.output.paired2
    output:
        read1 = config["dir_out"] + "/filt/{sample}_metagenome_1.fastq",
        read2 = config["dir_out"] + "/filt/{sample}_metagenome_2.fastq"
    shell:
        """
        fastq_quality_filter -q 20 -p 80 -i {input.read1} -o {output.read1} -Q33
        fastq_quality_filter -q 20 -p 80 -i {input.read2} -o {output.read2} -Q33
        """

rule cmp:
    input:
        read1 = rules.filt.output.read1,
        read2 = rules.filt.output.read2
    output:
        common1 = config["dir_out"] + "/filt/{sample}_metagenome_1.fastq-common.out",
        common2 = config["dir_out"] + "/filt/{sample}_metagenome_2.fastq-common.out",
        uniq1 = config["dir_out"] + "/filt/{sample}_metagenome_1.fastq-unique.out",
        uniq2 = config["dir_out"] + "/filt/{sample}_metagenome_2.fastq-unique.out"
    shell:
        "../tools/cmpfastq_pe {input.read1} {input.read2}"

rule host_genome_index:
    input:
        genome = config["host_genome"]["fasta"]
    output:
        index = config["host_genome"]["index_bowtie2"]
    threads:
        config["threads"]
    shell:
        "bowtie2-build {input.genome} {output.index} --threads {threads}"
        
rule host_genome_removal:
    input:
        read1 = rules.cmp.output.common1,
        read2 = rules.cmp.output.common2,
        index = rules.host_genome_index.output.index
    output:
        unconc = config["dir_out"] + "/host_genome_removal/{sample}_metagenome_unconc",
        sam = config["dir_out"] + "/host_genome_removal/{sample}_metagenome.sam"
    log:
        config["dir_out"] + "/host_genome_removal/{sample}_metagenome.log"
    threads:
        config["threads"]
    shell:
        """
        bowtie2 -x {input.index} -1 {input.read1} -2 {input.read2} \
        -S {output.sam} --un-conc {output.unconc} --threads {threads} 2> {log}
        """
        
rule assembly:
    input:
        read1 = config["dir_out"] + "/host_genome_removal/{sample}_metagenome_unconc.1",
        read2 = config["dir_out"] + "/host_genome_removal/{sample}_metagenome_unconc.2"
    output:
        config["dir_out"] + "/assembly/{sample}"
    log:
        config["dir_out"] + "/assembly/{sample}.log"
    params:
        k_min = config["megahit_param"]["k_min"],
        k_max = config["megahit_param"]["k_max"],
        k_step = config["megahit_param"]["k_step"],
        prune_depth = config["megahit_param"]["prune_depth"]
    threads:
        config["threads"]
    shell:
        """
        megahit -1 {input.read1} -2 {input.read2} -o {output} \
        --k-min {params.k_min} --k-max {params.k_max} --k-step {params.k_step} \
        --prune-depth {params.prune_depth} -t {threads}
        """

rule contig_index:
    input:
        contig = rules.assembly.output + "/final.contigs.fa" 
    output:
        index = config["dir_out"] + "/assembly/bowtie2_index/{sample}"
    threads:
        config["threads"]
    shell:
        "bowtie2-build {input.contig} {output.index} --threads {threads}"
        
rule map_contig:
    input:
        read1 = rules.assembly.input.read1,
        read2 = rules.assembly.input.read2,
        index = rules.contig_index.output.index
    output:
        bam_read1 = config["dir_out"] + "/assembly/{sample}_1.bam",
        bam_read2 = config["dir_out"] + "/assembly/{sample}_2.bam"
    threads:
        config["threads"]
    params:
        memory = str(round(float(config["memory_per_core"]) * 0.8 - 1)) + "G"
    shell:
        """
        bowtie2 -p {threads} -x {input.index} -U {input.read1} \
        | samtools view -bS - | samtools sort -m {params.memory} -@ {threads} -o {output.bam_read1} -
        bowtie2 -p {threads} -x {input.index} -U {input.read2} \
        | samtools view -bS - | samtools sort -m {params.memory} -@ {threads} -o {output.bam_read2} -
        """

rule bam_merge:
    input:
        bam_read1 = rules.map_contig.output.bam_read1,
        bam_read2 = rules.map_contig.output.bam_read2
    output:
        config["dir_out"] + "/assembly/{sample}_merge.bam"
    shell:
        "samtools merge -f {output} {input.bam_read1} {input.bam_read2}"

rule bam_sort:
    input:
        rules.bam_merge.output
    output:
        config["dir_out"] + "/assembly/{sample}_merge_sort.bam"
    params:
        memory = str(round(float(config["memory_per_core"]) * 0.8 - 1)) + "G"
    threads:
        config["threads"]
    shell:
        "samtools sort {input} -n -m {params.memory} -@ {threads} -o {output}"
        
rule bam_to_sam:
    input:
        rules.bam_sort.output
    output:
        config["dir_out"] + "/assembly/{sample}_merge_sort.sam"
    shell:
        "samtools view {input} > {output}"
        
rule scaffold:
    input:
        contig = rules.assembly.output + "/final.contigs.fa",
        sam = rules.bam_to_sam.output
    output:
        config["dir_out"] + "/scaffold/{sample}"
    shell:
        "../tools/OPERA-LG_v2.0.6/bin/OPERA-LG {input.contig} {input.sam} {output}"
        

