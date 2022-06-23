# Metatranscriptome
# - preprocessing

FILES = glob_wildcards(config["dir_raw_metatranscriptome"]+"/{sample}_1.fastq")
SAMPLES = sorted(FILES.sample)

rule all:
    input:
        expand(config["dir_out"] + "/host_genome_removal/{sample}_metatranscriptome_unconc", sample=SAMPLES),
        expand(config["dir_out"] + "/host_genome_removal/{sample}_metatranscriptome.sam", sample=SAMPLES),    

rule trim:
    input:
        read1 = config["dir_raw_metatranscriptome"] + "/{sample}_1.fastq",
        read2 = config["dir_raw_metatranscriptome"] + "/{sample}_2.fastq"
    output:
        paired1 = config["dir_out"] + "/trim/{sample}_metatranscriptome_paired_1.fastq",
        unpaired1 = config["dir_out"] + "/trim/{sample}_metatranscriptome_unpaired_1.fastq",
        paired2 = config["dir_out"] + "/trim/{sample}_metatranscriptome_paired_2.fastq",
        unpaired2 = config["dir_out"] + "/trim/{sample}_metatranscriptome_unpaired_2.fastq"
    log:
        config["dir_out"] + "/trim/trimmomatic.log"
    params:
        adapter = config["adapter"]
    threads:
        config["threads"]
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 {input.read1} {input.read2} \
        {output.paired1} {output.unpaired1} {output.paired2} {output.unpaired2} \
        ILLUMINACLIP:{params.adapter}:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 2> {log}
        """

rule filt:
    input:
        read1 = rules.trim.output.paired1,
        read2 = rules.trim.output.paired2
    output:
        read1 = config["dir_out"] + "/filt/{sample}_metatranscriptome_1.fastq",
        read2 = config["dir_out"] + "/filt/{sample}_metatranscriptome_2.fastq"
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
        common1 = config["dir_out"] + "/filt/{sample}_metatranscriptome_1.fastq-common.out",
        common2 = config["dir_out"] + "/filt/{sample}_metatranscriptome_2.fastq-common.out",
        uniq1 = config["dir_out"] + "/filt/{sample}_metatranscriptome_1.fastq-unique.out",
        uniq2 = config["dir_out"] + "/filt/{sample}_metatranscriptome_2.fastq-unique.out"
    shell:
        "../tools/cmpfastq_pe {input.read1} {input.read2}"

rule merge:
    input:
        read1 = rules.cmp.output.common1,
        read2 = rules.cmp.output.common2
    output:
        config["dir_out"] + "/sortmerna/{sample}_metatranscriptome_merge.fastq"
    shell:
        "merge-paired-reads.sh {input.read1} {input.read2} {output}"

rule sortmerna_index:
    input:
        rfam_5_8s = config["dir_sortmerna"] + "/rfam-5.8s-database-id98.fasta",
        rfam_5s = config["dir_sortmerna"] + "/rfam-5s-database-id98.fasta",
        silva_arc_16s = config["dir_sortmerna"] + "/silva-arc-16s-id95.fasta",
        silva_arc_23s = config["dir_sortmerna"] + "/silva-arc-23s-id98.fasta",
        silva_bac_16s = config["dir_sortmerna"] + "/silva-bac-16s-id90.fasta",
        silva_bac_23s = config["dir_sortmerna"] + "/silva-bac-23s-id98.fasta",
        silva_euk_18s = config["dir_sortmerna"] + "/silva-euk-18s-id95.fasta",
        silva_euk_28s = config["dir_sortmerna"] + "/silva-euk-28s-id98.fasta"
    output:
        rfam_5_8s = config["sortmerna"] + "/index/rfam-5.8s-database-id98",
        rfam_5s = config["sortmerna"] + "/index/rfam-5s-database-id98",
        silva_arc_16s = config["sortmerna"] + "/index/silva-arc-16s-id95",
        silva_arc_23s = config["sortmerna"] + "/index/silva-arc-23s-id98",
        silva_bac_16s = config["sortmerna"] + "/index/silva-bac-16s-id90",
        silva_bac_23s = config["sortmerna"] + "/index/silva-bac-23s-id98",
        silva_euk_18s = config["sortmerna"] + "/index/silva-euk-18s-id95",
        silva_euk_28s = config["sortmerna"] + "/index/silva-euk-28s-id98"
    shell:
        """
        indexdb_rna --ref {input.rfam_5_8s},{output.rfam_5_8s}:{input.rfam_5s},{output.rfam_5s}:\
        {input.silva_arc_16s},{output.silva_arc_16s}:{input.silva_arc_23s},{output.silva_arc_23s}:\
        {input.silva_bac_16s},{output.silva_bac_16s}:{input.silva_bac_23s},{output.silva_bac_23s}:\
        {input.silva_euk_18s},{output.silva_euk_18s}:{input.silva_euk_28s},{output.silva_euk_28s}
        """
        
rule sortmerna:
    input:
        merged_pair = rule.merge.output,
        rfam_5_8s = rules.sortmerna_index.output.rfam_5_8s,
        rfam_5s = rules.sortmerna_index.output.rfam_5s,
        silva_arc_16s = rules.sortmerna_index.output.silva_arc_16s,
        silva_arc_23s = rules.sortmerna_index.output.silva_arc_23s,
        silva_bac_16s = rules.sortmerna_index.output.silva_bac_16s,
        silva_bac_23s = rules.sortmerna_index.output.silva_bac_23s,
        silva_euk_18s = rules.sortmerna_index.output.silva_euk_18s,
        silva_euk_28s = rules.sortmerna_index.output.silva_euk_28s,
    output:
        ribo = config["dir_out"] + "/sortmerna/{sample}_metatranscriptome_merge_rRNA.fastq",
        rem_ribo = config["dir_out"] + "/sortmerna/{sample}_metatranscriptome_merge_rRNA_removed.fastq"
    params:
        evalue = config["sortmerna_param"]["evalue"]
    log:
        config["dir_out"] + "/sortmerna/{sample}_metatranscriptome.log"
    threads:
        config["threads"]
    shell:
        """
        sortmerna --ref {input.rfam_5_8s},{input.rfam_5s},{input.silva_arc_16s},{input.silva_arc_23s},\
        {input.silva_bac_16s},{input.silva_bac_23s},{input.silva_euk_18s},{input.silva_euk_28s} \
        -a {threads} --reads {input.merged_pair} --aligned {output.ribo} --other {output.rem_ribo} --fastx \
        --paired_in -e {params.evalue} --log 2> {log}
        """

rule unmerge:
    input:
        rules.sortmerna.output.rem_ribo
    output:
        read1 = config["dir_out"] + "/sortmerna/{sample}_metatranscriptome_merge_rRNA_removed_1.fastq",
        read2 = config["dir_out"] + "/sortmerna/{sample}_metatranscriptome_merge_rRNA_removed_2.fastq"
    shell:
        "unmerge-paired-reads.sh {input} {output.read1} {output.read2}"

rule host_genome_index:
    input:
        genome = config["host_genome"]["fasta"]
    output:
        index = config["host_genome"]["index_hisat2"]
    threads:
        config["threads"]
    shell:
        "hisat2-build {input.genome} {output.index} -p {threads}"
        
rule host_genome_removal:
    input:
        read1 = rules.unmerge.output.read1,
        read2 = rules.unmerge.output.read2,
        index = rules.host_genome_index.output.index
    output:
        unconc = config["dir_out"] + "/host_genome_removal/{sample}_metatranscriptome_unconc",
        sam = config["dir_out"] + "/host_genome_removal/{sample}_metatranscriptome.sam"
    log:
        config["dir_out"] + "/host_genome_removal/{sample}_metatranscriptome.log"
    threads:
        config["threads"]
    shell:
        """
        hisat2 -x {input.index} \
        -1 {input.read1} -2 {input.read2} \
        -p {threads} -S {output.sam} --un-conc {output.unconc} 2> {log}
        """        
