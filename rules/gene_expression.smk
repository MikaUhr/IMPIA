FILES = glob_wildcards(config["dir_metatranscriptome"]+"/{sample}_1.fastq")
SAMPLES = sorted(FILES.sample)

rule all:
    input:
        expand(config["dir_out"] + "/rna_mapping/{sample}_tpm.tsv", sample=SAMPLES),
        config["dir_out"] + "/rna_mapping/tpm_join_" + config["individual_name"] + ".txt",

rule convert_fasta:
    input: 
        config["reconstructed_metagenome"]["fasta"]
    output:
        config["dir_out"] + "/reconstructed_metagenome/merged_metagenome.fasta"
    shell:
        """
        cat {input} | awk '/^>/ {print ">Merged_" ++i; next}{print}' > {output}
        """
rule gene_prediction:
    input:
        reference = rules.convert_fasta.output
    output:
        gene_nucleotide = config["dir_out"] + "/gene_prediction/gene_nucleotide.fasta",
        gene_aminoacid = config["dir_out"] + "/gene_prediction/gene_aminoacid.fasta",
        gff = config["dir_out"] + "/gene_prediction/gene.gff"
    params: 
        metagememark_mod = config["metagenemark_mod"]
    log:
        config["dir_out"] + "/gene_prediction/metagenemark.log"
    shell:
        """
        gmhmmp -m {params.metagenemark_mod} {output.reference} -f G -D {output.gene_nucleotide} \
        -A {output.gene_aminoacid} -L {log} -o {output.gff}
        """

rule gene_prediction_convert:
    input:
        gene_nucleotide = rules.gene_prediction.output.gene_nucleotide,
        gene_aminoacid = rules.gene_prediction.output.gene_aminoacid
    output:
        gene_nucleotide = config["dir_out"] + "/gene_prediction/gene_nucleotide_convert.fasta",
        gene_aminoacid = config["dir_out"] + "/gene_prediction/gene_aminoacid_convert.fasta",
    shell:
        """
        cat {input.gene_nucleotide} | awk -F "|" '{print $1}' > {output.gene_nucleotide}
        cat {input.gene_aminoacid} | awk -F "|" '{print $1}' > {output.gene_nucleotide}
        """
        
rule gene_length:
    input:
        rules.gene_prediction.output.gff
    output:
        length = config["dir_out"] + "/gene_prediction/gene_length.tsv"
    shell:
        """
        cat {input.gff} | awk -F "\t" '$3=="CDS" {{FS=OFS="\t"} {print $9,$5-$4+1}} > {output.length}'
        """
        
rule make_cog_database:
    input:
        cog = config["COG"]["COG_fasta"]
    output:
        config["COG"]["COG_diamond"]
    threads:
        config["threads"]
    shell:
        "diamond makedb --in {input.cog} --db {output} --threads {threads}"
        
rule gene_annotation:
    input:
        cog_diamond = rules.make_cog_database.output,
        gene_aminoacid = rules.gene_prediction_convert.output.gene_aminoacid
    output:
        config["dir_out"] + "/gene_annotation/annotation.tsv"
    params:
        tmp = config["gene_annotation"]["directory_tmp"],
        evalue = config["gene_annotation"] ["evalue"],
        query_cover = config["gene_annotation"]["query_cover"]
        min_orf = config["gene_annotation"]["min_orf"]
    shell:
        """
        diamond blastp -t {params.tmp} --db {input.cog_diamond} \
        --query {input.gene_aminoacid} --out {output} \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps qlen slen \
        --evalue {params.evalue} --query-cover {params.query_cover} --max-target-seqs 1 --min-orf {params.min_orf} --threads {threads} 
        """ 
        
rule cogid:
    input:
        config["COG"]["COG_csv"]
    output:
        config["dir_out"] + "/gene_annotation/protid_cogid.tsv"
    shell:
        "cat {input} | awk -F "," '{print $3,$7}' | sort -k 1,1 > {output}"

rule gene_annotation_cogid:
    input:
        blastp = rules.gene_annotation.output,
        cog = rules.cogid.output,
        individual_name = config["individual_name"]
    output:
        config["dir_out"] + "/gene_annotation/annotation_cogid_" + config["individual_name"] + ".tsv"
    shell:
        """
        sample_name={input.individual_name}; \
        cat {input.blastp} | awk -F "\t" '{print $2,$1}' | uniq -f 1 | sed 's/ //g' \
        | awk -F "|" '{print $2,$5}' | sort -k 1,1 | join -a 1 -o 1.2 2.2 - {input.cog} | sort -k 1,1 | tr ' ' '\t' \
        | sed "s/^/${sample_name}_/g" > {output}
        """

rule unknowngene_list:
    input:
        gff = rules.gene_prediction.output.gff,
        gene_known = rules.gene_annotation_cogid.output
    output:
        config["dir_out"] + "/gene_annotation/unknown.list"
    shell:
        """
        cat {input.gff} | awk -F "\t" '$9~"gene_id=" {print $9}' | sed 's/gene_id=/gene_/g' \
        | sort -k 1,1 | join -v 1 - {input.gene_known} > {output}
        """

rule unknowngene_fastq:
    input:
        unknown_list = rules.unknowngene_list.output,
        gene_aminoacid = rules.rules.gene_prediction_convert.output.gene_aminoacid,
        individual_name = config["individual_name"]
    output:
        config["dir_out"] + "/gene_annotation/unknown_aminoacid_" + config["individual_name"] + ".fasta"
    shell:
        """
        sample_name={input.individual_name}; \
        seqtk subseq {input.gene_aminoacid} {input.unknown_list} | sed "s/>/>${sample_name}_/g" > {output}
        """

rule metagenome_index:
    input:
        genome = config["reconstructed_metagenome"]["fasta"]
    output:
        index = config["reconstructed_metagenome"]["index_bowtie2"]
    threads:
        config["threads"]
    shell:
        "bowtie2-build {input} {index} -p {threads}"

rule rna_mapping:
    input:
        read1 = config["dir_metatranscriptome"]+"/{sample}_1.fastq",
        read2 = config["dir_metatranscriptome"]+"/{sample}_2.fastq" ,
        index = rules.metagenome_index.output.index
    output:
        bam = config["dir_out"] + "/rna_mapping/{sample}.bam"
    params:
        memory = str(round(float(config["memory_per_core"]) * 0.8 - 1)) + "G"
    threads:
        config["threads"]
    log:
        config["dir_out"] + "/rna_mapping/{sample}.log"
    shell:
        """
        bowtie2 -p {threads} -x {input.index} -1 ${input.read1} -2 {input.read2} 2> {log}\
        | samtools view -bS - | samtools sort -n -m {params.memory} -@ {threads} -o {output.bam} -
        """
        
rule bam_to_sam:
    input:
        bam = rules.rna_mapping.output.bam
    output:
        sam = config["dir_out"] + "/rna_mapping/{sample}_merge_sort.sam"
    shell:
        "samtools view {input.bam} > {output.sam}"
        
rule count_rna:
    input:
        sam = rules.bam_to_sam.output.sam,
        gff = rules.gene_prediction.output.gff
    output:
        config["dir_out"] + "/rna_mapping/{sample}_readcount.tsv"
    shell:
        "htseq-count -s reverse --type CDS -f sam -i gene_id {input.sam} {input.gff} > {output}"

rule conv_tpm:
    input:
        length = rules.gene_length.output.length,
        read_count = rules.count_rna.output,
        sample_name = "{sample}",
        individual_name = config["individual_name"]
    output:
        tpm = config["dir_out"] + "/rna_mapping/{sample}_tpm.tsv"
    shell:
        "python ../scripts/TPM.py {input.length} {input.read_count} {output.tpm} {input.sample_name} {input.individual_name}"
        
rule join_multisample_tpm:
    input:
        expand(config["dir_out"] + "/rna_mapping/{sample}_tpm.txt", sample=SAMPLES),
    output:
        config["dir_out"] + "/rna_mapping/tpm_join_" + config["individual_name"] + ".txt"
    shell:
        "../scripts/multijoin.sh {input} > {output}"
        