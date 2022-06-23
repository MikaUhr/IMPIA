# IMPIA
Integrated Metagenome-metatranscriptome Pipeline for Intra-intestinal Analysis (IMPIA) reconstructs a common reference metagenome across multiple intestinal sites from metagenome read data, allowing for identifying the same genes across sites. 
This approach allows gene expression levels to be compared among multi-site samples including those of unknown genes from metatranscriptome data. 
In addition, spatial covariance analysis predicts the function of unknown genes.

## Usage
### Step1: Installation
To install IMPIA, you need conda.

 1. Clone this repository.
```
git clone https://github.com/MikaUhr/IMPIA.git ./IMPIA
```

 2. Change into the IMPIA directory.
```
cd IMPIA
```

 3.  Create the environment to be execute this pipeline.
```
conda create -n IMPIA -f IMPIA_env.yml
```

4. Download the following data:
- Adapter sequence used for adapter removal in [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- rRNA sequence used in [SortMeRNA](https://bioinfo.lifl.fr/RNA/sortmerna/)
  - rfam-5.8s-database-id98.fasta
  - rfam-5s-database-id98.fasta
  - silva-arc-16s-id95.fasta
  - silva-arc-23s-id98.fasta
  - silva-bac-16s-id90.fasta
  - silva-bac-23s-id98.fasta
  - silva-euk-18s-id95.fasta
  - silva-euk-28s-id98.fasta 
- Host genome fasta used for reference filtering
- [COG database](https://www.ncbi.nlm.nih.gov/research/cog-project/) used for gene annotation

5. Install [MetaGeneMark v3.38](http://topaz.gatech.edu/genemark/license_download.cgi) not included in bioconda and set the path.


### Step2: Configure
Adjust the four config files to your setting, e.g.:

1. `config_processing.yml` 
- Raw read directory
  - dir_raw_metagenome - set the directory where the metagenomic raw read data is stored. Performs pre-processing on all fastq files in this directory. The fastq file names should be '{sample}_1.fastq' and '{sample}_2.fastq'. 
  - dir_raw_metatranscriptome - set the directory where the metatranscriptomic raw read data is stored. Performs pre-processing on all fastq files in this directory. The fastq file names should be '{sample}_1.fastq' and '{sample}_2.fastq'.
- Output directory
  - dir_out - set output directory.
- Memory and core
  - memory_per_core - set the size of the RAM of one core of your nodes (GB) (default: 8).
  - threads - set this to the maximum number of cores you want to be using in a run (default: 4).
- Tool parameters
  - megahit_param - the assembly parameter in MEGAHIT
    - k_min - minimum kmer size (<= 255), must be odd number (default: 21).
    - k_max - maximum kmer size (<= 255), must be odd number (default: 141).
    - k_step - increment of kmer size of each iteration (<= 28), must be even number (default: 12).
    - prune_depth - remove unitigs with avg kmer depth less than this value (default: 20).
  - sortmerna_param - the rRNA filtering parameter in SortMeRNA
    - evalue - maximum e-value to report alignments (default: 1e-10).
- Database
  - adapter - set the fasta file of adapter sequence to be trimmed by Trimmomatic.
  - dir_sortmerna - set the directory containing rRNA fasta files.
  - host_genome
    - fasta - set the host genome you want to remove from the read data.
    - index_bowtie2 - set the Bowtie2 index of the host genome.
    - index_hisat2 - set the HISAT2 index of the host genome.

2. `config_merge.yml`
- Contig 
  - contig1 - set the first contig file to be merged.
  - contig2 - set the second contig file to be merged.
- Output directory
  - dir_out - set output directory.
- Tool parameter
  - merge - the merge parameter in quickmerge
    - hco - the quickmerge hco parameter (default: 50)
    - c - the quickmerge c parameter (default: 50)
    - l - minimum seed contig length to be merged (default: 1000)
    - mlength - set the merging length cutoff necessary for use in quickmerge (default: 1000)

3. `config_gene_expression.yml`
- Individual
  - individual_name - set the individual name.
- Reconstructed metagenome
  - reconstructed_metagenome
    - fasta - set the fasta file of metagenomic sequences reconstructed by merging.
    - index_bowtie2 - set the Bowtie2 index of the reconstructed metagenome.
- Output directory
  - dir_out - set output directory.
- Memory and core
  - memory_per_core - set the size of the RAM of one core of your nodes (GB) (default: 8).
  - threads - set this to the maximum number of cores you want to be using in a run (default: 4).
- Tool parameter
  - gene_annotation - diamond options
    - evalue - maximum e-value to report alignments (defaul: 1e-10)
    - query_cover - minimum query cover% to report an alignment (default: 85)
    - min_orf - ignore translated sequences without an open reading frame of at least this length (default: 100)
    - directory_tmp - set the directory for temporary files (defauld: "/dev/shm").
- Database
  - COG
    - COG_fasta - set the fasta file of COG database contains RefSeq accession codes for all proteins with assigned COG domains.
    - COG_csv - set the csv file of COG database. Contains list of orthology domains. Comma-delimited, format: `<domain-id>, <genome-name>, <protein-id>,<protein-length>, <domain-start>, <domain-end>, <COG-id>, <membership-class>,`
    - COG_diamond - set output path of the DIAMOND database to be built from the FASTA file.

4. `config_predict_genefunctions.yml`
- Set directory -
  - dir_input - set input directory. Move all three of the following files (All (the three files) x (the number of individuals)) to this directory: `{dir_our}/gene_annotation/annotation_cogid_{individual}.tsv`, `{dir_our}/gene_annotation/unknown_aminoacid_{individual}.fasta`, `{dir_our}/rna_mapping/tpm_join_{individual}.txt`. These files are included under the output directory `dir_out` configured in `config_gene_expression.yml`.
  - dir_output - set output directory.
- Memory and core
  - memory_per_core - set the size of the RAM of one core of your nodes (GB) (default: 8).
  - threads - set this to the maximum number of cores you want to be using in a run (default: 4).
- Model parameter 
  - gene_clustering - the gene clustering parameter in MMseqs2
    - cluster_mode - 0: Setcover, 1: connected component, 2: Greedy clustering by sequence length  3: Greedy clustering by sequence length (low mem) (default: 2)
    - cov_mode - 0: coverage of query and target, 1: coverage of target, 2: coverage of query 3: target seq. length needs be at least x% of query length, 4: query seq. length needs be at least x% of target length (default: 1)
    - c - list matches above this fraction of aligned (covered) residues (default: 0.9)
    - s - sensitivity will be automatically determined but can be adjusted (default: 7)
    - kmer - kmer per sequence (default: 20)
  - spatial_covariation - the spatial covariation model calculate the bivariate spatial association measure (L statistic) of expression levels between all gene cluster pairs.
    - pseudo_expression - this pseudo-expression value is added when the gene expression level is 0 (default: 0.001). If this value is set to 0, gene clusters containing samples with zero expression levels are removed.
    - variance - gene clusters where the variance of expression levels between sites is greater than this value are used in the analysis (default: 1.0).
    - b - the distance friction coefficient of weight to calculate L statistic (default: 2.0).
    - distance - the distance of the intestinal sites. Distance should be entered with the most proximal position as 0, separated by commas, e,g: "0,10,15" (do not include spaces).
    - L_threshold - when the bivariate spatial association measure (L statistic) exceeds this threshold, the gene cluster pair is defined as a covariate gene (default: 0.885)


### Step3: Execute the workflow

Note that a dry-run must be performed with `--dry-run` to test the configuration before running each workflow.


Execute the workflows from metagenome preprocessing to scaffolding with
```
snakemake -s rules/metagenome_processing.smk --use-conda --conda-frontend conda --configfile config/config_processing.yml --core 4
```

Execute the workflows for metatranscriptome preprocessing with
```
snakemake -s rules/metatranscriptome_processing.smk --use-conda --conda-frontend conda --configfile config/config_processing.yml --core 4
```

Execute the workflows for merging with
```
snakemake -s rules/merge.smk --use-conda --conda-frontend conda --configfile config/config_merge.yml --core 4
```

Execute the workflows from quantification of gene expression level to prediction of unknown function genes by covariation analysis with
```
snakemake -s rules/gene_expression.smk --use-conda --conda-frontend conda --configfile config/config_gene_expression.yml --core 4
```

Execute the workflows for predicting the functions of unknown genes by covariation analysis incorporating bivariate spatial relevance with
```
snakemake -s rules/predict_genefunctions.smk --use-conda --conda-frontend conda --configfile config/config_predict_genefunctions.yml --core 4
```


## Overview of IMPIA workflows

### Metagenome processing

1. Trimming - 
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) removes the defined adapter sequence and trim low-quality regions.
2. Quality filtering - 
[FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) eliminates low quality reads. 
3. Reference filtering - 
Potential host contaminants are filtered by removing reads with sequences aligned to the defined host genome using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
4. Assembly - 
The metagenome reads are assembled by [MEGAHIT](https://www.metagenomics.wiki/tools/assembly/megahit).
5. Scaffolding - 
The contigs are scaffolded by [OPERA-LG](https://sourceforge.net/projects/operasf/) using paired-end read information.

### Metatranscriptome processing

1. Trimming - 
Same as the metagenome processing.
2. Quality filtering - 
Same as the metagenome processing.
3. rRNA filtering - 
[SortMeRNA](https://bioinfo.lifl.fr/RNA/sortmerna/) filter rRNA reads from metatranscriptome data by aligning reads against rRNA database.
4. Reference filtering - 
Potential host contaminants are filtered by removing reads with sequences aligned to the defined host genome using [hisat2](http://daehwankimlab.github.io/hisat2/).

### Merging scaffolds between sites
the scaffolds are merged among sites using [quickmerge](https://github.com/mahulchak/quickmerge).

### Quantification of gene expression levels
1. Gene prediction - 
Gene-coding regions are predicted in the metagenomic sequences by [MetaGeneMark](http://topaz.gatech.edu/genemark/license_download.cgi).
2. Gene Annotation - 
The predicted genes are annotated according to orthologous groups in [COG database](https://www.ncbi.nlm.nih.gov/research/cog-project/) using [DIAMOND](https://github.com/bbuchfink/diamond).
3. Mapping RNA reads to the reconstructed metagenome - 
mRNA reads are mapped to the metagenomic reference sequences by [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
4. Quantification of gene expression levels - 
The number of mRNA reads was counted by [HTSeq](https://htseq.readthedocs.io/en/master/) to quantify the gene expression level. 
The read counts are converted to TPM by `scripts/TPM.py`.
Namely, The read counts are first normalized by the gene length (per kilobase), and then gene-length normalized values are divided by the sum of the gene-length normalized values and multiplied by 10^6.


### Covariation analysis to predict unknown gene functions
1. Generate gene clustering - 
For known genes, gene clusters are generated by COG ID.
The Unknown genes are grouped into gene clusters by protein sequence similarity using [MMSEQS2](https://github.com/soedinglab/MMseqs2).
2. Covariation analysis - 
The functions of unknown genes is predicted by covariation analysis incorporating spatial relevance.
The bivariate spatial association measure (L statistic value) is calculated to detect covarying gene pairs in expression levels in samples. 
This analysis is based on the assumption that functionally similar genes are covariant in their expression levels.
