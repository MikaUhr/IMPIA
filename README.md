# IMPIA
Integrated Metagenome-metatranscriptome Pipeline for Intra-intestinal Analysis (IMPIA) reconstructs a common reference metagenome across multiple intestinal sites from metagenome read data, allowing for identifying the same genes across sites. 
This approach allows gene expression levels to be compared among multi-site samples including those of unknown genes from metatranscriptome data. 
In addition, spatial covariance analysis predicts the function of unknown genes.

## Usage
### Step2: Installation
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

4. Download databases

The following data is required.
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

5. Install tools
Install [MetaGeneMark v3.38](http://topaz.gatech.edu/genemark/license_download.cgi) and set the path.


### Step2: Configure
Adjust the four config files to your setting, e.g.:

1. `config_processing.yml` 
- Raw read directory
  - dir_raw_metagenome - set the directory where the metagenomic raw read data is stored. Performs pre-processing on all fastq files in the directory.
  - dir_raw_metatranscriptome - set the directory where the metatranscriptomic raw read data is stored. Performs pre-processing on all fastq files in the directory.
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
  - contig1 - set the first contig to be merged.
  - contig2 - set the second contig to be merged.
- Output directory
  - dir_out - set output directory.
- Tool parameter
  - merge - the merge parameter in Quickmerge
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
  - dir_input - set input directory. Move all three of the following files (All (the three files) x (the number of individuals)) to the input directory: `{dir_our}/gene_annotation/annotation_cogid_{individual}.tsv`, `{dir_our}/gene_annotation/unknown_aminoacid_{individual}.fasta`, `{dir_our}/rna_mapping/tpm_join_{individual}.txt`. These files are included under the output directory `dir_out` configured in `config_gene_expression.yml`.
  - dir_output - set output directory.
- Memory and core
  - memory_per_core - set the size of the RAM of one core of your nodes (GB) (default: 8).
  - threads - set this to the maximum number of cores you want to be using in a run (default: 4).
- Model parameter 
  - gene_clustering - the gene clustering parameter in MMseqs2
    - cluster_mode: 0: Setcover, 1: connected component, 2: Greedy clustering by sequence length  3: Greedy clustering by sequence length (low mem) (default: 2)
    - cov_mode: 0: coverage of query and target, 1: coverage of target, 2: coverage of query 3: target seq. length needs be at least x% of query length, 4: query seq. length needs be at least x% of target length (default: 1)
    - c: list matches above this fraction of aligned (covered) residues (default: 0.9)
    - s: sensitivity will be automatically determined but can be adjusted (default: 7)
    - kmer: kmer per sequence (default: 20)
  - spatial_covariation
    - pseudo_expression: 0.001
    - variance: 1.0
    - b: 2.0
    - distance: ""
    - L_threshold: 0.885


### Step3: Execute the workflow

Note that a dry-run must be performed with `--dry-run` to test the configuration before running each workflow.


Execute the workflows from metagenome preprocessing to scaffolding with
```
snakemake -s metagenome_preprocessing_scaffolding.smk --use-conda --conda-frontend conda --configfile config/config_preprocess.yaml --core 4
```

Execute the workflows for metatranscriptome preprocessing with
```
snakemake -s metatranscriptome_preprocessing.smk --use-conda --conda-frontend conda --configfile config/config_preprocess.yaml --core 4
```

Execute the workflows for merging with
```
snakemake -s merge.smk --use-conda --conda-frontend conda --configfile config/config_merge.yaml --core 4
```

Execute the workflows from quantification of gene expression level to prediction of unknown function genes by covariation analysis with
```
snakemake -s gene_expression.smk --use-conda --conda-frontend conda --configfile config/config_gene_expression.yaml --core 4
```


## Overview of IMPIA steps

### Step1: Metagenome preprocessing

1. Trimming
2. Quality filtering
3. Reference filtering

- Trimmomatic v0.36
- FASTX-Toolkit version 0.0.14 
- cmpfastq_pe
- Bowtie2 version 2.3.4.3


### Step2: Metatranscriptome preprocessing

1. Trimming
2. Quality filtering
3. rRNA filtering
4. Reference filtering

- Trimmomatic v0.36
- FASTX-Toolkit version 0.0.14 
- cmpfastq_pe
- SortMeRNA version 2.1
- hisat2

### Step3: Reconstruct metagenomes

1. Assembly
2. Scaffolding
3. Merge
4. Gene prediction
5. Gene Annotation

### Step4: Quantification of gene expression levels

1. Mapping RNA reads to the reconstructed metagenome
2. Count reads

### Step5: Covariation analysis to predict unknown gene functions

1. Unknown gene clustering
2. Covariation analysis
