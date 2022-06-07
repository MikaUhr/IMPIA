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
conda create -n IMPIA -f IMPIA_envs.yaml
```

4. Download databases
The following data is required.
- Adapter sequence used for adapter removal in [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- rRNA sequence used in [SortMeRNA](https://bioinfo.lifl.fr/RNA/sortmerna/)
- Host genome used for reference filtering
- [COG database](https://www.ncbi.nlm.nih.gov/research/cog-project/) used for gene annotation

### Step2: Configure workflow
Adjust the file `config.yaml` to your setting.

### Step3: Execute workflow
Test your configuration by performing a dry-run via
```
snakemake --use-conda --conda-frontend conda --configfile config.yaml --dry-run
```

Execute the workflow with
```
snakemake --use-conda --conda-frontend conda --configfile config.yaml --core 4
```


## Overview of IMPIA steps

### Step1: Metagenome preprocessing

1. Trimming
Trimming is performed by Trimmomatic.  Trimmomatic trimming  
3. Quality filtering
4. Reference filtering


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
2. covariation analysis

## Overview of IMPIA output

