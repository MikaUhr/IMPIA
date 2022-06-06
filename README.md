# IMPIA
Integrated Metagenome-metatranscriptome Pipeline for Intra-intestinal Analysis (IMPIA) reconstructs a common reference metagenome across multiple intestinal sites from metagenome read data, allowing for identifying the same genes across sites. 
This approach allows gene expression levels to be compared among multi-site samples including those of unknown genes from metatranscriptome data. 
In addition, spatial covariance analysis predicts the function of unknown genes.


## Requirement
### Preprocessing for metagenome data

### Preprocessing for metatranscriptome data

### Metagenome reconstruction

### 

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
conda create -n IMPIA -f IMPIA_envs.yaml
```

4. Download databases

### Step2: Configure workflow
1. Adjust the file `config.yaml` to your setting.
2. 

### Step3: Execute workflow
```
snakemake --use-conda --conda-frontend conda --configfile config.yaml --core 4
```


## Overview of IMPIA steps

### Preprocessing

Metagenome data
- Trimmomatic v0.36
- FASTX-Toolkit version 0.0.14 
- cmpfastq_pe
- Bowtie2 version 2.3.4.3

Metatranscriptome data
- Trimmomatic v0.36
- FASTX-Toolkit version 0.0.14 
- cmpfastq_pe
- SortMeRNA version 2.1
- hisat2

### Reconstruct metagenomes


### Quantification of gene expression levels


### Covariation analysis to predict unknown gene functions

