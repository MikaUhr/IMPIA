#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import csv
import sys


path_length = sys.argv[1]
path_read_count = sys.argv[2]
path_tpm = sys.argv[3]
sample_name = sys.argv[4]
individual_name = sys.argv[5]

prefix_gene = individual_name + "_gene_"
header = individual_name + "_" + sample_name

df_length = pd.read_csv(path_length, sep="\t", header=None, names = ["gene_id", "length"])
df_read_count = pd.read_csv(path_read_count, sep="\t", header=None, names = ["gene_id", "read_count"])

# Remove the special counters, which count reads that were not counted for any feature
df_read_count = df_read_count[:-5]

# Merge dataframes
df_length.loc[:, "gene_id"] = df_length["gene_id"].str.replace("gene_id=", "")
df_read_count = pd.merge(df_length, df_read_count, on = "gene_id")
df_read_count.loc[:, "gene_id"] = prefix_gene + df_read_count["gene_id"]

# read counts to TPM
read_count = df_read_count.loc[:,'read_count']
gene_len = df_read_count.loc[:, 'length']
norm_by_genelength = read_count.values / gene_len.values
tpm = (norm_by_genelength * 1e6 / np.sum(norm_by_genelength, axis=0))
df_read_count.loc[:, header] = tpm
df_read_count = df_read_count.drop(['length', 'read_count'], axis=1)

df_read_count.to_csv(path_tpm, sep="\t", index=False)
