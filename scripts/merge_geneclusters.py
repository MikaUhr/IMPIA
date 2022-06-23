#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import csv
import sys

path_ugc = sys.argv[1]
path_tpm = sys.argv[2]
path_anotation = sys.argv[3]
individual_name = sys.argv[4]
path_out = sys.argv[5]

df_ugc = pd.read_csv(path_ugc, sep="\t", header=None, names=['cluster_id', 'gene_id'])
df_tpm = pd.read_csv(path_tpm, sep=" ", header=0)
df_annotation = pd.read_csv(path_anotation, sep="\t", header=None, names=['gene_id', 'cluster_id'])

# extract individual_name
df_ugc = df_ugc[df_ugc['gene_id'].str.contains(individual_name+"_")]
# Concat unknown-gene-cluster and COG
df_concat = pd.concat([df_annotation, df_ugc], axis=0, sort=False)
df_tpm = pd.merge(df_concat, df_tpm, on='gene_id').sort_values(by = 'cluster_id', axis=0)
# Sum of tpm per COG
df_tpm = df_tpm.groupby('cluster_id').sum()

df_tpm.to_csv(path_out, sep = "\t")