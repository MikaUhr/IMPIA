#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import csv
import sys

path_df1 = sys.argv[1]
path_df2 = sys.argv[2]
path_out = sys.argv[3]

df1 = pd.read_csv(path_df1, sep = "\t", header = 0)
df2 = pd.read_csv(path_df2, sep = "\t", header = 0)

df_merge = pd.merge(df1, df2, how = 'outer').fillna(0)

df_merge.to_csv(path_out, sep="\t", index=False)