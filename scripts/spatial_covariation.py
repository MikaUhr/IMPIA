#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import csv
import sys
from argparse import ArgumentParser
import statistics as stt
import math
import pysal as ps
from scipy.spatial import distance
from scipy.spatial.distance import cdist
import scipy.sparse
import multiprocessing as mp

#spatial correlation -> spcor

#Input
path_tpm = sys.argv[1]
path_spcor = sys.argv[2]
pseudo_expression = sys.argv[3]
variance = sys.argv[4]
b = sys.argv[5]
distance_x = sys.argv[6].split(',')
L_threshold = sys.argv[7]
df_tpm =  pd.read_csv(path_tpm, sep = "\t", header=0, index_col=0)

def remove_any_zero_row(df):
    df = df.copy()
    for row in df.index:
        if (df.loc[row] == 0).any():
            df.drop(row, axis=0, inplace=True)
    return df

def weight_matrix(num_col, distance_x, b):
    distance_y =  [0 for n in range(num_col)]
    distance_x_tmp = distance_x
    for i in range(num_individual-1):
        distance_x =  distance_x + distance_x_tmp
    distance = np.array(list(zip(distance_x, distance_y)))
    weight_matrix = cdist(distance, distance, metric='euclidean') ** (-b) 
    weight_matrix[np.isinf(weight_matrix)] = 1
    values = weight_matrix.reshape(num_col * num_col)[:-1].reshape(num_col - 1, num_col + 1)[:, 1:].reshape(num_col * (num_col-1))

    rows = np.array([], dtype=np.uint64)
    cols = np.array([], dtype=np.uint64)
    for i in range(num_col):
        for j in range(num_col):
            #　if i==j, 0
            if i != j:
                rows = np.append(rows, i)
                cols = np.append(cols, j)

    overap = np.unique(rows)
    overap = overap.tolist()
    rows_l = rows.tolist()
    values_l = values.tolist()

    for i in overap:
        location_i = [z for z, x in enumerate(rows_l) if x==i] 
        sumvalue_l = []
        for loci in location_i:   
            value = values_l[loci]
            sumvalue_l.append(value)
        sumvalue = sum(sumvalue_l)
        if sumvalue != 0:
            for loci in location_i:
                value = values_l[loci]
                values_l[loci] = value / sumvalue

    rows = np.array(rows_l)
    values = np.array(values_l)
    values = values.astype(np.float16)
    
    # convert to sparse matrix(csr)
    sparse = scipy.sparse.csr_matrix((values, (rows, cols)),shape=(num_col,num_col))
    wsp = ps.weights.WSP(sparse)
    w = ps.weights.WSP2W(wsp)
    w_array = w.full()[0]
    return w_array

# Check number of columns, sites, individuals
num_col = len(df_tpm.columns)
num_sites = len(df_tpm.columns[df_tpm.columns.str.contains('individual1')])
num_individual = int(num_col/num_sites)

# log2 transformed
df_tpm = df_tpm.replace(0,pseudo_expression)
df_tpm = remove_any_zero_row(df_tpm)
df_tpm = df_tpm.applymap(np.log2) #convert log2

#Extract gene clusters where the variance among sites > threshold
for column_st, column_end in zip(range(0,num_col,num_sites), range(2,num_col,num_sites)):
    df_tpm = df_tpm[df_tpm.iloc[:,column_st:column_end+1].var(axis=1)>variance] 

for column_st, column_end in zip(range(0,num_col,num_sites), range(2,num_col,num_sites)):
    stdval = df_tpm.iloc[:,column_st:column_end+1].std(axis=1)
    for column in range(column_st,column_end+1):
        df_tpm.iloc[:,column] = df_tpm.iloc[:,column]/stdval

# Generate spatial weight matrix
w_array = weight_matrix(num_col, distance_x, b)      

#Compute SSS and subtraction for all geneclusters
df_tpm.loc[:,'avg'] = df_tpm.iloc[:,0:num_col].mean(axis=1)
num_genecluster = len(df_tpm)
sbt_wx = [[0 for i in range(num_col)] for j in range(num_genecluster)] #Sigma(w_ij * x_j) - (Sigma(w_ij * x_j) / n)
sss = [0 for i in range(num_genecluster)] #Initialize SSS
sss_n = [0 for i in range(num_genecluster)] #Initialize SSS_numerator

for genecluster in range(num_genecluster):
    sigma_ij_wx = 0 #Initialize Sigma_i,j(w_ij * x_j)
    sigma_i_wx = [0 for i in range(num_col)] #Initialize Sigma_i(w_ij * x_j)
    avg_wx = 0 #Initialize Sigma_i,j(w_ij * x_j) / n
    sss_d = 0 #Initialize the denominator of SSS_x
    for i in range(num_col):    
        for j in range(num_col):
            sigma_i_wx[i] += (w_array[i][j] * df_tpm.iloc[genecluster, j]) #Sigma_i(w_ij * x_j)
        sigma_ij_wx += sigma_i_wx[i] #Sigma_i,j(w_ij * x_j)

    avg_wx = sigma_ij_wx/num_col #Sigma_i,j(w_ij * x_j) / n

    for i in range(num_col):
        sbt_wx[genecluster][i] = sigma_i_wx[i] - avg_wx
        sss_n[genecluster] += (sbt_wx[genecluster][i]) ** 2 #The numerator of SSS_x
        sss_d += ((df_tpm.iloc[genecluster, i] - df_tpm.iloc[genecluster, num_col])**2) #The denominator of SSS_x
    sss[genecluster] = (sss_n[genecluster] / sss_d)**(1/2)

#L = SSS_x * SSS_y * r
L = [[0 for i in range(num_genecluster)] for j in range(num_genecluster)]
for genecluster1 in range(num_genecluster):
    for genecluster2 in range(num_genecluster):
        if(genecluster2 > genecluster1):
            r_n = 0
            r = 0

            for i in range(num_col):
                r_n += sbt_wx[genecluster1][i] * sbt_wx[genecluster2][i]            

            r = r_n / ((sss_n[genecluster1] * sss_n[genecluster2]) ** (1/2))
            L[genecluster1][genecluster2] = sss[genecluster1] * sss[genecluster2] * r

#Create correlation list
list_spcor = []
for column in range(len(L)):
    for row in range(len(L)):
        if(column < row):
            add = [ df_tpm.index[column], df_tpm.index[row] , L[column][row]]
            list_spcor.append(add)

#Convert list to DataFrame
dt_list_spcor = pd.DataFrame(list_spcor, columns=['ortholog1', 'ortholog2', 'L'])

#Extract gene cluster pairs with L > L_threshold
dt_list_spcor.loc[:,'L'] = dt_list_spcor['L'].astype(float)
dt_list_spcor = dt_list_spcor[abs(dt_list_spcor['L']) > L_threshold]

dt_list_spcor.to_csv(path_spcor, sep = ",", index = False)

