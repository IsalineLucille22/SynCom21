#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 07:44:16 2024

@author: iguex
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import io
import seaborn as sns
import os
from sklearn.model_selection import train_test_split
import joblib
from sklearn.preprocessing import QuantileTransformer
#from hgboost import hgboost
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, classification_report
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.preprocessing import label_binarize
from sklearn.preprocessing import LabelEncoder
from scipy.integrate import odeint #To sole ODEs system using odeint function

directory = "/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data/"

df = pd.read_excel('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data/SynCom_growthrates_abundancesv2.xlsx', sheet_name = 'community growth soil')

data_name = df.iloc[2:23, 18].values
data_int = df.iloc[2:23, 19:21]
data_val = data_int.values
Sorted_name = df.iloc[2:22, 17].values
Sorted_name_Lysobacter = df.iloc[2:23, 21].values

fig, ax = plt.subplots()
len_specie = range(len(data_name))

label = ['2 months', '6 months']
bottom = [0]*len(label)
for i in len_specie:
    ax.bar(label, data_val[i,:],  bottom=bottom, label=data_name[i])
    bottom += data_val[i,:]


death_rate_vect = np.zeros([len(Sorted_name), 1])
int_time = 4*30*24
N_iter = range(0,len(data_val) - 1)
for i in N_iter:
    ind = np.where(Sorted_name == data_name[i])[0]
    if data_val[i,0] != 0:
        death_rate_vect[ind] = np.log(data_val[i,1]/data_val[i,0])/int_time
        
        
death_rate_vect_Lysobacter = np.zeros([len(Sorted_name_Lysobacter), 1])
int_time = 4*30*24
N_iter = range(0,len(data_val))
for i in N_iter:
    ind = np.where(Sorted_name_Lysobacter == data_name[i])[0]
    if data_val[i,0] != 0:
        death_rate_vect_Lysobacter[ind] = np.log(data_val[i,1]/data_val[i,0])/int_time
    