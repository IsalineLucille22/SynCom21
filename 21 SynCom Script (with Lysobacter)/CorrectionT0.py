#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 13:21:33 2025

@author: iguex
"""

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from scipy.optimize import least_squares #Least squares regression
from scipy.integrate import solve_ivp #Solve ODEs when initial point is known
from scipy.optimize import curve_fit


path = "/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/21 SynCom Script (with Lysobacter)/Data/"
name_file = path + "SSC21_genera_relative-abundances.xlsx"
Data_sheet = "Corr_Absolute"

def remove_last(text):
    text = float(text[:-2])
    return(text)

def linear_curve(a, b, x):
    y = a*x + b
    return(y)

nb_rep = 4
nb_time_point = 2
part_select = nb_rep*nb_time_point + 1

data = pd.read_excel(name_file, "Senka abundance")#pd.read_excel(name_file, 'Absolute_abundances')
time_init = data.columns[1:]
Names = data.time_days
data_part = data.iloc[:, :part_select] #data.iloc[25:46, :4] #
data_copy = pd.read_excel(name_file, Data_sheet)
data = data.iloc[:, :part_select] #data.iloc[25:46, :7]

nb_species = np.shape(data)[0] #Update it when absolute abundance
#time_init = data_part.columns[1:] #data.iloc[nb_species - 2, 1:np.shape(data)[1]]
time_step = list(map(remove_last, time_init)) #np.unique(list(map(remove_last, time_init)))#time_init.apply(remove_last) #pd.apply(remove_last, data.iloc[nb_species - 2, :])
time_step = np.array(time_step, dtype=np.float64).reshape(-1, 1) 
time_step_plt = time_step.flatten()
slopes = np.zeros((nb_species, 1))
intercepts = np.zeros((nb_species, 1))
time_plot = np.arange(0, 170, 1)
time_init_tot = data.columns[1:]
time_step_tot = list(map(remove_last, time_init_tot)) #time_init.apply(remove_last) #pd.apply(remove_last, data.iloc[nb_species - 2, :])
time_step_tot = np.array(time_step_tot, dtype=np.float64).reshape(-1, 1) 
time_plot_tot = np.arange(0, 500, 1)
intercepts_log = np.zeros((nb_species, 1))
for i in range(0, nb_species):
    temp = data_part.iloc[i,:].values
    temp_numeric = pd.to_numeric(temp[1:], errors = 'coerce')
    temp_numeric[temp_numeric == 0] = 1e-10
    log_temp = np.log(temp_numeric)
    min_first_rep = log_temp[0]#np.mean(log_temp[0:5])
    model = LinearRegression()
    model.fit(time_step[:nb_time_point*nb_rep,:], log_temp) 
    slopes[i, 0] = model.coef_[0]
    intercepts[i, 0] = np.min([model.intercept_, min_first_rep])
    y_val = linear_curve(model.coef_[0], model.intercept_, time_plot)
    plt.figure()
    plt.title(Names[i])
    plt.plot(time_plot, y_val, 'r')
    plt.scatter(time_step[:nb_time_point*nb_rep,:].flatten(), log_temp, color = 'b')
    plt.show()
    

Absolute_intercepts = 10*np.exp(intercepts) 
df_intercepts = pd.DataFrame(intercepts)
df_intercepts.insert(0, "Names", Names)
df_intercepts.insert(2, "Absolute_0", Absolute_intercepts)
df_intercepts.rename(columns = {0 : "log_abundances"}, inplace = True)

df_slopes = pd.DataFrame(slopes)
df_slopes.insert(0, "Names", Names)
df_slopes.rename(columns = {0 : "Exp_mu"}, inplace = True)
print(Absolute_intercepts)

with pd.ExcelWriter("Data/Senka_Corr_Initial_abundances.xlsx.xlsx") as writer:
    df_intercepts.to_excel(writer, sheet_name = "Intercepts", index = False)
    df_slopes.to_excel(writer, sheet_name = "Slopes", index = False)
    
data_copy['0_1'] = df_intercepts.iloc[:, 2].values  # or adjust depending on your DataFrame structure


