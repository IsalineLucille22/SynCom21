#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Logistic growth model with lag phase and extrapolated initial biomass (t=0)
Author: Isaline Guex
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares


# -------------------------------------------------------------------
# --- File paths
# -------------------------------------------------------------------
path = "/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/21 SynCom Script (with Lysobacter)/Data/"
name_file = path + "SSC21_genera_relative-abundances.xlsx"
Data_sheet = "Corr_Absolute"


# -------------------------------------------------------------------
# --- Helper functions
# -------------------------------------------------------------------
def remove_last(text):
    """Convert Excel column headers like '12d' -> 12.0"""
    return float(text[:-2])


def fun_Logistic_Exp_(t, mu_max, K_S, LT, x_0):
    """
    Logistic growth model with lag phase.
    t: time vector (in shifted coordinates)
    mu_max: max growth rate
    K_S: carrying capacity
    LT: lag time (in same units as t)
    x_0: initial biomass (array-like)
    """
    x_0 = x_0[0]

    # Ensure monotonic growth (avoid K_S < x_0)
    # if x_0 >= K_S:
    #     K_S = x_0 * 1.01

    A = (K_S - x_0) / x_0
    #A = max(A, 1e-6)

    x = K_S/(1 + A*np.exp(-mu_max*(t - 0*LT)))
    x[x < 1e-09] = 1e-09
    return x


def residuals(params, t_shift, y_obs, x_0):
    """
    Residuals for least_squares fitting.
    Fit is done in shifted time coordinates (first observation = 0)
    """
    mu_max, K_S, LT = params
    scale = np.max(y_obs)
    if scale == 0:
        return np.zeros_like(y_obs)
    y_obs_scaled = y_obs/scale
    x_0_scaled = np.array(x_0)/scale
    K_S_scaled = K_S/scale
    y_pred_scaled = fun_Logistic_Exp_(t_shift, mu_max, K_S_scaled, LT, x_0_scaled)
    return y_pred_scaled - y_obs_scaled


# -------------------------------------------------------------------
# --- Load data
# -------------------------------------------------------------------
data = pd.read_excel(name_file, "Senka abundance")#pd.read_excel(name_file, "Absolute_abundances")#
data_copy = pd.read_excel(name_file, Data_sheet)
nb_species = np.shape(data)[0]
Names = data.time_days
time_init = data.columns[1:]
time_init = list(map(remove_last, time_init)) 
time_step = np.array(time_init, dtype=np.float64).reshape(-1, 1)
time_step_plt = time_step.flatten()

# -------------------------------------------------------------------
# --- Storage matrices
# -------------------------------------------------------------------
mu_max_mat = np.zeros((nb_species, 1))
KS_mat = np.zeros((nb_species, 1))
LT_mat = np.zeros((nb_species, 1))
x_init = np.zeros((nb_species, 1))  # extrapolated biomass at t = 0

# -------------------------------------------------------------------
# --- Fit per species
# -------------------------------------------------------------------
for i in range(nb_species):
    Data = data.iloc[i, :].values
    Data = pd.to_numeric(Data[1:], errors="coerce")
    name_species = Names[i]
    rep_1 = Data

    diff = -rep_1[0] + rep_1[-1]
    max_diff = max(rep_1) - rep_1[0]
    ind = len(rep_1)  # you can change this index if needed

    if (ind <= 2 or diff < 0) and max_diff < 1e-8:
        mu_max_mat[i, 0] = 0
        KS_mat[i, 0] = 0
        LT_mat[i, 0] = 0
        plt.figure()
        plt.scatter(time_step, rep_1)
        plt.plot(time_step, np.full_like(time_step, rep_1[0]), "r")
        plt.title(name_species)
        plt.show()
    else:
        # --- Prepare data for fitting ---
        time_obs = time_step[0:(ind + 1)].flatten()
        rep_obs = rep_1[0:(ind + 1)]
        #rep_obs[rep_obs == 0] = min(rep_obs[rep_obs > 0])
        x_0 = [rep_obs[0]]

        # Shift time so first observation = 0
        t_shift = time_obs - time_obs[0]
        time_plot_shift = np.linspace(0, t_shift[-1], 300)

        # --- Fit model ---
        param0 = [0.2, 3 * np.mean(rep_obs), 0]  # μ, K_S, LT
        bounds = ([0, 1e-10, 0], [10, 1e-2, 200])

        param_sim = least_squares(
            residuals, x0 = param0, args = (t_shift, rep_obs, x_0), bounds = bounds
        )

        mu_max_temp, KS_temp, LT_temp = param_sim.x
        print(f"{name_species}: μ={mu_max_temp:.3f}, KS={KS_temp:.3e}, LT={LT_temp:.2f}")

        # --- Predict on shifted timeline ---
        y_sim_shift = fun_Logistic_Exp_(time_plot_shift, mu_max_temp, KS_temp, LT_temp, x_0)

        # --- Extrapolate back to real time t = 0 ---
        t_real_full = np.linspace(0, time_obs[-1], 300)
        t_for_extrapolation = t_real_full - time_obs[0]  # negative before first obs
        y_full = fun_Logistic_Exp_(t_for_extrapolation, mu_max_temp, KS_temp, LT_temp, x_0)

        # store extrapolated abundance at t = 0
        x_init[i, 0] = y_full[0]

        # --- Plot ---
        plt.figure()
        plt.scatter(time_obs, rep_obs, label="observed")
        plt.plot(t_real_full, y_full, "r", label="fit + extrapolation")
        plt.axvline(0, color="k", linestyle="--", alpha=0.5)
        plt.title(name_species)
        plt.xlabel("Time (days)")
        plt.ylabel("Biomass (g/mL)")
        plt.legend()
        plt.show()

        mu_max_mat[i, 0] = mu_max_temp
        KS_mat[i, 0] = KS_temp
        LT_mat[i, 0] = LT_temp


# -------------------------------------------------------------------
# --- Save estimated initial abundances
# -------------------------------------------------------------------
print(x_init)
Absolute_0 = 10*x_init
df_intercepts = pd.DataFrame(Absolute_0)
df_intercepts.insert(0, "Names", Names)
df_intercepts.rename(columns={0: "Senka_Corr_Initial_abundances"}, inplace=True)

with pd.ExcelWriter("Data/Senka_Corr_Initial_abundances.xlsx") as writer:
    df_intercepts.to_excel(writer, sheet_name="Initial_Abund", index=False)

