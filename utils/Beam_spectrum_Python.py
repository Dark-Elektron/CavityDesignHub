# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 11:36:38 2022

@author: shgorgiz
"""

import numpy as np
import matplotlib
from scipy.fft import fft, ifft
import matplotlib.pyplot as plt
import pickle

import scipy.io

# from scipy.io import savemat

# %reset -f
# Some constants
from scipy.io import savemat

eps0 = 8.854187817e-12
mue0 = 4 * np.pi * 1e-7
c0 = 1 / np.sqrt(eps0 * mue0)

# beam properties Z
circumference = 91.1e3
bunch_length = 14.5e-3  # Meter
bunch_length_s = bunch_length / c0  # bunch length in second
beam_current = 1.28  # Ampere
Bunch_charge = 2.43e11 * 1.60217662 * 1e-19 * 1  # Columb
Bunch_beam = 10000  # bunches per beam
Bunch_spacing = circumference / c0 / Bunch_beam  # average bunch spacing based on harmonic number=bucket_fill*dt_gap
Bunch_spacing_table = Bunch_spacing
current_scale = Bunch_spacing / Bunch_spacing_table
bucket_fill = 8

# beam properties H

# bunch_length=5.3e-3; #Meter
# bunch_length_s=bunch_length/c0;#bunch length in second
# beam_current=0.029; #Ampere
# Bunch_charge=1.8e11*1.60217662*1e-19*1; #Columb
# Bunch_spacing=994e-9*1;#average bunch spacing based on harmonic number=bucket_fill*dt_gap
# Bunch_spacing_table=994e-9*1;
# current_scale=Bunch_spacing/Bunch_spacing_table;
# Bunch_beam=328/1;#bunches per beam
# bucket_fill=796;#


# beam properties tt

# bunch_length=2.54e-3; #Meter
# bunch_length_s=bunch_length/c0;#bunch length in second
# beam_current=0.0054; #Ampere
# Bunch_charge=2.3e11*1.60217662*1e-19*1; #Columb
# Bunch_spacing=1*3396e-9*1;#average bunch spacing based on harmonic number=bucket_fill*dt_gap
# Bunch_spacing_table=1*3396e-9*1;
# current_scale=Bunch_spacing/Bunch_spacing_table;
# Bunch_beam=2*48/1;#bunches per beam
# bucket_fill=5440;#

# Periodicity
T = circumference / c0

# Table stores charge and length of the bunch
# odd buckets => long bunches
# even buckets => short bunches

# N_buckets = round(1*400.79e6/(299792458/97.756e3));
N_buckets = round(1 * 400.79e6 / (c0 / circumference))

# Gap size of the bunch
dt_gap = T / N_buckets  # Bucket size
dt_gap = round(dt_gap, 12)
# time scale
t_step = 5e-12
t = np.arange(0, Bunch_spacing, t_step)

N = t.shape[0]

# # %%-----------------------------------------------------------------------%%
# # % Beam filling
# # %%-----------------------------------------------------------------------%%
q_beam = np.zeros((1, N))
q = Bunch_charge * current_scale
sig_z = bunch_length
t_shift = dt_gap - dt_gap / 2
q_beam = q_beam + q / np.sqrt(2 * np.pi) / sig_z * np.exp(-np.power(c0 * (t - t_shift), 2) / (2 * sig_z ** 2))

# Bessy VSR beam current
i_beam = c0 * q_beam

# % Full beam spectrum
i_beam_full = np.kron(np.ones((1, Bunch_beam - 1)), i_beam[0, 1:])
i_beam_full = np.concatenate((i_beam, i_beam_full), axis=1)
N_full = i_beam_full.shape[1]
t_full = np.arange(0, N_full * t_step, t_step)

# # % Determine Fouriertransformation of periodic signals
I_beam = fft(i_beam_full) / N_full;
I_beam[1:] = 1 * I_beam[1:]
f_beam = 1 / (t_full[-1] - t_full[0]) * ((np.arange(1, N_full + 1, 1)) - 1);

f_harm = np.arange(0, 3e9, 1 / T)

spec_BESSY_VSR_beam = {"f": f_harm, "I": I_beam[0, 0:f_harm.shape[0]]}

# plt.stem(spec_BESSY_VSR_beam["f"],np.abs(spec_BESSY_VSR_beam["I"]))


savemat("Beam_Spectrum_save.mat", spec_BESSY_VSR_beam)
# spec_BESSY_VSR_beam = scipy.io.loadmat('Beam_Spectrum_save.mat')
