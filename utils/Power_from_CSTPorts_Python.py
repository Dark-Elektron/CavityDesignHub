# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 17:35:26 2022

@author: shgorgiz
"""

import numpy as np
import matplotlib
from numpy import loadtxt
from scipy.interpolate import interp1d
import scipy.io


def get_beam_spectrum(wp):
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
    # import dill  # pip install dill --user

    import scipy.io

    # from scipy.io import savemat

    # %reset -f
    # Some constants
    eps0 = 8.854187817e-12
    mue0 = 4 * np.pi * 1e-7
    c0 = 1 / np.sqrt(eps0 * mue0)

    # machine parameters
    circumference = 91.1e3
    # beam properties Z
    bunch_length = 14.5e-3  # Meter
    bunch_length_s = bunch_length / c0  # bunch length in second
    beam_current = 1.28  # Ampere
    Bunch_charge = 2.43e11 * 1.60217662 * 1e-19 * 1  # Columb
    Bunch_beam = 10000  # bunches per beam
    Bunch_spacing = circumference / c0 / Bunch_beam  # average bunch spacing based on harmonic number=bucket_fill*dt_gap
    Bunch_spacing_table = Bunch_spacing
    current_scale = Bunch_spacing / Bunch_spacing_table
    bucket_fill = 8

    Z = {"freq": 400.79e6,
         "sigma": bunch_length,
         "sigma_s": bunch_length_s,
         "I0": beam_current,
         "Qb": Bunch_charge,
         "bunches per beam": Bunch_beam,
         "bunch spacing": Bunch_spacing,
         "bunch spacing table": Bunch_spacing_table,
         "current scale": current_scale,
         "bucket_fill": bucket_fill}

    # beam properties H
    bunch_length = 5.3e-3  # Meter
    bunch_length_s = bunch_length / c0  # bunch length in second
    beam_current = 0.029  # Ampere
    Bunch_charge = 1.8e11 * 1.60217662 * 1e-19 * 1  # Columb
    Bunch_spacing = 994e-9 * 1  # average bunch spacing based on harmonic number=bucket_fill*dt_gap
    Bunch_spacing_table = 994e-9 * 1
    current_scale = Bunch_spacing / Bunch_spacing_table
    Bunch_beam = int(328 / 1)  # bunches per beam
    bucket_fill = 796  #

    H = {"freq": 400.79e6,
         "sigma": bunch_length,
         "sigma_s": bunch_length_s,
         "I0": beam_current,
         "Qb": Bunch_charge,
         "bunches per beam": Bunch_beam,
         "bunch spacing": Bunch_spacing,
         "bunch spacing table": Bunch_spacing_table,
         "current scale": current_scale,
         "bucket_fill": bucket_fill}

    # beam properties tt

    bunch_length = 2.54e-3  # Meter
    bunch_length_s = bunch_length / c0  # bunch length in second
    beam_current = 0.0054  # Ampere
    Bunch_charge = 2.3e11 * 1.60217662 * 1e-19 * 1  # Columb
    Bunch_spacing = 1 * 3396e-9 * 1  # average bunch spacing based on harmonic number=bucket_fill*dt_gap
    Bunch_spacing_table = 1 * 3396e-9 * 1
    current_scale = Bunch_spacing / Bunch_spacing_table
    Bunch_beam = int(2 * 48 / 1)  # bunches per beam
    bucket_fill = 5440  #

    ttbar = {"freq": 801.58e6,
             "sigma": bunch_length,
             "sigma_s": bunch_length_s,
             "I0": beam_current,
             "Qb": Bunch_charge,
             "bunches per beam": Bunch_beam,
             "bunch spacing": Bunch_spacing,
             "bunch spacing table": Bunch_spacing_table,
             "current scale": current_scale,
             "bucket_fill": bucket_fill}

    wp_dict = {"Z": Z, "H": H, "ttbar": ttbar}

    # Periodicity
    T = circumference / c0

    # Table stores charge and length of the bunch
    # odd buckets => long bunches
    # even buckets => short bunches

    # N_buckets = round(1*400.79e6/(299792458/97.756e3))
    N_buckets = round(1 * wp_dict[wp]['freq'] / (c0 / circumference))

    # Gap size of the bunch
    dt_gap = T / N_buckets  # Bucket size
    dt_gap = round(dt_gap, 12)

    # time scale
    t_step = 5e-12
    t = np.arange(0, wp_dict[wp]["bunch spacing"], t_step)

    N = t.shape[0]

    # # %%-----------------------------------------------------------------------%%
    # # % Beam filling
    # # %%-----------------------------------------------------------------------%%
    q_beam = np.zeros((1, N))
    q = wp_dict[wp]["Qb"] * wp_dict[wp]["current scale"]
    sig_z = wp_dict[wp]["sigma"]
    t_shift = dt_gap - dt_gap / 2
    q_beam = q_beam + q / np.sqrt(2 * np.pi) / sig_z * np.exp(-np.power(c0 * (t - t_shift), 2) / (2 * sig_z ** 2))

    # Bessy VSR beam current
    i_beam = c0 * q_beam

    # % Full beam spectrum
    i_beam_full = np.kron(np.ones((1, wp_dict[wp]["bunches per beam"] - 1)), i_beam[0, 1:])
    i_beam_full = np.concatenate((i_beam, i_beam_full), axis=1)
    N_full = i_beam_full.shape[1]
    t_full = np.arange(0, N_full * t_step, t_step)

    # # % Determine Fourier transformation of periodic signals
    I_beam = fft(i_beam_full) / N_full
    I_beam[1:] = 1 * I_beam[1:]
    f_beam = 1 / (t_full[-1] - t_full[0]) * ((np.arange(1, N_full + 1, 1)) - 1)

    f_harm = np.arange(0, 3e9, 1 / T)

    beam_spectrum = {"f": f_harm, "I": I_beam[0, 0:f_harm.shape[0]]}

    # plt.stem(beam_spectrum["f"],np.abs(beam_spectrum["I"]))

    # savemat("Beam_Spectrum_save.mat", beam_spectrum)
    # beam_spectrum = scipy.io.loadmat('Beam_Spectrum_save.mat')
    return beam_spectrum


result_folder = r'D:\CST Studio\LONGITIDUNAL_WAKE_TESLA_CAVITY_100m\Export'

# beam_spectrum = scipy.io.loadmat('Beam_Spectrum_save.mat')
beam_spectrum = get_beam_spectrum('ttbar')
f_harm = beam_spectrum["f"]
N_port = 20

Charge_dist = loadtxt(f"{result_folder}/Particle Beams_ParticleBeam1_Charge distribution spectrum.txt", comments="#", unpack=False)
Long_Imp = loadtxt(f"{result_folder}/Particle Beams_ParticleBeam1_Wake impedance_Z.txt", comments="#", unpack=False)

Charge_dist_interpolated_function = interp1d(Charge_dist[:, 0], Charge_dist[:, 1], 'cubic')

Charge_dist_interpolated = Charge_dist_interpolated_function(Long_Imp[:, 0])

Port_summation = 0
Port_signal = []
for i in range(N_port):
    filename = f"{result_folder}/Port signals_o{1}({i+1}),pb.txt"
    data = loadtxt(filename, comments="#", unpack=False)
    data[:, 1] = (2 / np.pi) * data[:, 1] / Charge_dist_interpolated ** 2
    Port_signal.append(data)
    Port_summation = Port_summation + abs(data[:, 1])

Long_Imp_interp_f = interp1d(Long_Imp[:, 0], Long_Imp[:, 1], 'cubic')
Long_Imp_interp = Long_Imp_interp_f(f_harm / 1e6)

Port_interp = np.empty((np.size(f_harm, 1), N_port))
power_port = np.empty((np.size(f_harm, 1), N_port))

for i in range(N_port):
    Port_interp_f = interp1d(Long_Imp[:, 0], Port_signal[i][:, 1], 'cubic')
    Port_interp[:, i] = Port_interp_f(f_harm / 1e6)
    power_port[:, i] = 2 * abs(beam_spectrum["I"] ** 2) * abs(Port_interp[:, i])

power_imp = 2 * abs(beam_spectrum["I"] ** 2) * Long_Imp_interp

ind1 = np.argwhere(f_harm > 0.45e9)[0, 1]
ind2 = np.argwhere(f_harm > 1.5e9)[0, 1]

power_imp_distribution = [sum(power_imp[0, ind1:ind2]), 0, 0]
power_imp_distribution[1] = sum(power_imp[0, ind2:])
power_imp_distribution[2] = power_imp_distribution[0] + power_imp_distribution[1] + power_imp_distribution[
    1] * 19.5 / 32.2

power_port_distribution = sum(power_port[ind1:ind2, :])
power_port_distribution = np.c_[power_port_distribution, sum(power_port[ind2:, :])]
power_port_distribution = np.c_[power_port_distribution, (
        power_port_distribution[:, 0]
        + power_port_distribution[:, 1]
        + power_port_distribution[:, 1] * 19.5 / 32.2)]

power_port_all = sum(power_port_distribution)
power_port_percent = 100 * power_port_distribution[:, 2] / power_port_all[2]
