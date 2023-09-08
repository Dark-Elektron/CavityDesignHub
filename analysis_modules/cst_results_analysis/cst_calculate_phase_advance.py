import math
import os

import pandas as pd
import numpy as np
from icecream import ic
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

from analysis_modules.cst_results_analysis.cst_result_compile import monopole_mode, dipole_mode

c0 = 299792458


# Define a cosine square window function
def cosine_square_window(length):
    window = np.cos(np.pi * np.arange(length) / (length - 1)) ** 2
    return window / np.sum(window)


def phase_advance(filepath, cell_length, n_cells):
    fig, axs = plt.subplot_mosaic([[0, 0, 1]], figsize=(15, 4))
    ax_twin = axs[0].twinx()
    PHI = {}
    phi_med = []
    mode_num = []
    # takes in arguments as filepath which is a cst field on axis plot and calculates the phase advane
    for path in os.listdir(filepath):
        if 'e_Z (Z)_' in path:
            mn = int(path.split(' ')[-1].split('.')[0])
            if mn < 500:
                # load file
                df = pd.read_csv(fr"{filepath}/{path}", header=None, sep='\s+', names=['s', 'E'])
                S, E = df['s'], df['E']

                # # check if mode is nonopole mode
                # if not monopole_mode(S, E, cell_length/2, n_cells):
                #     continue
                ic(mn)
                if not dipole_mode(S, E, cell_length/2, n_cells):
                    continue

                # interpolate points
                f = interp1d(S, E, kind='cubic')

                # new x
                N = 1000  # number of samples
                s_new = np.linspace(np.min(S), np.max(S), N)
                mode_num.append(mn)
                # sweep along z
                PHI[mn] = []
                for i, s in enumerate(s_new):
                    if np.min(S) + cell_length < s < np.max(S) - cell_length:
                        El = f(s - cell_length)
                        Er = f(s + cell_length)
                        Em = f(s)

                        phi = np.rad2deg(np.arccos((El + Er) / (2 * Em)))
                    else:
                        phi = 0

                    PHI[mn].append(phi)

                # plt.plot(s_new, PHI[mn])
                # plt.show()
                if n_cells <= 2:
                    ic("it's here now")
                    # fit x^2 to points
                    phi_med.append(fit_x2(s_new[(np.where((s_new > -cell_length) & (s_new < cell_length) & np.logical_not(np.isnan(np.array(PHI[mn])))))],
                           np.array(PHI[mn])[np.where((s_new > -cell_length) & (s_new < cell_length) & np.logical_not(np.isnan(np.array(PHI[mn]))))]))
                else:
                    # fit x^10 to points
                    phi_med.append(fit_x10(s_new[(np.where((s_new > -cell_length) & (s_new < cell_length) & np.logical_not(np.isnan(np.array(PHI[mn])))))],
                           np.array(PHI[mn])[np.where((s_new > -cell_length) & (s_new < cell_length) & np.logical_not(np.isnan(np.array(PHI[mn]))))]))

                # phi_med.append(np.nanmedian(np.array(PHI[mn])[np.where(
                #     (s_new > -cell_length) & (s_new < cell_length))]))  # take median in the neighbourhood of the center

                # make plot
                ax_twin.plot(s_new / cell_length, f(s_new), ls='--', c='k', lw=0.8)

                axs[0].plot(s_new / cell_length, PHI[mn])

    axs[0].set_ylim(0, 180)
    # axs[0].grid()

    # get mode frequencies
    mode_num = np.array(mode_num)
    phi_med = np.array(phi_med)[mode_num.argsort()]

    mode_num = np.sort(mode_num)
    ic(mode_num, phi_med)

    freq = pd.read_csv(fr"{filepath}/Frequency (Multiple Modes).txt", header=None, sep='\s+', names=['mode', 'freq'])
    mode_index = [np.where(freq['mode'] == element)[0].tolist()[0] for element in np.unique(mode_num)]
    ic(mode_index)
    axs[1].scatter(phi_med, np.array(freq['freq'])[mode_index], edgecolor='r', facecolor='none')
    # conf_intervals = [2.5]*len(mode_index)
    # axs[1].errorbar(phi_med, np.array(freq['freq'])[mode_index], xerr=conf_intervals, fmt='o', capsize=5, mec='r', mfc='none')

    # light line
    p = math.ceil(max(freq['freq']) / min(freq['freq']))  # number of light lines
    nn = np.tile(np.array([0, 180]), p)
    nn = [nn[i:i + 2] for i in range(0, len(nn) - 1)]

    op_freq = c0 / (2 * cell_length) * 1e-3
    freq_list = np.arange(0, p + 4) * op_freq
    freq_list = [freq_list[i:i + 2] for i in range(0, len(freq_list))]

    # light line
    for i, x in enumerate(nn):
        # ic(x, freq_list[i])
        axs[1].plot(x, freq_list[i], c='k', label='light line')

    # plt.ylim(min(freq['freq']) - 10, max(freq['freq'])+10)
    axs[1].set_xlim(-3, 183)

    plt.margins(x=0, y=0)
    plt.tight_layout()
    plt.show()


def x10(x, A, B):
    return A + B*x**10


def fit_x10(x, y):
    parameters, cov = curve_fit(x10, x, y)
    plt.plot(x, y)
    x_new = np.linspace(np.min(x), np.max(x), 1000)
    # ic(x10(0, *parameters))
    plt.plot(x_new, x10(x_new, *parameters))
    plt.show()
    return x10(0, *parameters)


def x2(x, A, B, C, D, E):
    return A + B*x + C*x**2 + D*x**3 + E*x**4


def fit_x2(x, y):
    parameters, cov = curve_fit(x2, x, y)
    # plt.plot(x, y)
    # x_new = np.linspace(np.min(x), np.max(x), 1000)
    # ic(x2(0, *parameters))
    # plt.plot(x_new, x2(x_new, *parameters))
    # plt.show()
    return x2(0, *parameters)


def plot_brillouin(op_freq, p, cst_result_folder=None):
    fig, axs = plt.subplots(1, p)

    # get dictionary instead from cst run folder
    d = pd.read_csv(fr"D:\CST Studio\5. tt\Eigenmode\E3795_PBC.csv", sep=",", skipfooter=1, engine='python').to_dict(
        orient='list')
    freq_dict = {key: val for key, val in d.items() if 'Mode' in key}
    phase = d['phase']

    nn = np.tile(np.array([0, 180]), p)
    nn = [nn[i:i + 2] for i in range(0, len(nn) - 1)]

    freq_list = np.arange(0, p + 3) * op_freq
    freq_list = [freq_list[i:i + 2] for i in range(0, len(freq_list))]

    # light line
    for i, x in enumerate(nn):
        for axs[0] in axs:
            axs[0].plot(x, freq_list[i], c='k', label='light line')

    for key, mode_freqs in enumerate(freq_dict.values()):
        for i, axs[0] in enumerate(axs):
            axs[0].plot(np.array(phase), mode_freqs, marker='o', mec='k', lw=3)
            axs[0].set_ylim(min(freq_list[i + 1]) - 50, max(freq_list[i + 1]))
            axs[i].set_xlim(-3, 183)

    plt.tight_layout()
    plt.legend()
    plt.show()


if __name__ == '__main__':
    # phase_advance(fr'D:\CST Studio\3. W\Eigenmode\E_C3794\Export', 2 * 187, 2)
    phase_advance(fr'D:\CST Studio\5. tt\Eigenmode\E_C3795\Export', 2*93.5, 5)
    # phase_advance(fr'D:\CST Studio\5. tt\Eigenmode\E_C3795_dipole\Export', 2 * 93.5, 5)
    # phase_advance(fr'D:\CST Studio\3. W\Eigenmode\E_C3794_390_665\Export', 2 * 187, 2)
    # phase_advance(fr'D:\CST Studio\Youtube Videos\TESLA\Export', 2 * 57.6524, 9)
    # phase_advance(fr'D:\CST Studio\Youtube Videos\TESLA_\Export', 2 * 57.6524, 9)
    # phase_advance(fr'D:\CST Studio\Youtube Videos\TESLA_hex\Export\e_Z (Z).txt', 2 * 57.6524, 9)
