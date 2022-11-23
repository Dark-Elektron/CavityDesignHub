import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mplcursors
from icecream import ic
rf0 = r'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\E_C3794\Export'
# rf0 = r'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\E_3794_4HC_1FPC_Optimize\Export'
# rf1 = r'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\E_3794_4HC_1FPC_Optimize_Reference\Export'
result_folders = [rf0]#, rf1

for nn, result_folder in enumerate(result_folders):
    c0 = 299792458
    mu0 = 4 * np.pi * 1e-7
    imp_cal = 1
    freq_all = pd.read_csv(fr"{result_folder}\Frequency (Multiple Modes).txt", sep='\t', header=None)
    mode_num = np.shape(freq_all)[0]
    freq_all = freq_all.iloc[:, 1] * 1e6  # second column is frequency: 1e6 for MHz and 1e9 for GHz

    VT_x, VT_y = [], []
    for i, freq in enumerate(freq_all):
        Ex = pd.read_csv(fr"{result_folder}\e_X (Z)_Mode {i + 1}.txt", sep='\t', header=None)
        Ey = pd.read_csv(fr"{result_folder}\e_Y (Z)_Mode {i + 1}.txt", sep='\t', header=None)
        Hx = pd.read_csv(fr"{result_folder}\h_X (Z)_Mode {i + 1}.txt", sep='\t', header=None)
        Hy = pd.read_csv(fr"{result_folder}\h_Y (Z)_Mode {i + 1}.txt", sep='\t', header=None)

        Z_axis = Ex.iloc[:, 0] * 1e-3  # 1e-3 to convert mm to m
        # ic(Z_axis)
        # ic(Ex)
        vt_x = np.trapz((Ex.iloc[:, 1] + 1j * Ex.iloc[:, 2] - c0 * mu0 * (Hy.iloc[:, 1] + 1j * Hy.iloc[:, 2])) * np.exp(
            1j * 2 * np.pi * freq * Z_axis / c0), Z_axis)
        vt_y = np.trapz((Ey.iloc[:, 1] + 1j * Ey.iloc[:, 2] + c0 * mu0 * (Hx.iloc[:, 1] + 1j * Hx.iloc[:, 2])) * np.exp(
            1j * 2 * np.pi * freq * Z_axis / c0), Z_axis)

        VT_x.append(vt_x)
        VT_y.append(vt_y)

    R_Q_x = abs(np.array(VT_x) ** 2) / (2 * np.pi * freq_all)
    R_Q_y = abs(np.array(VT_y) ** 2) / (2 * np.pi * freq_all)
    R_Q_T = R_Q_x + R_Q_y
    ic(R_Q_x)
    pp = plt.scatter(freq_all, VT_x, label=f'{nn}', zorder=10, edgecolor='k')
    mplcursors.cursor(pp)

    Qext_all = pd.read_csv(fr"{result_folder}\Q-Factor (lossy E) (Multiple Modes).txt", sep='\t', header=None)
    Qext_all = Qext_all.iloc[:, 1]
    Z_T = Qext_all * R_Q_T * 0.5 * 2 * np.pi * freq_all / c0
    Z_export = [freq_all, R_Q_T, Qext_all, Z_T]

    Z_x = Qext_all * R_Q_x * 0.5 * 2 * np.pi * freq_all / c0
    Z_y = Qext_all * R_Q_y * 0.5 * 2 * np.pi * freq_all / c0

    fs = np.linspace(410e6, 700e6, 2000)  # the frequency interval to calculate impedance from eigenmodes

    Z_sum_x, Z_sum_y = [], []
    for i, freq in enumerate(freq_all):
        z_sum_x = 1 * Z_x[i] / (1 + 1j * Qext_all[i] * (fs / freq - freq / fs))
        z_sum_y = 1 * Z_y[i] / (1 + 1j * Qext_all[i] * (fs / freq - freq / fs))

        Z_sum_x.append(z_sum_x)
        Z_sum_y.append(z_sum_y)

    Z_x_all = sum(Z_sum_x)*1e-3  # impedance in x direction, kOhm/m
    Z_y_all = sum(Z_sum_y)*1e-3  # impedance y direction, kOhm/m
    Z_avg = 0.5*(Z_x_all + Z_y_all)  # impedance y direction
    # print(Z_x_all)
    # po1 = plt.plot(fs, abs(Z_x_all), label='Z_x_all {nn}')
    # po2 = plt.plot(fs, abs(Z_y_all), label='Z_y_all {nn}')
    po3 = plt.plot(fs, abs(Z_avg), label=f'Z_avg {nn}')

    # mplcursors.cursor(po1)
    # mplcursors.cursor(po2)
    mplcursors.cursor(po3)
    # plt.scatter(freq_all, abs(Z_x))
    # plt.scatter(freq_all, abs(Z_y))
plt.legend()
plt.yscale('log')
plt.show()
