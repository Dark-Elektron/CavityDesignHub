import json

import mplcursors
import numpy as np
import matplotlib.pyplot as plt
from icecream import ic
from matplotlib import ticker
import matplotlib as mpl
from utils.file_reader import FileReader
from modules.plot_module.matplotlib_annotation_objects import DraggableText

fr = FileReader()




def plot_settings():
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16

    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['axes.titlesize'] = 16
    mpl.rcParams['legend.fontsize'] = 14
    mpl.rcParams['legend.title_fontsize'] = 14

    mpl.rcParams['figure.figsize'] = [10, 6]
    mpl.rcParams['figure.dpi'] = 100


plot_settings()


class Cavity:
    def __init__(self, n_cells, l_cell_mid, freq, vrf, R_Q, Q, G, Epk_Eacc, Bpk_Eacc):
        # geometric parameters
        # input
        self.n_cells = n_cells
        self.l_cell_mid = l_cell_mid  # m
        self.v_rf = vrf
        self.R_Q = R_Q
        self.Q = Q
        self.op_freq = freq  # Hz
        self.e = Epk_Eacc
        self.b = Bpk_Eacc
        self.G = G

        # calculated
        self.l_active = 2 * self.n_cells * self.l_cell_mid  # m
        # self.n_cav = self.v_rf / (self.E_acc * self.l_active)
        # self.l_cavity = self.n_cav * self.n_cells * self.l_active + (self.n_cav - 1) * 6 * self.l_active

        # self.E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m

        # cryo parameters
        # self.eta = 1/219  # 400 MHz
        self.eta = 1/745  # 800 MHz

        # power parameters
        # calculated
        self.p_sr = 50e6  # W

    def set_Eacc(self, Eacc):
        self.E_acc = Eacc
        self.n_cav = self.v_rf / (self.E_acc * self.l_active)
        self.l_cavity = self.l_active + 8 * self.l_cell_mid  # Lcavity = Lactive + 2 λ0

        Rs_NbCu_2k_400Mhz = 0.57 * (Eacc * self.b) + 28.4  # nOhm
        Rs_NbCu_4_5k_400Mhz = 39.5 * np.exp(0.014 * (Eacc * self.b)) + 27  # nOhm
        Rs_bulkNb_2k_400Mhz = (2.33 / 1000) * (Eacc * self.b) ** 2 + 26.24  # nOhm
        Rs_bulkNb_4_5k_400Mhz = 0.0123 * (Eacc * self.b) ** 2 + 62.53  # nOhm

        Rs_NbCu_2k_800Mhz = 1.45 * (Eacc * self.b / Eacc) + 92  # nOhm
        Rs_NbCu_4_5k_800Mhz = 50 * np.exp(0.033 * (Eacc * self.b / Eacc)) + 154  # nOhm
        Rs_bulkNb_2k_800Mhz = (16.4 + Eacc * self.b / Eacc * 0.092) * (800 / 704) ** 2  # nOhm
        Rs_bulkNb_4_5k_800Mhz = 4 * (62.7 + (Eacc * self.b / Eacc) ** 2 * 0.012)  # nOhm

        self.Q0 = self.G * 1e9 / Rs_bulkNb_2k_800Mhz  # c

    def get_power(self):
        self.p_in = 2 * self.p_sr / self.n_cav  # maximum synchrotron radiation per beam multipled by two for two beams
        self.p_cryo = 8 / (np.sqrt(self.op_freq / 500e6))
        self.pdyn = (1 / self.eta) * self.v_rf * (self.E_acc * self.l_active / (self.R_Q * self.Q0))
        self.pstat = (1 / self.eta) * (self.l_cavity * self.v_rf / (self.l_active * self.E_acc)) * self.p_cryo
        self.p_wp = self.pdyn + self.pstat

        print(self.p_cryo, self.l_active, self.l_cavity)
        return self.p_in, self.p_wp


def E_acc_Pin(cavity):
    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m

    cavity.set_Eacc(Eacc=E_acc)

    p_in, p_wp = cavity.get_power()

    # print("pin", p_in[np.where((E_acc >= 4.97) & (E_acc <= 5.01))])
    # print("pin", p_in[np.where((E_acc >= 7.05) & (E_acc <= 7.15))])
    # print("pin", p_in[np.where((E_acc >= 9.97) & (E_acc <= 10.01))])
    print("pin", p_in[np.where((E_acc >= 24.6* 1e6) & (E_acc <= 24.73* 1e6))])
    # print("pwp", p_wp[np.where((E_acc >= 4.97) & (E_acc <= 5.01))])
    # print("pwp", p_wp[np.where((E_acc >= 7.05) & (E_acc <= 7.15))])
    # print("pwp", p_wp[np.where((E_acc >= 9.97) & (E_acc <= 10.01))])
    print("pwp", p_wp[np.where((E_acc >= 24.6* 1e6) & (E_acc <= 24.73* 1e6))])
    # print("ncav", c3794.n_cav[np.where((E_acc >= 4.97) & (E_acc <= 5.01))])
    # print("ncav", c3794.n_cav[np.where((E_acc >= 7.05) & (E_acc <= 7.15))])
    # print("ncav", c3794.n_cav[np.where((E_acc >= 9.97) & (E_acc <= 10.01))])
    print("ncav", cavity.n_cav[np.where((E_acc >= 24.6* 1e6) & (E_acc <= 24.73* 1e6))])
    print("Q0", cavity.Q0[np.where((E_acc >= 24.6* 1e6) & (E_acc <= 24.73* 1e6))])
    print("dyn", (cavity.pdyn/cavity.n_cav)[np.where((E_acc >= 24.6* 1e6) & (E_acc <= 24.73* 1e6))])
    print("stat", (cavity.pstat/cavity.n_cav)[np.where((E_acc >= 24.6* 1e6) & (E_acc <= 24.73* 1e6))])

    fig, ax = plt.subplots()
    ax_right = ax.twinx()
    ax_right._get_lines.prop_cycler = ax._get_lines.prop_cycler
    ax_right2 = ax.twinx()
    ax_right2._get_lines.prop_cycler = ax._get_lines.prop_cycler
    ax_right2.spines["right"].set_position(("axes", 1.2))

    p1, = ax.plot(cavity.E_acc * 1e-6, p_wp/cavity.n_cav * 1e-3, lw=2, c='k', label="$P_\mathrm{dynamic} + P_\mathrm{static}$")
    ax.plot(cavity.E_acc * 1e-6, cavity.pdyn/cavity.n_cav * 1e-3, lw=2, label="$P_\mathrm{dynamic}$")
    ax.plot(cavity.E_acc * 1e-6, cavity.pstat/cavity.n_cav * 1e-3, lw=2, label="$P_\mathrm{static}$")
    p2, = ax_right.plot(cavity.E_acc * 1e-6, cavity.n_cav / 2, lw=2)
    p3, = ax_right2.plot(cavity.E_acc * 1e-6, p_in * 1e-3, lw=2)

    ax.set_xlabel("$E_\mathbf{acc}$ [MV/m]")
    ax.set_ylabel("$P_\mathbf{wp}$ [MW]")
    ax_right.set_ylabel("$N_\mathbf{cav}/beam$")
    ax_right2.set_ylabel("$P_\mathbf{in}$ [kW]")
    # ax.axvline(5, ls='--', c='k')
    # ax.axvline(7.13, ls='--', c='k')
    # ax.axvline(10, ls='--', c='k')
    # ax.axvline(15, ls='--', c='k')
    # ax_right2.axhline(500, ls='--', c='k')
    # ax_right2.axhline(1000, ls='--', c='k')
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
    # ax_right.yaxis.set_major_locator(ticker.MultipleLocator(100))
    # ax_right2.yaxis.set_major_locator(ticker.MultipleLocator(200))

    ax.yaxis.label.set_color(p1.get_color())
    ax_right.yaxis.label.set_color(p2.get_color())
    ax_right2.yaxis.label.set_color(p3.get_color())

    ax.set_xlim(0, 30)
    # ax.set_ylim(0, 50)
    ax_right.set_ylim(100, 400)
    ax_right2.set_ylim(0, 700)

    tkw = dict(size=4, width=1.5)
    ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
    ax_right.tick_params(axis='y', colors=p2.get_color(), **tkw)
    ax_right2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    ax.tick_params(axis='x', **tkw)

    ax.minorticks_on()
    mplcursors.cursor(ax)
    mplcursors.cursor(ax_right)
    mplcursors.cursor(ax_right2)
    ax.legend()
    # ax.grid(True, which='both', axis='both')
    plt.tight_layout()
    plt.show()


def QL_Pin(labels, geometry, RF, QOI, Machine, p_data=None):
    # check if entries are of same length

    it = iter(geometry)
    the_len = len(next(it))
    if not all(len(l) == the_len for l in it):
        raise ValueError('not all lists have same length!')

    it = iter(RF)
    the_len = len(next(it))
    if not all(len(l) == the_len for l in it):
        raise ValueError('not all lists have same length!')

    it = iter(QOI)
    the_len = len(next(it))
    if not all(len(l) == the_len for l in it):
        raise ValueError('not all lists have same length!')

    it = iter(Machine)
    the_len = len(next(it))
    if not all(len(l) == the_len for l in it):
        raise ValueError('not all lists have same length!')

    n_cells, l_cells, G, b = [np.array(x) for x in geometry]
    E_acc, Vrf = [np.array(x) for x in RF]

    fig, ax = plt.subplots()

    # QOI
    f0, R_Q = [np.array(x) for x in QOI]

    # Machine
    I0, rho, E0 = [np.array(x) for x in Machine]

    l_active = 2 * n_cells * l_cells
    l_cavity = l_active + 8 * l_cells

    # CALCULATED
    v_cav = E_acc * l_active

    U_loss = 88.46 * E0 ** 4 / rho * 1e-6  # GeV # energy lost per turn per beam
    v_loss = 2 * U_loss * 1e9  # V # v loss for two beams

    print(v_loss, Vrf, v_loss / Vrf)
    phi = np.arccos(v_loss / Vrf)
    delta_f = -R_Q * f0 * I0 * np.sin(phi) / (2 * v_cav)  # optimal df
    QL_0_x = v_cav / (R_Q * I0 * np.cos(phi))  # optimal Q loaded

    QL_0 = np.linspace(1e3, 1e8, 1000000)

    xy_list = [(0.3, 0.6), (0.35, 0.25), (0.6, 0.6), (0.77, 0.25)]
    for i in range(len(E_acc)):
        f1_2 = f0[i] / (2 * QL_0)  # 380.6
        ic(R_Q[i], v_cav[i], Vrf[i])
        pin = v_cav[i] ** 2 / (4 * R_Q[i] * QL_0) * (
                    (1 + ((R_Q[i] * QL_0 * I0[i]) / v_cav[i]) * np.cos(phi[i])) ** 2 + (
                    (delta_f[i] / f1_2) + ((R_Q[i] * QL_0 * I0[i]) / v_cav[i]) * np.sin(phi[i])) ** 2)

        p_cryo = 8 / (np.sqrt(f0[i] / 500e6))

        # material/ wall power
        e_acc = np.linspace(0.5, 25, 1000) * 1e6  # MV/m
        Rs_NbCu_4_5k_400Mhz = 39.5 * np.exp(0.014 * (E_acc[i] * 1e-6 * b[i])) + 27
        eta = 1 / 219  # c
        Q0 = G[i] * 1e9 / Rs_NbCu_4_5k_400Mhz  # c
        print(E_acc[i])
        p_wp = (1 / eta) * Vrf[i] * (E_acc[i] * l_active[i] / (R_Q[i] * Q0)) + (1 / eta) * (
                l_cavity[i] * Vrf[i] / (l_active[i] * E_acc[i])) * p_cryo

        if "*" in label[i]:
            l = ax.plot(QL_0, pin * 1e-3, label=f"${round(E_acc[i] * 1e-6, 2)}" + " ~[\mathrm{MV/m}]$", lw=4, ls='--')
        else:
            l = ax.plot(QL_0, pin * 1e-3, label=f"${round(E_acc[i] * 1e-6, 2)}" + " ~[\mathrm{MV/m}]$", lw=4)

        # add annotations
        txt = f"{labels[i]}, {n_cells[i]}-Cell 400 MHz \n {int(np.ceil(Vrf[i] / (2 * E_acc[i] * l_active[i]) / 2) * 2)} cav/beam" \
              + "\n V$_\mathbf{RF}\mathbf{/beam}$ =" + f"{round(Vrf[i] / 2 * 1e-9, 2)} GV " \
              + "\n V$_\mathbf{cav}$ =" + f"{round(v_cav[i] * 1e-6, 1)} MV \n " \
                                          "P$_{\mathrm{in}}$ = " + f"{round(min(pin) * 1e-3, 1)} kW \n" \
                                                                   "Q$_{\mathrm{L, 0}}^*$ = " + "{:.2e}".format(
            QL_0_x[i])

        annotext = ax.annotate(txt, xy=xy_list[i], xycoords='figure fraction', size=12, rotation=0, c=l[0].get_color(),
                               weight='bold')

        dt = DraggableText(annotext)
        dt.connect()

    if p_data:
        # plot QL with penetration
        ax_2 = ax.twinx()
        data = fr.excel_reader(p_data)
        data_ = data[list(data.keys())[0]]
        ax_2.plot(data_["QL"], data_["penetration"], lw=4)

    # plot decorations
    ax.set_xlabel("$Q_{L,0}$")
    ax.set_ylabel("$P_\mathrm{in} ~[\mathrm{kW}]$")
    ax.set_xscale('log')
    # ax.set_xlim(5e3, 5e7)
    # ax.set_ylim(100, 2000)
    ax.legend(loc='lower left', title="$E_\mathrm{acc}$")  #
    ax.minorticks_on()
    # ax.grid(which='both')
    fig.show()


def plot_beampipe_decay(freq: list, fc_list):
    wl = 4 * 187e-3
    c = 299792458
    z = np.linspace(0, 3 * wl, 100)  # m

    fig, ax = plt.subplots()
    for f, fc in zip(freq, fc_list):
        beta = -(2 * np.pi * f) / c * np.lib.scimath.sqrt(1 - (fc / f) ** 2)

        prop = np.exp(beta.imag * z)
        print(beta.imag)
        # print(prop)

        # convert prop to dB
        pdB = 20 * np.log10(prop)
        ax.plot(z, pdB, lw=4)

    ax.axhline(-120, ls='--', c='k')
    ax.axvline(1.12, ls='--', c='k')
    ax.set_xlabel('$z ~[\mathrm{m}]$')
    ax.set_ylabel(r'$\alpha ~[\mathrm{dB}]$')

    # xticks format
    ax.set_xlim(0, 3 * wl)
    ax.set_ylim(-250, 0)

    x = [wl * i / 2 for i in range(7)]
    labels = [fr"{round(i / wl, 2)}" + r"$\lambda$" for i in np.linspace(0, 3 * wl, 7)]

    ax.set_xticks(x)

    ax.set_xticklabels(labels)
    # ax.set_yticks(list(plt.yticks()[0]) + [-120])

    fig.show()


def plot_surface_resistance():
    Eacc = np.linspace(0.5, 30, 100)  # MV/m
    b = 4.88  # mT/MV/m

    Rs_NbCu_2k_400Mhz = 0.57 * (Eacc * b) + 28.4  # nOhm
    Rs_NbCu_4_5k_400Mhz = 39.5 * np.exp(0.014 * (Eacc * b)) + 27  # nOhm
    Rs_bulkNb_2k_400Mhz = (2.33 / 1000) * (Eacc * b) ** 2 + 26.24  # nOhm
    Rs_bulkNb_4_5k_400Mhz = 0.0123 * (Eacc * b) ** 2 + 62.53  # nOhm

    Rs_NbCu_2k_800Mhz = 1.45 * (Eacc * b) + 92  # nOhm
    Rs_NbCu_4_5k_800Mhz = 50 * np.exp(0.033 * (Eacc * b)) + 154  # nOhm
    Rs_bulkNb_2k_800Mhz = (16.4 + Eacc * b * 0.092) * (800 / 704) ** 2  # nOhm
    Rs_bulkNb_4_5k_800Mhz = 4 * (62.7 + (Eacc * b) ** 2 * 0.012)  # nOhm

    sr_list = [Rs_NbCu_2k_400Mhz, Rs_NbCu_4_5k_400Mhz, Rs_bulkNb_2k_400Mhz, Rs_bulkNb_4_5k_400Mhz,
               Rs_NbCu_2k_800Mhz, Rs_NbCu_4_5k_800Mhz, Rs_bulkNb_2k_800Mhz, Rs_bulkNb_4_5k_800Mhz]

    sr_labels = ['Rs_NbCu_2k_400Mhz', 'Rs_NbCu_4_5k_400Mhz', 'Rs_bulkNb_2k_400Mhz', 'Rs_bulkNb_4_5k_400Mhz',
                 'Rs_NbCu_2k_800Mhz', 'Rs_NbCu_4_5k_800Mhz', 'Rs_bulkNb_2k_800Mhz', 'Rs_bulkNb_4_5k_800Mhz']

    for i, sr in enumerate(sr_list):
        if '400Mhz' in sr_labels[i]:
            lw = 2
            ls = '--'

        if '800Mhz' in sr_labels[i]:
            lw = 2
            ls = '-'
        if '2k' in sr_labels[i]:
            c = 'tab:blue'

        if '4_5k' in sr_labels[i]:
            c = 'tab:orange'

        plt.plot(Eacc, sr, label=sr_labels[i], c=c, lw=lw, ls=ls)

    plt.legend()
    plt.xlabel("$E$ [MV/m]")
    plt.ylabel(r"$R_s [\mathrm{n\Omega}]$")
    plt.tight_layout()
    plt.yscale('log')
    plt.show()

plot_surface_resistance()
# we are working with two beams
# l_cell = 0.187  #
# n_cells = 2  #4 #
# freq = 400.79e6  #
# Vrf = 2*0.44e9  #2*0.75e9  #  RF voltage per beam multiplied by two beams
# PSR = 50e6  # synchrotron power/beam
# G = 198.42  # 273.2 # 0.00948*Q_factor*(f/1300)**0.5
# R_Q = 152.8 # 411  #  c Ohm linac definition
# eta = 1 / 219  # c

folder_name = "D:\Dropbox\CEMCodesHub\Cavity800\SimulationData\SLANS\CavityC3795"
# folder_name = "D:\Dropbox\CEMCodesHub\Cavity800\SimulationData\SLANS\CavityFCC_UROS5"
qois = 'qois.json'
geom_params = 'geometric_parameters.json'

d_qois = {}
with open(fr"{folder_name}\{qois}") as json_file:
    d_qois.update(json.load(json_file))

d_geom_params = {}
with open(fr"{folder_name}\{geom_params}") as json_file:
    d_geom_params.update(json.load(json_file))

print(d_qois)
print(d_geom_params)
n_cells = 5
l_cell_mid = d_geom_params['IC'][5] * 1e-3
# l_cell_end = d_geom_params['OC'][5] * 1e-3
freq = d_qois['freq [MHz]'] * 1e6
R_Q = d_qois['R/Q [Ohm]']
GR_Q = d_qois['GR/Q [Ohm^2]']
G = GR_Q/R_Q
print(G)
Q = d_qois['Q []']
Epk_Eacc = d_qois['Epk/Eacc []']
Bpk_Eacc = d_qois['Bpk/Eacc [mT/MV/m]']

vrf = 4.4e9  # per beam
Vrf = 2 * vrf
eta = 1 / 219  # c

c3795 = Cavity(n_cells, l_cell_mid, freq, Vrf, R_Q, Q, G, Epk_Eacc, Bpk_Eacc)

cavity = c3795
# E_acc_Pin(cavity)

label = ["$\mathbf{Z^*}$", 'Z', "$\mathbf{W^*}$", 'W']

# # geometry
# n_cells = [5, 5, 5]  # 4
# l_cell = [0.187, 0.187, 0.187]  # m
# G = [181.91, 215.47, 206.84]  # 273.2 # 0.00948*Q_factor*(f/1300)**0.5
# b = [5.32, 4.41, 4.76]
# geometry = [n_cells, l_cell, G, b]
#
# # QOI
# R_Q = [425.56, 519.26, 560.47]  # 411   # c Ohm linac definition
# f0 = [801.58e6, 801.58e6, 801.58e6]
# QOI = [f0, R_Q]
#
# # RF
# # Vrf = [2*0.1e9, 2*0.1e9, 2*0.75e9]  #   #2*0.75e9
# # Vrf = [2*0.12e9, 2*0.12e9, 2*1e9, 2*0.44e9]  #   #2*0.75e9 update
# Vrf = [9.19e9, 9.19e9, 9.19e9]  #   #ttbar
# # Eacc = [20e6, 20e6, 20e6]
# Eacc = [24.52e6, 24.52e6, 24.52e6]  # update
# RF = [Eacc, Vrf]
#
# # MACHINE
# # I0 = [1390e-3, 1390e-3, 147e-3, 147e-3]  # mA
# I0 = [0.01, 0.01, 0.01]  # mA parameter update
# # rho = [10.76e3, 10.76e3, 10.76e3, 10.76e3]  # bending radius
# rho = [9.937e3, 9.937e3, 9.937e3]  # bending radius
# E0 = [182.5, 182.5, 182.5]  # Beam energy GeV
# machine = [I0, rho, E0]
#
# QL_Pin(label, geometry, RF, QOI, machine)
# # QL_Pin([10e6], p_data=r"D:\Dropbox\penetration_QL.xlsx")
# # QL_Pin()
#
# freq = [801.58e6]
# # plot_beampipe_decay(freq, [765e6])
#
# plt.tight_layout()
# plt.show()
