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
    mpl.rcParams['legend.fontsize'] = 16
    mpl.rcParams['legend.title_fontsize'] = 14

    mpl.rcParams['figure.figsize'] = [12, 6]
    mpl.rcParams['figure.dpi'] = 100


plot_settings()


class Cavities:
    def __init__(self, cavities_list=None):
        if cavities_list is None:
            self.cavities_list = []

        self.cavities_list = cavities_list

        self.returned_results = None
        self.ls = ['solid', 'dashed', 'dashdot']

        self.E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
        self.set_cavities_field()

    def set_cavities_field(self):
        for cav in self.cavities_list:
            cav.set_Eacc(Eacc=self.E_acc)

    def set_cavities_slans(self, dir_list):
        for i, dirc in enumerate(dir_list):
            print(dirc)
            self.cavities_list[i].set_slans_qois(dirc)

    def set_cavities_abci(self, dir_list):
        for i, dirc in enumerate(dir_list):
            self.cavities_list[i].set_abci_qois(dirc, working_point=wp, bunch_length=sigma)

    def compare_power(self, E_acc=None):
        if E_acc is not None:
            self.E_acc = E_acc
            self.set_cavities_field()

        results = []
        for i, cav in enumerate(self.cavities_list):
            # E_acc_Pin(cavity, op_field[i], ls[i], fig, ax, ax_right, ax_right2)
            results.append(self.qois(cav, cav.op_field * 1e-6))

        self.returned_results = results

    def qois(self, cavity, op_field):
        ind = np.where((E_acc >= 0.99 * op_field * 1e6) & (E_acc <= 1.01 * op_field * 1e6))
        qois = {
            r"$N_\mathrm{cav}$/beam": np.average(cavity.n_cav[ind]),
            r"$P_\mathrm{stat}$/cav [kW]": np.average((cavity.pstat / cavity.n_cav)[ind]) * 1e-3,
            r"$P_\mathrm{dyn}$/cav [kW]": np.average((cavity.pdyn / cavity.n_cav)[ind]) * 1e-3,
            r"$P_\mathrm{wp}$/cav [kW]": np.average((cavity.p_wp / cavity.n_cav)[ind]) * 1e-3,
            r"$P_\mathrm{in}$ [kW]": np.average(cavity.p_in[ind]) * 1e-3,
            r"$Q_\mathrm{0} \mathrm{[10^7]}$": np.average(cavity.Q0[ind] * 1e-7)
        }

        return qois

    def qois_hom(self):
        results = []
        for cavity in self.cavities_list:
            results.append({
                r"$|k_\parallel| \mathrm{[V/pC]}$": cavity.k_loss,
                r"$|k_\perp| \mathrm{[V/pC/m]}$": cavity.k_kick,
                r"$P_\mathrm{HOM} \mathrm{[kW]}$": cavity.phom
            })

        return results

    def plot_power_comparison(self, fig=None, ax_list=None):
        if fig is not None:
            fig = fig
            ax1, ax2, ax3 = ax_list
        else:
            # create figure
            fig = plt.figure()
            gs = fig.add_gridspec(2, 2)
            ax1 = fig.add_subplot(gs[:, 0])
            ax2 = fig.add_subplot(gs[0, 1])
            ax3 = fig.add_subplot(gs[1, 1])
    # def E_acc_Pin(self, cavity, E_acc, op_field, ls='-', ):

        # ax_right2._get_lines.prop_cycler = ax._get_lines.prop_cycler
        # ax_right2.spines["right"].set_position(("axes", 1.2))
        for i, cavity in enumerate(self.cavities_list):
            ax1.plot(cavity.E_acc * 1e-6, cavity.pstat / cavity.n_cav * 1e-3,
                     ls=self.ls[i], lw=2, c='tab:orange',
                     label=r"$P_\mathrm{static}$" + fr"{cavity.name}")

            ax1.plot(cavity.E_acc * 1e-6, cavity.pdyn / cavity.n_cav * 1e-3,
                    ls=self.ls[i], lw=2, c='tab:blue', label=r"$P_\mathrm{dynamic}$" + fr"{cavity.name}")

            p1, = ax1.plot(cavity.E_acc * 1e-6, cavity.p_wp / cavity.n_cav * 1e-3,
                          ls=self.ls[i], lw=2, c='k', label=r"$P_\mathrm{wp}$" + fr"{cavity.name}")

            p2, = ax2.plot(cavity.E_acc * 1e-6, cavity.n_cav / 2, ls=self.ls[i], lw=2, c='tab:red',
                           label=fr"{cavity.name}")
            p3, = ax3.plot(cavity.E_acc * 1e-6, cavity.p_in * 1e-3, ls=self.ls[i], lw=2, c='tab:purple',
                           label=fr"{cavity.name}")

            ax1.set_xlabel(r"$E_\mathrm{acc}$ [MV/m]")
            ax1.set_ylabel(r"$P_\mathrm{wp}$ [kW]")
            ax2.set_xlabel(r"$E_\mathrm{acc}$ [MV/m]")
            ax2.set_ylabel(r"$N_\mathrm{cav}$/beam")
            ax3.set_xlabel(r"$E_\mathrm{acc}$ [MV/m]")
            ax3.set_ylabel(r"$P_\mathrm{in}$ [kW]")

            ax1.axvline(cavity.op_field*1e-6, ls=':', c='k')
            ax1.text(cavity.op_field*1e-6 - 1, 0.3, f"{cavity.op_field*1e-6} MV/m",
                     size=14, rotation=90, transform=ax1.get_xaxis_transform())
            ax2.axvline(cavity.op_field*1e-6, ls=':', c='k')
            ax2.text(cavity.op_field*1e-6 - 1, 0.5, f"{cavity.op_field*1e-6} MV/m",
                     size=14, rotation=90, transform=ax2.get_xaxis_transform())
            ax3.axvline(cavity.op_field*1e-6, ls=':', c='k')
            ax3.text(cavity.op_field*1e-6 - 1, 0.3, f"{cavity.op_field*1e-6} MV/m",
                     size=14, rotation=90,
                     transform=ax3.get_xaxis_transform())

            # ax.axvline(7.13, ls='--', c='k')
            # ax.axvline(10, ls='--', c='k')
            # ax.axvline(15, ls='--', c='k')
            # ax_right2.axhline(500, ls='--', c='k')
            # ax_right2.axhline(1000, ls='--', c='k')
            # ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
            # ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
            # ax_right.yaxis.set_major_locator(ticker.MultipleLocator(100))
            # ax_right2.yaxis.set_major_locator(ticker.MultipleLocator(200))

            # ax.yaxis.label.set_color(p1.get_color())
            # ax_right.yaxis.label.set_color(p2.get_color())
            # ax_right2.yaxis.label.set_color(p3.get_color())

            ax1.set_xlim(min(E_acc)*1e-6, max(E_acc)*1e-6)
            ax2.set_xlim(min(E_acc)*1e-6, max(E_acc)*1e-6)
            ax3.set_xlim(min(E_acc)*1e-6, max(E_acc)*1e-6)
            # # ax.set_ylim(0, 50)
            # ax_right.set_ylim(100, 400)f
            # ax_right2.set_ylim(0, 700)
            ax1.set_yscale('log')
            ax2.set_yscale('log')
            ax3.set_yscale('log')

            # tkw = dict(size=4, width=1.5)
            # ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
            # ax_right.tick_params(axis='y', colors=p2.get_color(), **tkw)
            # ax_right2.tick_params(axis='y', colors=p3.get_color(), **tkw)
            # ax.tick_params(axis='x', **tkw)

            ax1.minorticks_on()
            mplcursors.cursor(ax1)
            mplcursors.cursor(ax2)
            mplcursors.cursor(ax3)
            # ax.grid(True, which='both', axis='both')

        # dummy lines with NO entries, just to create the black style legend
        dummy_lines = []
        for b_idx, b in enumerate(self.cavities_list):
            dummy_lines.append(ax1.plot([], [], c="gray", ls=self.ls[b_idx])[0])

        lines = ax1.get_lines()
        legend1 = ax1.legend([lines[i] for i in range(3)],
                             [r"$P_\mathrm{stat}$", r"$P_\mathrm{dyn}$", r"$P_\mathrm{stat}+P_\mathrm{dyn}$"], loc=3)
        legend2 = ax1.legend([dummy_lines[i] for i in range(len(self.cavities_list))], [cavity.name for cavity in self.cavities_list],
                             loc=0)
        ax1.add_artist(legend1)

        # ax1.legend(ncol=len(cavities))
        ax2.legend(loc='upper left')
        ax3.legend(loc=3)

        label = [r"$\mathbf{Z^*}$", 'Z', r"$\mathbf{W^*}$", 'W']
        plt.tight_layout()
        plt.show()

    def plot_compare_bar(self):
        # plot barchart
        data = [list(d.values()) for d in self.returned_results]
        x = list(self.returned_results[0].keys())
        X = np.arange(len(x))

        fig, ax = plt.subplots()
        width = 1 / len(x)
        for i, cav in enumerate(self.cavities_list):
            print(cav.name)
            ax.bar(X + i * width, data[i], width=width, label=self.cavities_list[i].name)

        ax.set_xticks([r + width for r in range(len(x))], x)
        # label = ["C3794_H (2-Cell)", "C3795_H (5-Cell)"]
        ax.legend(loc="upper left")

        plt.show()

    def plot_compare_hom_bar(self):
        # plot barchart
        self.hom_results = self.qois_hom()
        data = [list(d.values()) for d in self.hom_results]
        x = list(self.hom_results[0].keys())
        X = np.arange(len(x))

        fig, ax = plt.subplots()
        width = 1 / (len(x)+10)
        for i, cav in enumerate(self.cavities_list):
            print(type(X), type(i), type(width), type(data[i]), data)
            ax.bar(X + i * width, data[i], width=width, label=self.cavities_list[i].name)

        ax.set_xticks([r+width for r in range(len(x))], x)
        # label = ["C3794_H (2-Cell)", "C3795_H (5-Cell)"]
        # label = ["C3795_ttbar (5-Cell)", "FCCUROS5_ttbar (5-Cell)", "TELSA_ttbar (5-Cell)"]
        # ax.legend(label, loc="upper left")
        ax.legend(loc="upper right")

        plt.show()

    def add_cavity(self, cav):
        self.cavities_list.append(cav)

    def remove_cavity(self, cav):
        self.cavities_list.remove(cav)


class Cavity:
    def __init__(self, n_cells, l_cell_mid, freq, vrf, R_Q, G, Epk_Eacc, Bpk_Eacc, inv_eta=219, name="Unnamed", op_field=1e6):
        # geometric parameters
        # input
        self.k_loss = None
        self.k_kick = None
        self.phom = None
        self.inv_eta = None
        self.name = name
        self.op_field = op_field
        self.p_wp = None
        self.pstat = None
        self.pdyn = None
        self.p_cryo = None
        self.p_in = None
        self.Q0 = None
        self.l_cavity = None
        self.n_cav = None
        self.E_acc = None
        self.n_cells = n_cells
        self.l_cell_mid = l_cell_mid  # m
        self.v_rf = vrf
        self.R_Q = R_Q
        self.op_freq = freq  # Hz
        self.e = Epk_Eacc
        self.b = Bpk_Eacc
        self.G = G

        # calculated
        self.l_active = 2 * self.n_cells * self.l_cell_mid  # m
        # self.n_cav = self.v_rf / (self.E_acc * self.l_active)
        # self.l_cavity = self.n_cav * self.n_cells * self.l_active + (self.n_cav - 1) * 6 * self.l_active

        self.E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m

        # cryo parameters
        # self.eta = 1/219  # 400 MHz
        self.eta = 1 / inv_eta  # 800 MHz

        # power parameters
        # calculated
        self.p_sr = 50e6  # W

    def set_Eacc(self, Eacc):
        self.E_acc = Eacc
        self.n_cav = self.v_rf / (self.E_acc * self.l_active)
        self.l_cavity = self.l_active + 8 * self.l_cell_mid  # L_cavity = L_active + 2 λ0

        Rs_NbCu_2k_400Mhz = 0.57 * (Eacc * 1e-6 * self.b) + 28.4  # nOhm
        Rs_NbCu_4_5k_400Mhz = 39.5 * np.exp(0.014 * (Eacc * 1e-6 * self.b)) + 27  # nOhm
        Rs_bulkNb_2k_400Mhz = (2.33 / 1000) * (Eacc * 1e-6 * self.b) ** 2 + 26.24  # nOhm
        Rs_bulkNb_4_5k_400Mhz = 0.0123 * (Eacc * 1e-6 * self.b) ** 2 + 62.53  # nOhm

        Rs_NbCu_2k_800Mhz = 1.45 * (Eacc * 1e-6 * self.b) + 92  # nOhm
        Rs_NbCu_4_5k_800Mhz = 50 * np.exp(0.033 * (Eacc * 1e-6 * self.b)) + 154  # nOhm
        Rs_bulkNb_2k_800Mhz = (16.4 + Eacc * 1e-6 * self.b * 0.092) * (800 / 704) ** 2  # nOhm
        Rs_bulkNb_4_5k_800Mhz = 4 * (62.7 + (Eacc * 1e-6 * self.b) ** 2 * 0.012)  # nOhm

        if np.isclose(self.op_freq, 801.58e6):
            self.Q0 = self.G * 1e9 / Rs_bulkNb_2k_800Mhz  # c
            # ic("800 MHz")
        elif np.isclose(self.op_freq, 400.79e6):
            self.Q0 = self.G * 1e9 / Rs_NbCu_4_5k_400Mhz  # c
            # ic("400 MHz")

        self.get_power()

    def set_op_field(self, op_field):
        self.op_field = op_field

    def get_op_field(self):
        return self.op_field

    def set_inv_eta(self, inv_eta):
        self.inv_eta = inv_eta
        self.eta = 1/inv_eta

    def get_inv_eta(self):
        return self.inv_eta

    def set_vrf(self, vrf):
        self.v_rf = vrf

    def get_vrf(self):
        return self.v_rf

    def get_power(self):
        self.p_in = 2 * self.p_sr / self.n_cav  # maximum synchrotron radiation per beam multipled by two for two beams
        self.p_cryo = 8 / (np.sqrt(self.op_freq / 500e6))  # W/m

        self.pdyn = (1 / self.eta) * self.v_rf * (self.E_acc * self.l_active / (self.R_Q * self.Q0))
        self.pstat = (1 / self.eta) * (self.l_cavity * self.v_rf / (self.l_active * self.E_acc)) * self.p_cryo
        self.p_wp = self.pdyn + self.pstat

        return self.p_in, self.p_wp

    def set_slans_qois(self, folder_name):

        qois = 'qois.json'
        geom_params = 'geometric_parameters.json'

        d_qois = {}
        with open(fr"{folder_name}\{qois}") as json_file:
            d_qois.update(json.load(json_file))

        d_geom_params = {}
        with open(fr"{folder_name}\{geom_params}") as json_file:
            d_geom_params.update(json.load(json_file))

        # ic(d_qois)
        self.l_cell_mid = d_geom_params['IC'][5] * 1e-3
        self.op_freq = d_qois['freq [MHz]'] * 1e6
        self.R_Q = d_qois['R/Q [Ohm]']
        self.GR_Q = d_qois['GR/Q [Ohm^2]']
        self.G = self.GR_Q/self.R_Q
        # print(G)
        # self.Q = d_qois['Q []']
        self.e = d_qois['Epk/Eacc []']
        self.b = d_qois['Bpk/Eacc [mT/MV/m]']

        vrf = 4.4e9  # per beam
        Vrf = 2 * vrf

    def set_abci_qois(self, folder_name, working_point='', bunch_length=''):

        qois = 'qois.json'
        geom_params = 'geometric_parameters.json'

        d_qois = {}
        with open(fr"{folder_name}\{qois}") as json_file:
            d_qois.update(json.load(json_file))

        d_geom_params = {}
        with open(fr"{folder_name}\{geom_params}") as json_file:
            d_geom_params.update(json.load(json_file))

        # ic(d_geom_params)

        d_qois = d_qois[f'{working_point}_{bunch_length}']
        self.k_loss = d_qois['|k_loss| [V/pC]']
        self.k_kick = d_qois['|k_kick| [V/pC/m]']
        self.phom = d_qois['P_HOM [kW]']

        # ic(d_qois)

    def write_cavity_for_multipac(self):
        pass


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
    Eacc = np.linspace(0.5, 30, 50)  # MV/m
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
            ls = '--'

        if '800Mhz' in sr_labels[i]:
            ls = '-'

        if '2k' in sr_labels[i]:
            marker = ''

        if '4_5k' in sr_labels[i]:
            marker = '+'

        if 'NbCu' in sr_labels[i]:
            c = 'tab:orange'

        if 'bulkNb' in sr_labels[i]:
            c = 'tab:blue'

        plt.xlim(min(Eacc), max(Eacc))
        plt.plot(Eacc, sr, label=sr_labels[i], c=c, ls=ls, marker=marker)

    plt.legend()
    plt.xlabel("$E$ [MV/m]")
    plt.ylabel(r"$R_s [\mathrm{n\Omega}]$")
    plt.tight_layout()
    plt.yscale('log')
    plt.show()


if __name__ == '__main__':

    c3794_H = Cavity(2, l_cell_mid=187e-3, freq=400.79e6, vrf=2.1e9, R_Q=152.8,G=198.42,
                     Epk_Eacc=2.05, Bpk_Eacc=6.39, inv_eta=219, name="C3794_H", op_field=11.87e6)

    c3795_H = Cavity(5, l_cell_mid=93.5e-3, freq=801.58e6, vrf=2.1e9, R_Q=448.12, G=261.63,
                     Epk_Eacc=2.43, Bpk_Eacc=4.88, inv_eta=745, name="C3795_H", op_field=24.72e6)

    c3795_tt = Cavity(5, l_cell_mid=93.5e-3, freq=801.58e6, vrf=8.8e9, R_Q=448.12, G=261.63,
                      Epk_Eacc=2.43, Bpk_Eacc=4.88, inv_eta=745, name="C3795_ttbar", op_field=24.72e6)

    cFCCUROS5 = Cavity(5, l_cell_mid=93.5e-3, freq=801.58e6, vrf=8.8e9, R_Q=521.06, G=272.93,
                       Epk_Eacc=2.05, Bpk_Eacc=4.33, inv_eta=745, name="FCCUROS5_ttbar", op_field=24.72e6)

    cTESLA = Cavity(5, l_cell_mid=93.5e-3, freq=801.58e6, vrf=8.8e9, R_Q=558.684, G=271.72,
                    Epk_Eacc=2.14, Bpk_Eacc=4.54, inv_eta=745, name="TESLA_ttbar", op_field=24.72e6)

    parent_dir_slans = r"D:\Dropbox\CEMCodesHub\Cavity800\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CEMCodesHub\Cavity800\SimulationData\ABCI"
    wp = 'H'  # working point
    sigma = 'SR_2.5mm'

    # slans_dirs = [fr"{parent_dir_slans}\CavityC3795", fr"{parent_dir_slans}\CavityFCC_UROS5", fr"{parent_dir_slans}\CavityTESLA_800MHZ"]
    # abci_dirs = [fr"{parent_dir_abci}\CavityC3795", fr"{parent_dir_abci}\CavityFCC_UROS5", fr"{parent_dir_abci}\CavityTESLA_800MHZ"]
    slans_dirs = [fr"{parent_dir_slans}\Cavity3794", fr"{parent_dir_slans}\CavityC3795"]
    abci_dirs = [fr"{parent_dir_abci}\Cavity3794", fr"{parent_dir_abci}\CavityC3795"]
    cavities = Cavities([c3794_H, c3795_H])
    # cavities = Cavities([c3795_tt, cFCCUROS5, cTESLA])
    cavities.set_cavities_slans(slans_dirs)
    cavities.set_cavities_abci(abci_dirs)

    # op_field = [11.87, 24.72, 24.72]  # in MV/m

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_hom_bar()
