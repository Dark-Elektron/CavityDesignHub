import json
import os.path
import scipy.io as spio
import scipy.interpolate as sci
import mplcursors
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from icecream import ic
import matplotlib as mpl
from scipy.optimize import fsolve
from utils.file_reader import FileReader
from analysis_modules.plot_module.matplotlib_annotation_objects import DraggableText
from utils.shared_functions import ellipse_tangent, writeCavityForMultipac_flat_top, writeCavityForMultipac

fr = FileReader()

m_e = 9.1093879e-31  # mass of electron
m_mu = 1.883531627e-28  # mass of muon
m_p = 1.67262192e-27  # mass of proton
q0 = 1.60217663e-19  # charge of electron, muon and proton
mu0 = 4 * np.pi * 1e-7  # permeability of free space
eps0 = 8.85418782e-12  # permittivity of free space
c0 = 2.99792458e8
eta0 = 376.7303134111465

PARTICLE = {
    "electron": {
        "m [kg]": m_e,
        "q [C]": q0,
        "rc [m]": 1 / (4 * np.pi * eps0) * q0 ** 2 / (m_e * c0 ** 2)  # classical radius rc = 1/(4pi*e0) * Q^2/(m c^2)
    },
    "muon": {
        "m [kg]": m_mu,
        "q [C]": q0,
        "rc [m]": 1 / (4 * np.pi * eps0) * q0 ** 2 / (m_mu * c0 ** 2)},
    "proton": {
        "m [kg]": m_p,
        "q [C]": q0,
        "rc [m]": 1 / (4 * np.pi * eps0) * q0 ** 2 / (m_p * c0 ** 2)
    }
}

WP = {
    'Z_2023': {'E [GeV]': 45.6,
               'I0 [mA]': 1280
               },
    'W_2023': {'E [GeV]': 80,
               'I0 [mA]': 135
               },
    'H_2023': {'E [GeV]': 120,
               'I0 [mA]': 54.3
               },
    'ttbar_2023': {'E [GeV]': 182.5,
                   'I0 [mA]': 10
                   },
    'MuCol RCS Stage 1': {"E [GeV]": 250.835,
                          "I0 [mA]": 20.38
                          }
}

MACHINE = {"FCC-ee": {'rho [m]': 9937},
           "MuCol": {'rho [m]': 581.8}}


def plot_settings():
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16

    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['axes.titlesize'] = 16
    mpl.rcParams['legend.fontsize'] = 16
    mpl.rcParams['legend.title_fontsize'] = 14

    mpl.rcParams['figure.figsize'] = [12, 6]
    mpl.rcParams['figure.dpi'] = 100

    # Set the desired colormap
    plt.rcParams['axes.prop_cycle'] = plt.cycler('color', plt.cm.Set2.colors)


plot_settings()


class Cavities:
    """
    Cavities object is an object containing several Cavity objects.
    """

    def __init__(self, cavities_list=None, save_folder='None', machine="FCC-ee", particle='electron'):
        """Constructs all the necessary attributes of the Cavity object

        Parameters
        ----------
        cavities_list: list, array like
            List containing Cavity objects.

        save_folder: str
            Folder to save generated images, latex, excel, and text files.
        """

        if cavities_list is None:
            self.cavities_list = []

        self.p_qois = None
        self.fm_results = None
        self.hom_results = None
        self.save_folder = save_folder
        self.machine = machine
        self.particle = particle

        self.cavities_list = cavities_list

        self.returned_results = None
        self.ls = ['solid', 'dashed', 'dashdot', 'dotted',
                   'solid', 'dashed', 'dashdot', 'dotted',
                   'solid', 'dashed', 'dashdot', 'dotted']

        self.E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
        self.set_cavities_field(self.E_acc)

    def set_cavities_field(self, E_acc):
        """
        Sets cavities analysis field range.

        Returns
        -------

        """
        for cav in self.cavities_list:
            self.E_acc = E_acc
            cav.set_Eacc(Eacc=E_acc)

    # def set_cavities_slans(self, dir_list):
    #     """
    #     Set folders to read SLANS results from
    #
    #     Parameters
    #     ----------
    #     dir_list: list, array like
    #         List of SLANS directories for the corresponding cavities
    #
    #     Returns
    #     -------
    #
    #     """
    #     for i, dirc in enumerate(dir_list):
    #         print(dirc)
    #         self.cavities_list[i].set_eigenmode_qois(dirc)
    #
    # def set_cavities_abci(self, dir_list):
    #     """
    #     Set folders to read ABCI results from
    #
    #     Parameters
    #     ----------
    #     dir_list: list, array like
    #         List of ABCI directories for the corresponding cavities
    #
    #     Returns
    #     -------
    #
    #     """
    #     for i, dirc in enumerate(dir_list):
    #         self.cavities_list[i].set_abci_qois(dirc, working_point=wp, bunch_length=sigma)

    def compare_power(self, E_acc=None):
        if E_acc is not None:
            self.E_acc = E_acc
            self.set_cavities_field(E_acc)

        self.p_qois = []
        results = []
        for i, cav in enumerate(self.cavities_list):
            # E_acc_Pin(cavity, op_field[i], ls[i], fig, ax, ax_right, ax_right2)
            results.append(self.qois(cav, cav.op_field * 1e-6, E_acc))

        self.returned_results = results

    def qois(self, cavity, op_field, E_acc):
        """

        Parameters
        ----------
        cavity: object
            Cavity object
        op_field: float
            Cavity operating field

        Returns
        -------
        Dictionary containing quantities of interest (normed optional).
        """

        qois = {
            r"$N_\mathrm{cav}$": np.average(cavity.n_cav_op_field),
            r"$Q_0 ~\mathrm{[10^{10}}]}$": cavity.Q0_desired * 1e-10,
            # r"Rs [Ohm]$": np.average(cavity.Rs[ind]),
            r"$P_\mathrm{stat}$/cav [W]": cavity.pstat,
            r"$P_\mathrm{dyn}$/cav [W]": cavity.pdyn,
            # r"P_stat [W]": np.average(cavity.pstat[ind]),
            # r"P_dyn [W]": np.average(cavity.pdyn[ind]),
            r"$P_\mathrm{wp}$/cav [kW]": cavity.p_wp * 1e-3,
            r"$P_\mathrm{in}$/cav [kW]": cavity.p_in * 1e-3,
            r"$P_\mathrm{HOM}$/cav [kW]": cavity.phom[0],
            # r"$P_\mathrm{HOM}$ [kW]": cavity.phom * np.average(cavity.n_cav[ind]),
            # r"P_tot_loss [kW]": np.average(cavity.pstat[ind]) * 1e-3
            #                     + np.average(cavity.pdyn[ind]) * 1e-3
            #                     + cavity.phom * np.average(cavity.n_cav[ind]),
            # r"$Q_\mathrm{0} \mathrm{[10^8]}$": np.average(cavity.Q0[ind] * 1e-8),
            # r"$Rs_\mathrm{0} \mathrm{[10^7]}$": np.average(cavity.Rs[ind])
        }
        self.p_qois.append(qois)

        # ind = np.where((E_acc >= 0.99 * op_field * 1e6) & (E_acc <= 1.01 * op_field * 1e6))
        # ic(op_field, np.average(cavity.n_cav[ind]))
        # qois_norm_units = {
        #     # r"$n_\mathrm{cav}$": np.average(cavity.n_cav[ind]),
        #     # r"$q_\mathrm{0}$": np.average(cavity.Q0[ind]),
        #     # r"$r_\mathrm{s}$": np.average(cavity.Rs[ind]),
        #     # r"$p_\mathrm{stat/cav}$": np.average(cavity.pstat[ind] / cavity.n_cav[ind]),
        #     # r"$p_\mathrm{dyn/cav}$": np.average(cavity.pdyn[ind] / cavity.n_cav[ind]),
        #     r"$p_\mathrm{in}/\mathrm{cav}$": np.average(cavity.p_in[ind]),
        #     r"$p_\mathrm{stat}$": np.average(cavity.pstat[ind]),
        #     r"$p_\mathrm{dyn}$": np.average(cavity.pdyn[ind]),
        #     r"$p_\mathrm{wp/cav}$": np.average(cavity.p_wp[ind] / cavity.n_cav[ind]),
        #     r"$p_\mathrm{wp}$": np.average(cavity.p_wp[ind]),
        #     r"$p_\mathrm{HOM}$": cavity.phom * np.average(cavity.n_cav[ind]),
        #     # r"$p_\mathrm{loss, tot}$": np.average(cavity.pstat[ind]) * 1e-3
        #     #                            + np.average(cavity.pdyn[ind]) * 1e-3
        #     #                            + np.average(cavity.phom * cavity.n_cav[ind]),
        #     # r"$Q_\mathrm{0} \mathrm{[10^8]}$": np.average(cavity.Q0[ind] * 1e-8),
        #     # r"$Rs_\mathrm{0} \mathrm{[10^7]}$": np.average(cavity.Rs[ind])
        # }

        ic(qois)
        return qois

    def qois_fm(self):
        """
        Retrieves the fundamental mode quantities of interest

        Returns
        -------
        Dictionary containing fundamental mode quantities of interest (normed optional).
        """
        results = []
        for cav in self.cavities_list:
            results.append({
                r"$E_\mathrm{pk}/E_\mathrm{acc} [\cdot]$": cav.e,
                r"$B_\mathrm{pk}/E_\mathrm{acc} \mathrm{[mT/MV/m]}$": cav.b,
                r"$k_\mathrm{cc}$": cav.k_cc,
                r"$R/Q \mathrm{[10^2\Omega]}$": cav.R_Q * 1e-2,
                r"$G \mathrm{[10^{2}\Omega]}$": cav.G * 1e-2,
                # r"$G\cdot R/Q \mathrm{[10^{5}\Omega^2]}$": cav.GR_Q * 1e-5
            })

        results_norm_units = []
        for cav in self.cavities_list:
            results_norm_units.append({
                r"$e_\mathrm{pk}/e_\mathrm{acc}$": cav.e,
                r"$b_\mathrm{pk}/e_\mathrm{acc}$": cav.b,
                r"$k_\mathrm{cc}$": cav.k_cc,
                r"$r/q$": cav.R_Q,
                r"$g$": cav.G,
                # r"$g\cdot r/q $": cav.GR_Q
            })
        ic(results)
        return results_norm_units

    def qois_hom(self):
        """
        Retrieves the higher-order modes quantities of interest


        Returns
        -------
        Dictionary containing higher-order modes quantities of interest (normed optional).
        """

        results = []
        for cavity in self.cavities_list:
            results.append({
                r"$|k_\parallel| \mathrm{[V/pC]}$": cavity.k_loss[0],
                r"$|k_\perp| \mathrm{[V/pC/m]}$": cavity.k_kick[0],
                r"$P_\mathrm{HOM}/cav \mathrm{[kW]}$": cavity.phom[0]
            })

        results_norm_units = []
        for cavity in self.cavities_list:
            results_norm_units.append({
                r"$k_\parallel$": cavity.k_loss[0],
                r"$k_\perp$": cavity.k_kick[0],
                r"$p_\mathrm{HOM}/cav$": cavity.phom[0]
            })
        ic(results)
        return results_norm_units

    def qois_all(self):
        """
        Retrieves the fundamental mode quantities of interest

        Returns
        -------
        Dictionary containing fundamental mode quantities of interest (normed optional).
        """
        results = []
        for cav in self.cavities_list:
            results.append({
                r"$E_\mathrm{pk}/E_\mathrm{acc} ~[\cdot]$": cav.e,
                r"$B_\mathrm{pk}/E_\mathrm{acc} ~\mathrm{[mT/MV/m]}$": cav.b,
                r"$k_\mathrm{cc}$": cav.k_cc,
                r"$R/Q~ \mathrm{[\Omega]}$": cav.R_Q,
                r"$G ~\mathrm{[\Omega]}$": cav.G,
                r"$G\cdot R/Q ~\mathrm{[10^4\Omega^2]}$": cav.GR_Q * 1e-4,
                r"$|k_\parallel| ~\mathrm{[V/pC]}$": cav.k_loss[0],
                r"$|k_\perp| ~\mathrm{[V/pC/m]}$": cav.k_kick[0],
                r"$P_\mathrm{HOM}/cav~ \mathrm{[kW]}$": cav.phom[0]
            })

        results_norm_units = []
        for cav in self.cavities_list:
            results_norm_units.append({
                r"$e_\mathrm{pk}/e_\mathrm{acc}$": cav.e,
                r"$b_\mathrm{pk}/e_\mathrm{acc}$": cav.b,
                r"$k_\mathrm{cc}$": cav.k_cc,
                r"$r/q$": cav.R_Q,
                r"$g$": cav.G,
                # r"$g\cdot r/q $": cav.GR_Q,
                # r"$|k_\mathrm{FM}|$": cav.k_fm,
                r"$|k_\parallel|$": cav.k_loss[0],
                r"$k_\perp$": cav.k_kick[0],
                r"$p_\mathrm{HOM}/cav$": cav.phom[0]
            })
        ic(results)
        return results

    def plot_power_comparison(self, fig=None, ax_list=None):
        """
        Can be called using ``cavities.plot_power_comparison()``

        .. math::

           W^{3 \\beta}_{\delta}

        Parameters
        ----------
        fig: matplotlib figure
        ax_list: list of matplotlib axes object

        Returns
        -------

        """
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
            # # per cavity
            # ax1.plot(cavity.E_acc * 1e-6, cavity.pstat / cavity.n_cav,
            #          ls=self.ls[i], lw=2, c='tab:orange',
            #          label=r"$P_\mathrm{stat$/cav}" + fr"{cavity.name}")
            #
            # ax1.plot(cavity.E_acc * 1e-6, cavity.pdyn / cavity.n_cav,
            #          ls=self.ls[i], lw=2, c='tab:blue', label=r"$P_\mathrm{dyn}/cav$" + fr"{cavity.name}")

            # total
            ax1.plot(cavity.E_acc * 1e-6, cavity.pstat,
                     ls=self.ls[i], lw=2, c='tab:orange',
                     label=r"$P_\mathrm{stat$}" + fr"{cavity.plot_label}")

            ax1.plot(cavity.E_acc * 1e-6, cavity.pdyn,
                     ls=self.ls[i], lw=2, c='tab:blue', label=r"$P_\mathrm{dyn}$" + fr"{cavity.plot_label}")

            # p1, = ax1.plot(cavity.E_acc * 1e-6, cavity.p_wp/cavity.n_cav,
            #                ls=self.ls[i], lw=2, c='k', label=r"$P_\mathrm{wp}$" + fr"{cavity.name}")

            p2, = ax2.plot(cavity.E_acc * 1e-6, cavity.n_cav, ls=self.ls[i], lw=2, c='tab:red',
                           label=fr"{cavity.plot_label}")

            p3, = ax3.plot(cavity.E_acc * 1e-6, cavity.p_in * 1e-3, ls=self.ls[i], lw=2, c='tab:purple',
                           label=fr"{cavity.plot_label}")

            ax1.set_xlabel(r"$E_\mathrm{acc}$ [MV/m]")
            ax1.set_ylabel(r"$P_\mathrm{stat, dyn}$ [W]")
            ax2.set_xlabel(r"$E_\mathrm{acc}$ [MV/m]")
            ax2.set_ylabel(r"$N_\mathrm{cav}$")
            ax3.set_xlabel(r"$E_\mathrm{acc}$ [MV/m]")
            ax3.set_ylabel(r"$P_\mathrm{in}/\mathrm{cav}$ [kW]")

            ax1.axvline(cavity.op_field * 1e-6, ls=':', c='k')
            ax1.text(cavity.op_field * 1e-6 - 1.5, 0.3, f"{round(cavity.op_field * 1e-6, 2)} MV/m",
                     size=14, rotation=90, transform=ax1.get_xaxis_transform())
            ax2.axvline(cavity.op_field * 1e-6, ls=':', c='k')
            ax2.text(cavity.op_field * 1e-6 - 1.5, 0.5, f"{round(cavity.op_field * 1e-6, 2)} MV/m",
                     size=14, rotation=90, transform=ax2.get_xaxis_transform())
            ax3.axvline(cavity.op_field * 1e-6, ls=':', c='k')
            ax3.text(cavity.op_field * 1e-6 - 1.5, 0.3, f"{round(cavity.op_field * 1e-6, 2)} MV/m",
                     size=14, rotation=90,
                     transform=ax3.get_xaxis_transform())

            ax1.set_xlim(min(cavity.E_acc) * 1e-6, max(cavity.E_acc) * 1e-6)
            ax2.set_xlim(min(cavity.E_acc) * 1e-6, max(cavity.E_acc) * 1e-6)
            ax3.set_xlim(min(cavity.E_acc) * 1e-6, max(cavity.E_acc) * 1e-6)
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
                             [r"$P_\mathrm{stat}$", r"$P_\mathrm{dyn}$"], loc=3)
        legend2 = ax1.legend([dummy_lines[i] for i in range(len(self.cavities_list))],
                             [cavity.name for cavity in self.cavities_list],
                             loc=0)
        ax1.add_artist(legend1)

        # ax1.legend(ncol=len(cavities))
        ax2.legend(loc='upper left')
        ax3.legend(loc=3)

        label = [r"$\mathbf{Z^*}$", 'Z', r"$\mathbf{W^*}$", 'W']
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_power_comparison.png")

        fig.show()

    def plot_compare_bar_(self):
        """
        Plots bar chart of power quantities of interest

        Returns
        -------

        """

        plt.rcParams["figure.figsize"] = (12, 4)
        # plot barchart
        data = np.array([list(d.values()) for d in self.returned_results])
        data_col_max = data.max(axis=0)

        x = list(self.returned_results[0].keys())
        X = np.arange(len(x))

        fig, ax = plt.subplots()
        width = min(0.15, 1 / (len(x) + 10))
        for i, cav in enumerate(self.cavities_list):
            ax.bar(X + i * width, data[i] / data_col_max, width=width, label=cav.plot_label, edgecolor='k')

        ax.set_xticks([r + width for r in range(len(x))], x)

        ax.axhline(1.05, c='k')
        ax.set_ylim(0.0, 1)
        ax.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                  mode="expand", borderaxespad=0, ncol=min(3, len(self.cavities_list)))
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_power_comparison_bar.png")

        fig.show()

    def plot_compare_bar(self, ncols=3):
        """
        Plots bar chart of power quantities of interest

        Returns
        -------

        """

        plt.rcParams["figure.figsize"] = (12, 4)
        # plot barchart
        ic(self.returned_results)
        df = pd.DataFrame.from_dict(self.returned_results)
        fig, axd = plt.subplot_mosaic([list(df.columns)], layout='constrained')

        labels = [cav.plot_label for cav in self.cavities_list]
        # Plot each column in a separate subplot
        for key, ax in axd.items():
            ax.bar(df.index, df[key], label=labels, color=plt.cm.get_cmap('Set2').colors[:len(df)], edgecolor='k',
                   width=1)
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_ylabel(key)
            h, l = ax.get_legend_handles_labels()

        # fig.set_tight_layout(True)
        if not ncols:
            ncols = min(4, len(self.cavities_list))
        fig.legend(h, l, loc='outside upper center', borderaxespad=0, ncol=ncols)

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_power_comparison_bar.png")

        fig.show()

    def plot_compare_hom_bar_(self):
        """
        Plot bar chart of higher-order mode's quantities of interest

        Returns
        -------

        """
        plt.rcParams["figure.figsize"] = (12, 4)
        # plot barchart
        self.hom_results = self.qois_hom()
        data = np.array([list(d.values()) for d in self.hom_results])
        data_col_max = data.max(axis=0)
        x = list(self.hom_results[0].keys())
        X = np.arange(len(x))

        fig, ax = plt.subplots()
        width = min(0.15, 1 / (len(x) + 10))
        for i, cav in enumerate(self.cavities_list):
            print(type(X), type(i), type(width), type(data[i]), data)
            ax.bar(X + i * width, data[i] / data_col_max, width=width, label=cav.plot_label, edgecolor='k')

        ax.set_xticks([r + width for r in range(len(x))], x)
        # label = ["C3794_H (2-Cell)", "C3795_H (5-Cell)"]
        # label = ["C3795_ttbar (5-Cell)", "FCCUROS5_ttbar (5-Cell)", "TELSA_ttbar (5-Cell)"]
        # ax.legend(label, loc="upper left")
        ax.axhline(1.05, c='k')
        ax.set_ylim(0.0, 1)
        ax.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                  mode="expand", borderaxespad=0, ncol=min(4, len(self.cavities_list)))
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_hom_bar.png")

        fig.show()

    def plot_compare_hom_bar(self):
        """
        Plot bar chart of higher-order mode's quantities of interest

        Returns
        -------

        """
        plt.rcParams["figure.figsize"] = (12, 4)
        # plot barchart
        self.hom_results = self.qois_hom()
        df = pd.DataFrame.from_dict(self.hom_results)
        fig, axd = plt.subplot_mosaic([list(df.columns)])
        # Plot each column in a separate subplot
        for key, ax in axd.items():
            ax.bar(df.index, df[key], color=plt.cm.get_cmap('Set2').colors[:len(df)], edgecolor='k', width=1)
            ax.set_xticklabels([])
            ax.set_ylabel(key)
        plt.tight_layout()
        fig.show()

    def plot_compare_fm_bar(self):
        """
        Plot bar chart of fundamental mode quantities of interest

        Returns
        -------

        """
        plt.rcParams["figure.figsize"] = (12, 4)
        # plot barchart
        self.fm_results = self.qois_fm()
        data = np.array([list(d.values()) for d in self.fm_results])
        data_col_max = data.max(axis=0)
        x = list(self.fm_results[0].keys())
        X = np.arange(len(x))

        fig, ax = plt.subplots()
        width = min(0.15, 1 / (len(x) + 10))
        width = 0.15
        for i, cav in enumerate(self.cavities_list):
            print(type(X), type(i), type(width), type(data[i]), data)
            ax.bar(X + i * width, data[i] / data_col_max, width=width, label=self.cavities_list[i].plot_label,
                   edgecolor='k')

        ax.set_xticks([r + width for r in range(len(x))], x)
        # label = ["C3794_H (2-Cell)", "C3795_H (5-Cell)"]
        # label = ["C3795_ttbar (5-Cell)", "FCCUROS5_ttbar (5-Cell)", "TELSA_ttbar (5-Cell)"]
        # ax.legend(label, loc="upper left")
        ax.axhline(1.05, c='k')
        ax.set_ylim(0.0, 1)
        ax.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                  mode="expand", borderaxespad=0, ncol=min(4, len(self.cavities_list)))
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_fm_bar.png")

        fig.show()

    def plot_compare_all_bar_(self, ncols=3):
        """
        Plot bar chart of fundamental mode quantities of interest

        Returns
        -------

        """
        plt.rcParams["figure.figsize"] = (10, 3)
        # plot barchart
        self.fm_results = self.qois_all()
        data = np.array([list(d.values()) for d in self.fm_results])
        data_col_max = data.max(axis=0)
        x = list(self.fm_results[0].keys())
        X = np.arange(len(x))

        fig, ax = plt.subplots()
        width = 0.15
        for i, cav in enumerate(self.cavities_list):
            # print(type(X), type(i), type(width), type(data[i]), data)
            ax.bar(X + i * width, data[i] / data_col_max, width=width, label=self.cavities_list[i].plot_label,
                   edgecolor='k')

        ax.set_xticklabels(ax.get_xticks(), rotation=0)
        ax.set_xticks([r + width for r in range(len(x))], x)
        # label = ["C3794_H (2-Cell)", "C3795_H (5-Cell)"]
        # label = ["C3795_ttbar (5-Cell)", "FCCUROS5_ttbar (5-Cell)", "TELSA_ttbar (5-Cell)"]
        # ax.legend(label, loc="upper left")
        # ax.axhline(1.05, c='k')
        ax.set_ylim(0.0, 1)
        ax.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                  mode="expand", borderaxespad=0, ncol=min(ncols, len(self.cavities_list)))
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_all_bar.png")

        fig.show()

    def plot_compare_all_bar(self, ncols=3):
        """
        Plot bar chart of fundamental mode quantities of interest

        Returns
        -------

        """
        plt.rcParams["figure.figsize"] = (15, 3)
        # plot barchart
        self.all_results = self.qois_all()

        df = pd.DataFrame.from_dict(self.all_results)
        fig, axd = plt.subplot_mosaic([list(df.columns)], layout='constrained')

        labels = [cav.plot_label for cav in self.cavities_list]
        # Plot each column in a separate subplot
        for key, ax in axd.items():
            ax.bar(df.index, df[key], label=labels, color=plt.cm.get_cmap('Set2').colors[:len(df)], edgecolor='k',
                   width=1)
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_ylabel(key)
            h, l = ax.get_legend_handles_labels()

        # fig.set_tight_layout(True)
        if not ncols:
            ncols = min(4, len(self.cavities_list))
        fig.legend(h, l, loc='outside upper center', borderaxespad=0, ncol=ncols)

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_all_bar.png")

        fig.show()

    def plot_cryomodule_comparison(self):
        """
        Plot cryomodule power comparison

        Returns
        -------

        """
        plt.rcParams["figure.figsize"] = (9, 3)
        fig, axs = plt.subplots(1, 2)
        n_cav_per_cryomodule = np.arange(1, 11)
        for cav in self.cavities_list:
            n_cryomodules_list = []
            cryomodules_len_list = []
            for ncpc in n_cav_per_cryomodule:
                cryomodule_len = cav.n_cells * (2 * cav.l_cell_mid) * ncpc + (ncpc + 1) * 8 * cav.l_cell_mid
                cryomodules_len_list.append(cryomodule_len)

                n_cryomodules = cav.n_cav_op_field / ncpc
                n_cryomodules_list.append(n_cryomodules)

            # ic(cryomodules_len_list)
            axs[0].plot(n_cav_per_cryomodule, cryomodules_len_list, marker='o', mec='k', label=f'{cav.plot_label}')
            axs[1].plot(n_cav_per_cryomodule, n_cryomodules_list, marker='o', mec='k', label=f'{cav.plot_label}')
            # ic(n_cav_per_cryomodule)
            # ic(cryomodules_len_list, n_cryomodules_list)
            # ic(n_cav_per_cryomodule)
            # ic()
        axs[0].set_xlabel("$N_\mathrm{cav}$/mod.")
        axs[0].set_ylabel("$L_\mathrm{cryo}$ [m]")
        axs[1].set_xlabel("$N_\mathrm{cav}$/mod.")
        axs[1].set_ylabel("$N_\mathrm{cryo}$")
        axs[0].legend()
        axs[1].legend()
        mplcursors.cursor(axs[0])
        mplcursors.cursor(axs[1])
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_cryo.png")

        plt.show()

    def plot_cavities_contour(self, opt='mid'):
        """Plot geometric contour of Cavity objects

        Parameters
        ----------
        opt: {"mid", "end", "all"}
            Either plot contour for only mid cells or end cells or the entire cavity
        n_cells: int
            Option used only when opt is set to "all"

        Returns
        -------

        """
        min_x, max_x, min_y, max_y = [], [], [], []

        if opt.lower() == 'mid' or opt.lower() == 'end':
            plt.rcParams["figure.figsize"] = (4, 5)
        else:
            plt.rcParams["figure.figsize"] = (10, 4)

        fig, ax = plt.subplots()

        for i, cav in enumerate(self.cavities_list):
            # write contour
            # self.write_contour(cav, opt)

            if opt.lower() == 'mid':
                mid_cell = np.array(cav.d_geom_params['IC']) * 1e-3
                end_cell_left = np.array(cav.d_geom_params['IC']) * 1e-3
                end_cell_right = np.array(cav.d_geom_params['IC']) * 1e-3
            elif opt.lower() == 'end':
                mid_cell = np.array(cav.d_geom_params['IC']) * 1e-3
                end_cell_left = np.array(cav.d_geom_params['OC']) * 1e-3
                end_cell_right = np.array(cav.d_geom_params['OC']) * 1e-3
            else:
                mid_cell = np.array(cav.d_geom_params['IC']) * 1e-3
                end_cell_left = np.array(cav.d_geom_params['OC']) * 1e-3
                end_cell_right = np.array(cav.d_geom_params['OC']) * 1e-3

            scale = cav.op_freq / c0
            if cav.cell_type == 'flat top':
                writeCavityForMultipac_flat_top(fr'{cav.slans_dir}\contour.txt', 1, mid_cell, end_cell_left,
                                                end_cell_right, beampipe='none', unit=1, scale=scale, plot=True)
            else:
                writeCavityForMultipac(fr'{cav.slans_dir}\contour.txt', 1, mid_cell, end_cell_left, end_cell_right,
                                       beampipe='none', unit=1, scale=scale)

            data = pd.read_csv(fr"{cav.slans_dir}\contour.txt", sep=r'\s+', header=None, skiprows=3)

            # ax.plot(data[1] * 1000, data[0] * 1000, lw=3., label=cav.plot_label)
            ax.plot(data[1], data[0], lw=3., label=cav.plot_label)
            ax.legend(loc='lower left')

            x_label = r"$z/\lambda$"
            y_label = r"$r/\lambda$"
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            min_x.append(min(data[1]))
            min_y.append(min(data[0]))
            max_x.append(max(data[1]))
            max_y.append(max(data[0]))

        if opt.lower() == 'mid' or opt.lower() == 'end':
            ax.set_xlim(0, max(max_x))
            ax.set_ylim(0, max(max_y))
        else:
            ax.set_xlim(min(min_x), max(max_x))
            ax.set_ylim(min(min_y), max(max_y))

        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_contour_{opt}.png")

        fig.show()

    def plot_axis_fields(self):
        """
        Plot axis fields of cavities

        Returns
        -------

        """
        for cav in self.cavities_list:
            # normalize fields
            e_axis = np.abs(cav.axis_field['1'])
            e_axis_norm = e_axis / e_axis.max()

            # shift to mid
            z = cav.axis_field['0']
            z_shift = z - z.max() / 2
            plt.plot(z_shift, e_axis_norm, label=cav.plot_label)

        plt.xlabel('$z$ [mm]')
        plt.ylabel('$|E_\mathrm{axis}|/|E_\mathrm{axis}|_\mathrm{max}$')
        plt.axhline(1.02, c='k')
        plt.ylim(-0.01, 1.5)
        plt.legend(loc='upper center', ncol=len(self.cavities_list))
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_axis_fields.png")

        plt.show()

    def plot_surface_fields(self):
        """
        Plot surface fields of cavities

        Returns
        -------

        """
        for cav in self.cavities_list:
            # normalize fields
            e_surf = np.abs(cav.surface_field['0'])
            e_surf_norm = e_surf / e_surf.max()

            plt.plot(e_surf_norm, label=cav.plot_label)

        plt.axhline(1.02, c='k')
        plt.ylim(-0.01, 1.5)
        plt.xlabel('$L_\mathrm{surf}$ [mm]')
        plt.ylabel('$|E_\mathrm{surf}|/|E_\mathrm{surf}|_\mathrm{max}$')
        plt.legend(loc='upper center', ncol=len(self.cavities_list))
        plt.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_surface_fields.png")

        plt.show()

    def plot_multipac_triplot(self, folders, kind='triplot'):
        """
        Plot Multipac triplot

        Parameters
        ----------
        folders: list, array like
            List of folder to read multipacting results from
        kind

        Notes
        -----
        This will be changed later so that the multipac results will be in the same location as the SLANS and ABCI
        results

        Returns
        -------

        """

        if kind == 'triplot':
            # create figure
            fig = plt.figure()
            gs = fig.add_gridspec(3, 1)
            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[1, 0])
            ax3 = fig.add_subplot(gs[2, 0])
            axs = [ax1, ax2, ax3]

        else:
            print("in here")
            fig, axs = plt.subplots(1, 1)
            axs = [axs]

        mpl.rcParams['figure.figsize'] = [6, 10]

        Eacc_list = [cav.op_field * 1e-6 for cav in self.cavities_list]
        Epk_Eacc_list = [cav.e for cav in self.cavities_list]
        labels = [cav.name for cav in self.cavities_list]
        for Eacc, Epk_Eacc, folder, label in zip(Eacc_list, Epk_Eacc_list, folders, labels):
            # load_output_data
            # files
            fnames = ["Ccounter.mat", "Acounter.mat", "Atcounter.mat", "Efcounter.mat", "param",
                      "geodata.n", "secy1", "counter_flevels.mat", "counter_initials.mat"]
            data = {}
            # files_folder = "D:\Dropbox\multipacting\MPGUI21"
            for f in fnames:
                if ".mat" in f:
                    data[f] = spio.loadmat(fr"{folder}\\{f}")
                else:
                    data[f] = pd.read_csv(fr"{folder}\\{f}", sep='\s+', header=None)

            A = data["Acounter.mat"]["A"]
            At = data["Atcounter.mat"]["At"]
            C = data["Ccounter.mat"]["C"]
            Ef = data["Efcounter.mat"]["Ef"]
            flevel = data["counter_flevels.mat"]["flevel"]
            initials = data["counter_initials.mat"]["initials"]
            secy1 = data["secy1"].to_numpy()
            Pow = flevel
            n = len(initials[:, 0]) / 2  # number of initials in the bright set
            N = int(data["param"].to_numpy()[4])  # number of impacts
            U = flevel
            Efl = flevel
            q = 1.6021773e-19
            Efq = Ef / q

            e1 = np.min(np.where(secy1[:, 1] >= 1))  # lower threshold
            e2 = np.max(np.where(secy1[:, 1] >= 1))  # upper threshold
            val, e3 = np.max(secy1[:, 1]), np.argmax(secy1[:, 1])  # maximum secondary yield

            cl = 0
            ok, ok1, ok2 = 1, 1, 1
            if ok > 0:
                if n == 0:
                    ic('Unable to plot the counters. No initial points.')
                    return

                if ok1 * ok2 == 0:
                    cl = ic('Counter functions or impact energy missing.')
                else:
                    # if ss > 0:
                    #     cl = ic(np.array(['Plotting the triplot (counter, enhanced ', 'counter and impact energy).']))

                    if kind == 'counter function' or kind == 'triplot':
                        # fig, axs = plt.subplots(3)
                        axs[0].plot(Efl / 1e6, C / n, lw=2, label=label)
                        axs[0].set_ylabel("$c_" + "{" + f"{N}" + "}/ c_0 $")
                        axs[0].set_xlabel(r'$E_\mathrm{pk}$ [MV/m]')
                        # axs[0].set_title(r'$\mathbf{MultiPac 2.1~~~~~Counter function~~~~}$')
                        axs[0].set_xlim(np.amin(Efl) / 1e6, np.amax(Efl) / 1e6)
                        axs[0].set_ylim(0, np.max([0.1, axs[0].get_ylim()[1]]))

                        # plot peak operating field
                        axs[0].axvline(Eacc * Epk_Eacc, c='k', ls='--', lw=2)
                        axs[0].text(np.round(Eacc * Epk_Eacc, 2) - 1.5, 0.1,
                                    f"{label[0]}: {np.round(Eacc * Epk_Eacc, 2)} MV/m",
                                    size=12, rotation=90,
                                    transform=axs[0].get_xaxis_transform())

                        axs[0].minorticks_on()
                        axs[0].set_ylim(axs[0].get_ylim()[0], 1.1 * axs[0].get_ylim()[-1])
                        axs[0].legend(loc='upper center', ncol=len(self.cavities_list))

                    if kind == 'final impact energy' or kind == 'triplot':
                        s = 0
                        if kind == 'final impact energy':
                            s = 1
                        axs[1 - s].semilogy(Efl / 1e6, Efq, lw=2, label=label)

                        # axs[1-s].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [secy1[e1, 0], secy1[e1, 0]], '-r')
                        e0 = sci.interp1d(secy1[0:e1 + 1, 1], secy1[0:e1 + 1, 0])(1)
                        axs[1 - s].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [e0, e0], '-r')
                        axs[1 - s].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [secy1[e2, 0], secy1[e2, 0]], '-r')
                        axs[1 - s].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [secy1[e3, 0], secy1[e3, 0]], '--r')

                        axs[1 - s].set_ylabel("$Ef_" + "{" + f"{N}" + "}$")
                        axs[1 - s].set_xlabel(r'$E_\mathrm{pk}$ [MV/m]')
                        # axs[1-s].set_title('$\mathbf{Final~Impact~Energy~in~eV}$')
                        axs[1 - s].set_xlim(np.min(Efl) / 1e6, np.max(Efl) / 1e6)
                        axs[1 - s].set_ylim(0, axs[1 - s].get_ylim()[1])

                        axs[1 - s].axvline(Eacc * Epk_Eacc, c='k', ls='--', lw=2)
                        axs[1 - s].text(np.round(Eacc * Epk_Eacc, 2) - 1.5, 0.1,
                                        f"{label[0]}: {np.round(Eacc * Epk_Eacc, 2)} MV/m",
                                        size=12, rotation=90,
                                        transform=axs[1 - s].get_xaxis_transform())

                        axs[1 - s].minorticks_on()
                        axs[1 - s].set_ylim(axs[1 - s].get_ylim()[0], 10 * axs[1 - s].get_ylim()[-1])
                        axs[1 - s].legend(loc='upper center', ncol=len(self.cavities_list))
                    if kind == 'enhanced counter function' or kind == 'triplot':
                        s = 0
                        if kind == 'enhanced counter function':
                            s = 2
                        axs[2 - s].semilogy(Efl / 1e6, (A + 1) / n, lw=2, label=label)
                        axs[2 - s].set_xlabel('$V$ [MV]')
                        axs[2 - s].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [1, 1], '-r')
                        axs[2 - s].set_xlim(np.min(Efl) / 1e6, np.max(Efl) / 1e6)
                        axs[2 - s].set_ylim(np.min((A + 1) / n), axs[2 - s].get_ylim()[1])
                        axs[2 - s].set_ylabel("$e_" + "{" + f"{N}" + "}" + "/ c_0$")
                        axs[2 - s].set_xlabel(r'$E_\mathrm{pk}$ [MV/m]')
                        # axs[2-s].set_title('$\mathbf{Enhanced~counter~function}$')

                        axs[2 - s].axvline(Eacc * Epk_Eacc, c='k', ls='--', lw=2)
                        axs[2 - s].text(np.round(Eacc * Epk_Eacc, 2) - 1, 0.1,
                                        f"{label[0]}: {np.round(Eacc * Epk_Eacc, 2)} MV/m",
                                        size=12, rotation=90,
                                        transform=axs[2 - s].get_xaxis_transform())

                        axs[2 - s].minorticks_on()
                        axs[2 - s].set_ylim(axs[2 - s].get_ylim()[0], 10 * axs[2 - s].get_ylim()[-1])
                        axs[2 - s].legend(loc='upper center', ncol=len(self.cavities_list))

        fig.tight_layout()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_{kind.replace(' ', '_')}.png")

        plt.show()

    def plot_dispersion(self):
        """
        Plot dispersion curve for the cavities

        Returns
        -------

        """
        fig, ax = plt.subplots()
        for cav in self.cavities_list:
            x = range(1, cav.n_cells + 1)
            ax.plot(x, cav.d_slans_all_results['FREQUENCY'][0:cav.n_cells], marker='o', mec='k',
                    label=f'{cav.plot_label} (kcc={round(cav.k_cc, 2)} %)')
            ax.set_xlabel('Mode Number')
            ax.set_ylabel('Frequency [MHz]')

        plt.legend()

        # save plots
        fname = [cav.name for cav in self.cavities_list]
        fname = '_'.join(fname)

        self.save_all_plots(f"{fname}_dispersion.png")

        plt.show()

    def write_contour(self, cav, opt='mid', n_cells=1):
        """
        Write geometric contour for cavities

        Parameters
        ----------
        cav: Cavity object
            Cavity object
        opt: str

        n_cells: int
            Number of cavity cells

        Returns
        -------

        """

        if opt.lower() == 'mid':
            A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, _ = np.array(cav.d_geom_params['IC']) * 1e-3
            A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, _ = np.array(cav.d_geom_params['IC']) * 1e-3
            A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, _ = np.array(cav.d_geom_params['IC']) * 1e-3
            n_cell = 1
            L_bp_l = 0.001
            L_bp_r = 0.001

            # calculate shift
            shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er) / 2

        elif opt.lower() == 'end':
            A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, _ = np.array(cav.d_geom_params['IC']) * 1e-3
            A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, _ = np.array(cav.d_geom_params['IC']) * 1e-3
            A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, _ = np.array(cav.d_geom_params['OC']) * 1e-3
            L_bp_l = 0.001
            L_bp_r = 1 * L_m

            n_cell = 1

            # calculate shift
            shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m) / 2
        else:
            A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, _ = np.array(cav.d_geom_params['IC']) * 1e-3
            A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, _ = np.array(cav.d_geom_params['OC']) * 1e-3
            try:
                A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, _ = np.array(cav.d_geom_params['OC_R']) * 1e-3
            except KeyError:
                A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, _ = np.array(cav.d_geom_params['OC']) * 1e-3

            L_bp_l = 4 * L_m
            L_bp_r = 4 * L_m

            n_cell = n_cells

            # calculate shift
            shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er) / 2

        step = 0.02  # step in boundary points in mm
        # shift = 0
        # shift = L_m  # for end cell

        # calculate angles outside loop
        # CALCULATE x1_el, y1_el, x2_el, y2_el
        data = ([0 + L_bp_l, Ri_el + b_el, L_el + L_bp_l, Req_el - B_el],
                [a_el, b_el, A_el, B_el])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])

        x1el, y1el, x2el, y2el = fsolve(ellipse_tangent, np.array(
            [a_el + L_bp_l, Ri_el + 0.85 * b_el, L_el - A_el + L_bp_l, Req_el - 0.85 * B_el]),
                                        args=data,
                                        xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

        # CALCULATE x1, y1, x2, y2
        data = ([0 + L_bp_l, Ri_m + b_m, L_m + L_bp_l, Req_m - B_m],
                [a_m, b_m, A_m, B_m])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        x1, y1, x2, y2 = fsolve(ellipse_tangent,
                                np.array([a_m + L_bp_l, Ri_m + 0.85 * b_m, L_m - A_m + L_bp_l, Req_m - 0.85 * B_m]),
                                args=data, xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

        # CALCULATE x1_er, y1_er, x2_er, y2_er
        data = ([0 + L_bp_r, Ri_er + b_er, L_er + L_bp_r, Req_er - B_er],
                [a_er, b_er, A_er, B_er])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        x1er, y1er, x2er, y2er = fsolve(ellipse_tangent, np.array(
            [a_er + L_bp_r, Ri_er + 0.85 * b_er, L_er - A_er + L_bp_r, Req_er - 0.85 * B_er]),
                                        args=data,
                                        xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

        with open(fr'{cav.slans_dir}\contour.txt', 'w') as fil:
            # SHIFT POINT TO START POINT
            start_point = [-shift, 0]
            fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   3.0000000e+00   0.0000000e+00\n")

            self.lineTo(start_point, [-shift, Ri_el], step)
            pt = [-shift, Ri_el]
            fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

            # ADD BEAM PIPE LENGTH
            self.lineTo(pt, [L_bp_l - shift, Ri_el], step)
            pt = [L_bp_l - shift, Ri_el]
            # print(pt)
            fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

            for n in range(1, n_cell + 1):
                if n == 1:
                    # DRAW ARC:
                    pts = self.arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el])
                    pt = [-shift + x1el, y1el]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW LINE CONNECTING ARCS
                    self.lineTo(pt, [-shift + x2el, y2el], step)
                    pt = [-shift + x2el, y2el]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                    pts = self.arcTo(L_el + L_bp_l - shift, Req_el - B_el, A_el, B_el, step, pt,
                                     [L_bp_l + L_el - shift, Req_el])
                    pt = [L_bp_l + L_el - shift, Req_el]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    if n_cell == 1:
                        # EQUATOR ARC TO NEXT POINT
                        # half of bounding box is required,
                        # start is the lower coordinate of the bounding box and end is the upper
                        # ic(pt, 1)
                        pts = self.arcTo(L_el + L_bp_l - shift, Req_er - B_er, A_er, B_er, step, [pt[0], Req_er - B_er],
                                         [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, Req_er])
                        pt = [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                        # ic(pt, 2)
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        # ic(pt, 3)

                        # STRAIGHT LINE TO NEXT POINT
                        self.lineTo(pt, [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step)
                        pt = [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        # ic(pt, 4, L_el + L_er - x1er + L_bp_l + L_bp_r - shift)

                        # ARC
                        # half of bounding box is required,
                        # start is the lower coordinate of the bounding box and end is the upper
                        # ic(shift)
                        pts = self.arcTo(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                         [L_bp_l + L_el + L_er - shift, y1er])
                        # ic(pt, 5, L_el + L_er + L_bp_l - shift)

                        pt = [L_bp_l + L_el + L_er - shift, Ri_er]
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")
                        # ic(pt, 6)

                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # calculate new shift
                        shift = shift - (L_el + L_er)
                        # ic(shift)
                    else:
                        # print("if else")
                        # EQUATOR ARC TO NEXT POINT
                        # half of bounding box is required,
                        # start is the lower coordinate of the bounding box and end is the upper
                        pts = self.arcTo(L_el + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, [pt[0], Req_m - B_m],
                                         [L_el + L_m - x2 + 2 * L_bp_l - shift, Req_m])
                        pt = [L_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # STRAIGHT LINE TO NEXT POINT
                        self.lineTo(pt, [L_el + L_m - x1 + 2 * L_bp_l - shift, y1], step)
                        pt = [L_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # ARC
                        # half of bounding box is required,
                        # start is the lower coordinate of the bounding box and end is the upper
                        pts = self.arcTo(L_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                         [L_bp_l + L_el + L_m - shift, y1])
                        pt = [L_bp_l + L_el + L_m - shift, Ri_m]
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # calculate new shift
                        shift = shift - (L_el + L_m)
                        # ic(shift)

                elif n > 1 and n != n_cell:
                    # print("elif")
                    # DRAW ARC:
                    pts = self.arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                    pt = [-shift + x1, y1]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW LINE CONNECTING ARCS
                    self.lineTo(pt, [-shift + x2, y2], step)
                    pt = [-shift + x2, y2]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                    pts = self.arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt,
                                     [L_bp_l + L_m - shift, Req_m])
                    pt = [L_bp_l + L_m - shift, Req_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = self.arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, [pt[0], Req_m - B_m],
                                     [L_m + L_m - x2 + 2 * L_bp_l - shift, Req_m])
                    pt = [L_m + L_m - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    self.lineTo(pt, [L_m + L_m - x1 + 2 * L_bp_l - shift, y1], step)
                    pt = [L_m + L_m - x1 + 2 * L_bp_l - shift, y1]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = self.arcTo(L_m + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                     [L_bp_l + L_m + L_m - shift, y1])
                    pt = [L_bp_l + L_m + L_m - shift, Ri_m]
                    # ic(pt)
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - 2 * L_m
                else:
                    # print("else")
                    # DRAW ARC:
                    pts = self.arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                    pt = [-shift + x1, y1]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW LINE CONNECTING ARCS
                    self.lineTo(pt, [-shift + x2, y2], step)
                    pt = [-shift + x2, y2]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                    pts = self.arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt,
                                     [L_bp_l + L_m - shift, Req_m])
                    pt = [L_bp_l + L_m - shift, Req_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = self.arcTo(L_m + L_bp_l - shift, Req_er - B_er, A_er, B_er, step, [pt[0], Req_er - B_er],
                                     [L_m + L_er - x2er + 2 * L_bp_l - shift, Req_er])
                    pt = [L_m + L_er - x2er + 2 * L_bp_l - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    self.lineTo(pt, [L_m + L_er - x1er + 2 * L_bp_l - shift, y1er], step)
                    pt = [L_m + L_er - x1er + 2 * L_bp_l - shift, y1er]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = self.arcTo(L_m + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                     [L_bp_l + L_m + L_er - shift, y1er])
                    pt = [L_bp_l + L_m + L_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

            # BEAM PIPE
            # reset shift
            # print("pt before", pt)
            shift = (L_bp_r + L_bp_l + (n_cell - 1) * 2 * L_m + L_el + L_er) / 2
            self.lineTo(pt, [L_bp_r + L_bp_l + 2 * (n_cell - 1) * L_m + L_el + L_er - shift, Ri_er], step)
            pt = [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, Ri_er]
            fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   3.0000000e+00   0.0000000e+00\n")
            # print("pt after", pt)

            # END PATH
            self.lineTo(pt, [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0],
                        step)  # to add beam pipe to right
            pt = [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0]
            # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
            # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
            fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

            # CLOSE PATH
            self.lineTo(pt, start_point, step)
            fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

        # plt.show()

    def write_contour_flat_top(self, cav, opt='mid', n_cells=1):
        plt.rcParams["figure.figsize"] = (12, 3)
        midJlab = np.array([64.453596, 54.579114, 19.1, 25.922107, 65, 83.553596, 163.975, 90,
                            20]) * 1e-3  # [A, B, a, b, Ri, L, Req, alpha, l]
        endJlab = np.array([64.453596, 54.579114, 19.1, 25.922107, 65, 83.553596, 163.975, 90, 11.187596]) * 1e-3
        endJlab_r = np.array([64.453596, 54.579114, 19.1, 25.922107, 65, 83.553596, 163.975, 90, 11.187596]) * 1e-3

        # cepc
        midCEPC = np.array(
            [94.4, 94.4, 20.1, 22.1, 78, 114.8873, 204.95, 90, 0.8254]) * 1e-3  # [A, B, a, b, Ri, L, Req, alpha, l]
        endCEPC = np.array([94.4, 94.4, 20.1, 22.1, 78, 114.8873, 204.95, 90, 0.8254]) * 1e-3
        endCEPC_r = np.array([94.4, 94.4, 20.1, 22.1, 78, 114.8873, 204.95, 90, 0.8254]) * 1e-3

        # # TESLA end cell 2
        A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, alpha, lft = midCEPC
        A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, alpha_el, lft_el = endCEPC
        A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, alpha_er, lft_er = midCEPC

        n_cell = 2
        step = 0.02  # step in boundary points in mm
        L_bp_l = 4 * L_m  # 0.0000  #
        L_bp_r = 4 * L_m  # 0.0000  #

        # calculate shift
        shift = (L_bp_r + L_bp_l + L_el + lft_el + (n_cell - 1) * 2 * L_m + (n_cell - 2) * lft + L_er + lft_er) / 2
        # shift = 0
        # shift = L_m  # for end cell

        # calculate angles outside loop
        # CALCULATE x1_el, y1_el, x2_el, y2_el
        data = ([0 + L_bp_l, Ri_el + b_el, L_el + L_bp_l, Req_el - B_el],
                [a_el, b_el, A_el, B_el])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])

        x1el, y1el, x2el, y2el = fsolve(f, np.array(
            [a_el + L_bp_l, Ri_el + 0.85 * b_el, L_el - A_el + L_bp_l, Req_el - 0.85 * B_el]),
                                        args=data, fprime=jac,
                                        xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

        # CALCULATE x1, y1, x2, y2
        data = ([0 + L_bp_l, Ri_m + b_m, L_m + L_bp_l, Req_m - B_m],
                [a_m, b_m, A_m, B_m])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        x1, y1, x2, y2 = fsolve(f, np.array([a_m + L_bp_l, Ri_m + 0.85 * b_m, L_m - A_m + L_bp_l, Req_m - 0.85 * B_m]),
                                args=data, fprime=jac, xtol=1.49012e-12)

        # CALCULATE x1_er, y1_er, x2_er, y2_er
        data = ([0 + L_bp_r, Ri_er + b_er, L_er + L_bp_r, Req_er - B_er],
                [a_er, b_er, A_er, B_er])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        x1er, y1er, x2er, y2er = fsolve(f, np.array(
            [a_er + L_bp_r, Ri_er + 0.85 * b_er, L_er - A_er + L_bp_r, Req_er - 0.85 * B_er]),
                                        args=data, fprime=jac,
                                        xtol=1.49012e-12)

        with open(r'D:\Dropbox\multipacting\MPGUI21\geodata.n', 'w') as fil:
            fil.write("   2.0000000e-03   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")
            fil.write(
                "   1.25000000e-02   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")  # a point inside the structure
            fil.write(
                "  -3.1415927e+00  -2.7182818e+00   0.0000000e+00   0.0000000e+00\n")  # a point outside the structure

            # SHIFT POINT TO START POINT
            start_point = [-shift, 0]
            fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   3.0000000e+00   0.0000000e+00\n")

            lineTo(start_point, [-shift, Ri_el], step)
            pt = [-shift, Ri_el]
            fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

            # ADD BEAM PIPE LENGTH
            if L_bp_l != 0:
                lineTo(pt, [L_bp_l - shift, Ri_el], step)
                pt = [L_bp_l - shift, Ri_el]
                print(pt)
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

            for n in range(1, n_cell + 1):
                if n == 1:
                    # DRAW ARC:
                    pts = arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el])
                    pt = [-shift + x1el, y1el]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW LINE CONNECTING ARCS
                    lineTo(pt, [-shift + x2el, y2el], step)
                    pt = [-shift + x2el, y2el]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                    pts = arcTo(L_el + L_bp_l - shift, Req_el - B_el, A_el, B_el, step, pt,
                                [L_bp_l + L_el - shift, Req_el])
                    pt = [L_bp_l + L_el - shift, Req_el]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # flat top
                    lineTo(pt, [L_bp_l + L_el + lft_el - shift, Req_el], step)
                    pt = [L_bp_l + L_el + lft_el - shift, Req_el]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    if n_cell == 1:
                        # EQUATOR ARC TO NEXT POINT
                        # half of bounding box is required,
                        # start is the lower coordinate of the bounding box and end is the upper
                        pts = arcTo(L_el + L_bp_l + lft_el - shift, Req_er - B_er, A_er, B_er, step,
                                    [pt[0], Req_er - B_er],
                                    [L_el + lft_el + L_er - x2er + 2 * L_bp_l - shift, Req_er])
                        pt = [L_el + lft_el + L_er - x2er + 2 * L_bp_l - shift, y2er]
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # STRAIGHT LINE TO NEXT POINT
                        lineTo(pt, [L_el + lft_el + L_er - x1er + 2 * L_bp_l - shift, y1er], step)
                        pt = [L_el + lft_el + L_er - x1er + 2 * L_bp_l - shift, y1er]
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # ARC
                        # half of bounding box is required,
                        # start is the lower coordinate of the bounding box and end is the upper
                        pts = arcTo(L_el + lft_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step,
                                    [pt[0], Ri_er],
                                    [L_bp_l + L_el + lft_el + L_er - shift, y1er])

                        pt = [L_bp_l + lft_el + L_el + L_er - shift, Ri_er]
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")

                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # calculate new shift
                        shift = shift - (L_el + L_er + lft_el)
                    else:
                        # EQUATOR ARC TO NEXT POINT
                        # half of bounding box is required,
                        # start is the lower coordinate of the bounding box and end is the upper
                        pts = arcTo(L_el + L_bp_l + lft_el - shift, Req_m - B_m, A_m, B_m, step, [pt[0], Req_m - B_m],
                                    [L_el + L_m + lft_el - x2 + 2 * L_bp_l - shift, Req_m])
                        pt = [L_el + lft_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # STRAIGHT LINE TO NEXT POINT
                        lineTo(pt, [L_el + lft_el + L_m - x1 + 2 * L_bp_l - shift, y1], step)
                        pt = [L_el + lft_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # ARC
                        # half of bounding box is required,
                        # start is the lower coordinate of the bounding box and end is the upper
                        pts = arcTo(L_el + lft_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                    [L_bp_l + L_el + lft_el + L_m - shift, y1])
                        pt = [L_bp_l + L_el + lft_el + L_m - shift, Ri_m]
                        for pp in pts:
                            if (np.around(pp, 12) != np.around(pt, 12)).all():
                                fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                            else:
                                print("Found one")
                        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                        # calculate new shift
                        shift = shift - (L_el + L_m + lft_el)
                        # ic(shift)

                elif n > 1 and n != n_cell:
                    print("elif")
                    # DRAW ARC:
                    pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                    pt = [-shift + x1, y1]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW LINE CONNECTING ARCS
                    lineTo(pt, [-shift + x2, y2], step)
                    pt = [-shift + x2, y2]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                    pts = arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req_m])
                    pt = [L_bp_l + L_m - shift, Req_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # flat top
                    lineTo(pt, [L_bp_l + L_m + lft - shift, Req_m], step)
                    pt = [L_bp_l + L_el + lft - shift, Req_el]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_m + L_bp_l + lft - shift, Req_m - B_m, A_m, B_m, step, [pt[0], Req_m - B_m],
                                [L_m + L_m + lft - x2 + 2 * L_bp_l - shift, Req_m])
                    pt = [L_m + L_m + lft - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_m + L_m + lft - x1 + 2 * L_bp_l - shift, y1], step)
                    pt = [L_m + L_m + lft - x1 + 2 * L_bp_l - shift, y1]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_m + L_m + lft + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                [L_bp_l + L_m + L_m + lft - shift, y1])
                    pt = [L_bp_l + L_m + L_m + lft - shift, Ri_m]
                    ic(pt)
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - 2 * L_m - lft
                else:
                    print("else")
                    # DRAW ARC:
                    pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                    pt = [-shift + x1, y1]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW LINE CONNECTING ARCS
                    lineTo(pt, [-shift + x2, y2], step)
                    pt = [-shift + x2, y2]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                    pts = arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req_m])
                    pt = [L_bp_l + L_m - shift, Req_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # flat top
                    lineTo(pt, [L_bp_l + L_m + lft_er - shift, Req_m], step)
                    pt = [L_bp_l + L_el + lft_er - shift, Req_el]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_m + lft_er + L_bp_l - shift, Req_er - B_er, A_er, B_er, step, [pt[0], Req_er - B_er],
                                [L_m + L_er + lft_er - x2er + L_bp_l + L_bp_r - shift, Req_er])
                    pt = [L_m + L_er + lft_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_m + L_er + lft_er - x1er + L_bp_l + L_bp_r - shift, y1er], step)
                    pt = [L_m + L_er + lft_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_m + L_er + lft_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                [L_bp_l + L_m + L_er + lft_er - shift, y1er])
                    pt = [L_bp_l + L_m + L_er + lft_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

            # BEAM PIPE
            # reset shift
            print("pt before", pt)
            shift = (L_bp_r + L_bp_l + L_el + lft_el + (n_cell - 1) * 2 * L_m + (n_cell - 2) * lft + L_er + lft_er) / 2
            lineTo(pt, [
                L_bp_r + L_bp_l + 2 * (n_cell - 1) * L_m + (n_cell - 2) * lft + lft_el + lft_er + L_el + L_er - shift,
                Ri_er], step)

            if L_bp_r != 0:
                pt = [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r + (
                        n_cell - 2) * lft + lft_el + lft_er - shift, Ri_er]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   3.0000000e+00   0.0000000e+00\n")
                print("pt after", pt)

            # END PATH
            lineTo(pt, [
                2 * (n_cell - 1) * L_m + L_el + L_er + (n_cell - 2) * lft + lft_el + lft_er + L_bp_l + L_bp_r - shift,
                0], step)  # to add beam pipe to right
            pt = [2 * (n_cell - 1) * L_m + L_el + L_er + (n_cell - 2) * lft + lft_el + lft_er + L_bp_l + L_bp_r - shift,
                  0]
            # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
            # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
            fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

            # CLOSE PATH
            lineTo(pt, start_point, step)
            fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

        plt.show()

    def plot_ql_vs_pin(self, Eacc_list=None, mesh=False, Vrfin=None, hist=None):
        """
        Plot loaded quality factor versus input power for fundamental power coupler (FPC)

        Returns
        -------

        """
        label = [cav.name for cav in self.cavities_list]

        # geometry
        n_cells = [cav.n_cells for cav in self.cavities_list]  # 4
        l_cell = [cav.l_cell_mid for cav in self.cavities_list]  # m
        G = [cav.G for cav in self.cavities_list]  # 273.2 # 0.00948*Q_factor*(f/1300)**0.5
        b = [cav.b for cav in self.cavities_list]
        geometry = [n_cells, l_cell, G, b]

        # QOI
        R_Q = [cav.R_Q for cav in self.cavities_list]  # 411   # c Ohm linac definition
        f0 = [cav.op_freq for cav in self.cavities_list]
        QOI = [f0, R_Q]

        # RF
        # Vrf = [2*0.1e9, 2*0.1e9, 2*0.75e9]  #   #2*0.75e9
        # Vrf = [2*0.12e9, 2*0.12e9, 2*1e9, 2*0.44e9]  #   #2*0.75e9 update
        Vrf = [cav.v_rf for cav in self.cavities_list]  # #ttbar
        # Eacc = [20e6, 20e6, 20e6]
        if Eacc_list is None:
            Eacc = [cav.op_field for cav in self.cavities_list]  # update
        else:
            Eacc = Eacc_list
        RF = [Eacc, Vrf]

        # MACHINE
        # I0 = [1390e-3, 1390e-3, 147e-3, 147e-3]  # mA
        I0 = [WP[cav.wp]['I0 [mA]'] * 1e-3 for cav in self.cavities_list]
        # rho = [10.76e3, 10.76e3, 10.76e3, 10.76e3]  # bending radius
        rho = [MACHINE[self.machine]['rho [m]'] for cav in self.cavities_list]  # bending radius
        E0 = [WP[cav.wp]['E [GeV]'] for cav in self.cavities_list]  # Beam energy GeV
        machine = [I0, rho, E0]

        if mesh:
            self.ql_pin_mesh(label, geometry, RF, QOI, machine, Eacc=Eacc_list, Vrf=Vrfin, hist=hist)
            # self.ql_pwp_mesh(label, geometry, RF, QOI, machine, Eacc=Eacc_list, Vrf=Vrfin, hist=hist)
            # self.ql_pin_pwp_ql_mesh(label, geometry, RF, QOI, machine, Eacc=Eacc_list, Vrf=Vrfin, hist=hist)
        else:
            self.ql_pin(label, geometry, RF, QOI, machine)

    def ql_pin(self, labels, geometry, RF, QOI, Machine, p_data=None):
        """
        Calculate the value of input power as a function of loaded quality factor

        Parameters
        ----------
        labels: list, array like
            Descriptive labels on matplotlib plot
        geometry: list, array like
            List of grouped geometric input parameters
        RF: list, array like
            List of grouped radio-frequency (RF) properties
        QOI:
            List of quantities of interest for cavities
        Machine:
            List of grouped machine related materials
        p_data:


        Returns
        -------

        """
        # check if entries are of same length

        it = iter(geometry)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        it = iter(RF)
        the_len = len(next(it))
        # if not all(len(l) == the_len for l in it):
        #     raise ValueError('not all lists have same length!')

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

        fig, ax = plt.subplots(figsize=(12, 4))

        # QOI
        f0, R_Q = [np.array(x) for x in QOI]

        # Machine
        I0, rho, E0 = [np.array(x) for x in Machine]

        l_active = 2 * n_cells * l_cells
        l_cavity = l_active + 8 * l_cells

        # CALCULATED
        v_cav = E_acc * l_active

        rc = PARTICLE[self.particle]['rc [m]']
        m = PARTICLE[self.particle]['m [kg]']
        m_GeV = m * c0 ** 2 * 6241506479.9632  # [GeV]

        C_y = (4 * np.pi / 3) * rc / m_GeV ** 3
        ic(C_y)
        U_loss = C_y * E0 ** 4 / rho  # GeV # energy lost per turn per beam
        v_loss = U_loss / 1 * 1e9  # V # v loss per beam
        ic(U_loss, v_loss, Vrf)

        phi = np.arccos(v_loss / Vrf)
        delta_f = -R_Q * f0 * I0 * np.sin(phi) / (2 * v_cav)  # optimal df
        QL_0_x = v_cav / (R_Q * I0 * np.cos(phi))  # optimal Q loaded
        Pin_x = 50 / np.array([cav.n_cav[np.in1d(self.E_acc, E_acc).nonzero()[0]] for cav in self.cavities_list]) * 1e6
        # ax.plot(QL_0_x, Pin_x[0] * 1e-3, marker='o')
        ic(E0, rho, v_loss, Vrf, v_loss / Vrf, f0, delta_f, I0, np.sin(phi), v_cav, QL_0_x, Pin_x * 1e-3)

        QL_0 = np.linspace(1e4, 1e8, 1000000)

        xy_list = [(0.15, 0.13), (0.1, 0.16), (0.1, 0.19), (0.1, 0.21)]

        if len(E_acc) == len(f0):
            for i in range(len(E_acc)):
                f1_2 = f0[i] / (2 * QL_0)  # 380.6
                ic(R_Q[i], v_cav[i], Vrf[i])
                pin = v_cav[i] ** 2 / (4 * R_Q[i] * QL_0) \
                      * ((1 + ((R_Q[i] * QL_0 * I0[i]) / v_cav[i]) * np.cos(phi[i])) ** 2
                         + ((delta_f[i] / f1_2) + ((R_Q[i] * QL_0 * I0[i]) / v_cav[i]) * np.sin(phi[i])) ** 2)

                # material/ wall power
                e_acc = np.linspace(0.5, 25, 1000) * 1e6  # MV/m

                txt = labels[i]

                if "*" in labels[i]:
                    l = ax.plot(QL_0, pin * 1e-3,
                                label=f"{txt}" + "($Q_\mathrm{L}^* = $" + f"{QL_0[np.argmin(pin)]:.2e})" + f"{round(min(pin) * 1e-3, 1)} kW",
                                lw=4,
                                ls='--')
                else:
                    l = ax.plot(QL_0, pin * 1e-3,
                                label=f"{txt}" + "($Q_\mathrm{L}^* = $" + f"{QL_0[np.argmin(pin)]:.2e})" + f"{round(min(pin) * 1e-3, 1)} kW")

                # add annotations
                ic(l_active)

                # annotext = ax.annotate(txt, xy=xy_list[i], xycoords='figure fraction', size=8, rotation=0,
                #                        c=l[0].get_color())
        else:
            for i in range(len(E_acc)):
                f1_2 = f0[0] / (2 * QL_0)  # 380.6
                # ic(R_Q[0], v_cav[i], Vrf[0])
                pin = v_cav[i] ** 2 / (4 * R_Q[0] * QL_0) \
                      * ((1 + ((R_Q[0] * QL_0 * I0[0]) / v_cav[i]) * np.cos(phi[0])) ** 2
                         + ((delta_f[i] / f1_2) + ((R_Q[0] * QL_0 * I0[0]) / v_cav[i]) * np.sin(phi[0])) ** 2)

                # material/ wall power
                e_acc = np.linspace(0.5, 25, 1000) * 1e6  # MV/m

                txt = labels[0]

                if "*" in labels[0]:
                    l = ax.plot(QL_0, pin * 1e-3,
                                label=f"{txt}" + "($Q_\mathrm{L}^* = $" + f"{QL_0[np.argmin(pin)]:.2e})" + f"{round(min(pin) * 1e-3, 1)} kW",
                                lw=4,
                                ls='--')
                    # plot min q
                    pmin = min(pin)
                    qmin = QL_0[np.argmin(pin)]
                    ax.scatter(qmin, pmin)
                else:
                    # l = ax.plot(QL_0, pin * 1e-3,
                    #             label=f"{txt}" + "($Q_\mathrm{L}^* = $" + f"{QL_0[np.argmin(pin)]:.2e})" + f"{round(min(pin) * 1e-3, 1)} kW")
                    # plot min q
                    pmin = min(pin)
                    qmin = QL_0[np.argmin(pin)]
                    # print(E_acc[i] * 1e-6, v_cav[i] * 1e-6, qmin * 1e-5, pmin * 1e-3)
                    # print([cav.n_cav[np.where(self.E_acc == E_acc[i])] for cav in self.cavities_list])
                    # print(self.E_acc, E_acc[i])
                    ax.scatter(qmin, pmin * 1e-3, ec='k')

                # add annotations
                # ic(l_active)

                # annotext = ax.annotate(txt, xy=xy_list[i], xycoords='figure fraction', size=8, rotation=0,
                #                        c=l[0].get_color())

            # fill between
            pin_low = v_cav[0] ** 2 / (4 * R_Q[0] * QL_0) \
                      * ((1 + ((R_Q[0] * QL_0 * I0[0]) / v_cav[0]) * np.cos(phi[0])) ** 2
                         + ((delta_f[0] / f1_2) + ((R_Q[0] * QL_0 * I0[0]) / v_cav[0]) * np.sin(phi[0])) ** 2)
            pin_high = v_cav[-1] ** 2 / (4 * R_Q[0] * QL_0) \
                       * ((1 + ((R_Q[0] * QL_0 * I0[0]) / v_cav[-1]) * np.cos(phi[0])) ** 2
                          + ((delta_f[-1] / f1_2) + ((R_Q[0] * QL_0 * I0[0]) / v_cav[-1]) * np.sin(phi[0])) ** 2)
            ax.fill_between(QL_0, pin_low * 1e-3, pin_high * 1e-3, zorder=0)

            # input power technological limit
            ax.axhline(1e3, lw=3)
        if p_data:
            # plot QL with penetration
            ax_2 = ax.twinx()
            data = fr.excel_reader(p_data)
            data_ = data[list(data.keys())[0]]
            ax_2.plot(data_["QL"], data_["penetration"], lw=4)

        # plot decorations
        ax.set_xlabel(r"$Q_\mathrm{L}$")
        ax.set_ylabel(r"$P_\mathrm{in} ~[\mathrm{kW}]$")
        ax.set_xscale('log')
        # ax.set_xlim(5e3, 1e10)
        # ax.set_ylim(0, 1200)
        ax.legend(loc='upper left')  #
        ax.minorticks_on()
        # ax.grid(which='both')
        fig.tight_layout()
        fig.show()

    def ql_pin_mesh(self, labels, geometry, RF, QOI, Machine, p_data=None, Vrf=None, Eacc=None, hist=None):
        """
        Calculate the value of input power as a function of loaded quality factor

        Parameters
        ----------
        labels: list, array like
            Descriptive labels on matplotlib plot
        geometry: list, array like
            List of grouped geometric input parameters
        RF: list, array like
            List of grouped radio-frequency (RF) properties
        QOI:
            List of quantities of interest for cavities
        Machine:
            List of grouped machine related materials
        p_data:


        Returns
        -------

        """
        # check if entries are of same length

        it = iter(geometry)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        it = iter(RF)
        the_len = len(next(it))
        # if not all(len(l) == the_len for l in it):
        #     raise ValueError('not all lists have same length!')

        it = iter(QOI)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        it = iter(Machine)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        n_cells, l_cells, G, b = [np.array(x) for x in geometry]
        E_acc, Vrf_ = [np.array(x) for x in RF]

        fig, ax = plt.subplots()

        # QOI
        f0, R_Q = [np.array(x) for x in QOI]

        # Machine
        I0, rho, E0 = [np.array(x) for x in Machine]

        l_active = 2 * n_cells * l_cells
        l_cavity = l_active + 8 * l_cells

        rc = PARTICLE[self.particle]['rc [m]']
        m = PARTICLE[self.particle]['m [kg]']
        m_GeV = m * c0 ** 2 * 6241506479.9632  # [GeV]

        C_y = (4 * np.pi / 3) * rc / m_GeV ** 3
        U_loss = C_y * E0 ** 4 / rho  # GeV # energy lost per turn per beam
        v_loss = U_loss / 1 * 1e9  # V # v loss per beam

        ic(Eacc)
        # filter Vrf, for values less than vloss, the computation does not make sense
        Vrf = Vrf[np.where(Vrf >= v_loss[0])]
        # add vloss to Vrf
        if v_loss[0] not in Vrf:
            Vrf = np.insert(Vrf, 0, v_loss[0])
        EACC, VRF = np.meshgrid(Eacc, Vrf)
        ic(EACC)

        # CALCULATED
        v_cav = EACC * l_active
        n_cav = VRF / v_cav

        phi = np.arccos(v_loss / VRF)
        delta_f = -R_Q * f0 * I0 * np.sin(phi) / (2 * v_cav)  # optimal df
        ic(I0, I0 * np.cos(phi), R_Q)
        QL_0_x = v_cav / (R_Q * I0 * np.cos(phi))  # optimal Q loaded
        Pin_x = 50 / n_cav * 1e6
        # ax.plot(QL_0_x, Pin_x[0] * 1e-3, marker='o')
        # ic(EACC, VRF, QL_0_x, Pin_x)

        # cvrf = ax.contour(QL_0_x, Pin_x * 1e-3, VRF*1e-9, levels=Vrf*1e-9, colors='r')
        # ceacc = ax.contour(QL_0_x, Pin_x * 1e-3, EACC*1e-6, levels=Eacc*1e-6, colors='k')

        cvrf = ax.contour(EACC * 1e-6, VRF * 1e-9, Pin_x * 1e-3, levels=100, colors='r')
        ceacc = ax.contour(EACC * 1e-6, VRF * 1e-9, QL_0_x * 1e-6, levels=100, colors='k')
        ax.clabel(cvrf, inline=True, fontsize=10, colors='r', zorder=50, fmt='%.0f', use_clabeltext=True, manual=True)
        ax.clabel(ceacc, inline=True, fontsize=10, colors='k', zorder=50, use_clabeltext=True, manual=True)

        h1, _ = cvrf.legend_elements()
        h2, _ = ceacc.legend_elements()

        # plot vloss line
        # ax.plot(QL_0_x[0, :], Pin_x[0, :]*1e-3, lw=3, zorder=46, label='$V_\mathrm{RF} = V_\mathrm{loss}$')
        #
        if hist:
            E_acc_hist, Vrf_hist = hist
            # selected points history
            v_cav_hist = E_acc_hist * l_active
            # n_cav_hist = Vrf_hist / v_cav_hist
            #
            # phi_hist = np.arccos(v_loss / Vrf_hist)
            # delta_f = -R_Q * f0 * I0 * np.sin(phi_hist) / (2 * v_cav_hist)  # optimal df
            # QL_0_x = v_cav_hist / (R_Q * I0 * np.cos(phi_hist))  # optimal Q loaded
            # Pin_x = 50/n_cav_hist * 1e6

            hist_labels = ['2018', '2020', '2022']
            handles = []
            for i, (ee, vv) in enumerate(zip(E_acc_hist, Vrf_hist)):
                ax.scatter(ee * 1e-6, vv * 1e-9, s=80, marker='o', ec='k', lw=2, label=hist_labels[i], zorder=100)

        handles, hlabels = plt.gca().get_legend_handles_labels()
        leg = ax.legend([h1[0], h2[0], *handles],
                        [r"$P^*_\mathrm{in} ~[\mathrm{kW}]$", r"$Q^*_\mathrm{L, FM}$ [1E6]", *hlabels],
                        bbox_to_anchor=(1.0, 1), loc='upper left')
        # leg = ax.legend([h1[0], h2[0], *handles], ['$V_\mathrm{RF}$ [GV]', '$E_\mathrm{acc}$ [MV/m]', *hlabels])
        leg.get_frame().set_edgecolor('k')
        # plot decorations
        ax.set_xlabel('$E_\mathrm{acc}$ [MV/m]')
        ax.set_ylabel('$V_\mathrm{RF}$ [GV]')
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        # ax.set_xlim(5e3, 5e6)
        # ax.set_ylim(0, 1200)
        # ax.legend(loc='upper left')  #
        ax.minorticks_on()
        # ax.tick_params(which='minor', length=5, grid_color='k', grid_alpha=0.3)
        # ax.tick_params(which='major', length=10, grid_color='k', grid_alpha=0.3)
        ax.grid(which='both', alpha=0.2)

        fig.tight_layout()
        fig.show()

    def ql_pin_pwp_ql_mesh(self, labels, geometry, RF, QOI, Machine, p_data=None, Vrf=None, Eacc=None, hist=None):
        """
        Calculate the value of input power as a function of loaded quality factor

        Parameters
        ----------
        labels: list, array like
            Descriptive labels on matplotlib plot
        geometry: list, array like
            List of grouped geometric input parameters
        RF: list, array like
            List of grouped radio-frequency (RF) properties
        QOI:
            List of quantities of interest for cavities
        Machine:
            List of grouped machine related materials
        p_data:


        Returns
        -------

        """
        # check if entries are of same length

        it = iter(geometry)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        it = iter(RF)
        the_len = len(next(it))
        # if not all(len(l) == the_len for l in it):
        #     raise ValueError('not all lists have same length!')

        it = iter(QOI)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        it = iter(Machine)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        n_cells, l_cells, G, b = [np.array(x) for x in geometry]
        E_acc, Vrf_ = [np.array(x) for x in RF]

        fig, ax = plt.subplots()

        # QOI
        f0, R_Q = [np.array(x) for x in QOI]

        # Machine
        I0, rho, E0 = [np.array(x) for x in Machine]

        l_active = 2 * n_cells * l_cells
        l_cavity = l_active + 8 * l_cells

        rc = PARTICLE[self.particle]['rc [m]']
        m = PARTICLE[self.particle]['m [kg]']
        m_GeV = m * c0 ** 2 * 6241506479.9632  # [GeV]

        C_y = (4 * np.pi / 3) * rc / m_GeV ** 3
        U_loss = C_y * E0 ** 4 / rho  # GeV # energy lost per turn per beam
        v_loss = U_loss / 1 * 1e9  # V # v loss per beam

        ic(Eacc)
        # filter Vrf, for values less than vloss, the computation does not make sense
        Vrf = Vrf[np.where(Vrf >= v_loss[0])]
        # add vloss to Vrf
        if v_loss[0] not in Vrf:
            Vrf = np.insert(Vrf, 0, v_loss[0])
        EACC, VRF = np.meshgrid(Eacc, Vrf)
        ic(EACC)

        # CALCULATED
        v_cav = EACC * l_active
        n_cav = VRF / v_cav

        phi = np.arccos(v_loss / VRF)
        delta_f = -R_Q * f0 * I0 * np.sin(phi) / (2 * v_cav)  # optimal df
        QL_0_x = v_cav / (R_Q * I0 * np.cos(phi))  # optimal Q loaded
        Pin_x = 50 / n_cav * 1e6
        QL_0 = np.linspace(1e4, 1e8, 1000000)

        xy_list = [(0.15, 0.13), (0.1, 0.16), (0.1, 0.19), (0.1, 0.21)]

        f1_2 = f0[0] / (2 * QL_0)  # 380.6

        # cvrf = ax.contour(QL_0_x, Pin_x * 1e-3, VRF*1e-9, levels=Vrf*1e-9, colors='r')
        # ceacc = ax.contour(QL_0_x, Pin_x * 1e-3, EACC*1e-6, levels=Eacc*1e-6, colors='k')
        NCAV = VRF / v_cav
        Q0 = 2.7e9
        RQ = 152.8
        PWP = 745 * NCAV * v_cav ** 2 / (RQ * Q0) + NCAV * l_active * 8 / (400 / 500) ** 2

        cvrf = ax.contour(EACC * 1e-6, VRF * 1e-9, Pin_x * 1e-3, levels=50, colors='r')
        ceacc = ax.contour(EACC * 1e-6, VRF * 1e-9, QL_0_x * 1e-6, levels=50, colors='k')
        cpwp = ax.contour(EACC * 1e-6, VRF * 1e-9, PWP * 1e-6, levels=60, colors='b', linestyles='-.')
        ax.clabel(cvrf, inline=True, fontsize=10, colors='r', zorder=50, fmt='%.0f', use_clabeltext=True)
        ax.clabel(ceacc, inline=True, fontsize=10, colors='k', zorder=50, use_clabeltext=True)
        ax.clabel(cpwp, inline=True, fontsize=10, colors='b', zorder=50, use_clabeltext=True)

        h1, _ = cvrf.legend_elements()
        h2, _ = ceacc.legend_elements()
        h3, _ = cpwp.legend_elements()

        # plot vloss line
        # ax.plot(QL_0_x[0, :], Pin_x[0, :]*1e-3, lw=3, zorder=46, label='$V_\mathrm{RF} = V_\mathrm{loss}$')
        #
        if hist:
            E_acc_hist, Vrf_hist = hist
            # selected points history
            v_cav_hist = E_acc_hist * l_active
            # n_cav_hist = Vrf_hist / v_cav_hist
            #
            # phi_hist = np.arccos(v_loss / Vrf_hist)
            # delta_f = -R_Q * f0 * I0 * np.sin(phi_hist) / (2 * v_cav_hist)  # optimal df
            # QL_0_x = v_cav_hist / (R_Q * I0 * np.cos(phi_hist))  # optimal Q loaded
            # Pin_x = 50/n_cav_hist * 1e6

            hist_labels = ['2018', '2020', '2022']
            handles = []
            for i, (ee, vv) in enumerate(zip(E_acc_hist, Vrf_hist)):
                ax.scatter(ee * 1e-6, vv * 1e-9, s=80, marker='o', ec='k', lw=2, label=hist_labels[i], zorder=100)

        handles, hlabels = plt.gca().get_legend_handles_labels()
        leg = ax.legend([h1[0], h2[0], h3[0], *handles],
                        [r"$P^*_\mathrm{in} ~[\mathrm{kW}]$", r"$Q^*_\mathrm{L, FM}$ [1E6]", r"$P_\mathrm{wp}$ [MW]",
                         *hlabels],
                        bbox_to_anchor=(1.0, 1), loc='upper left')
        # leg = ax.legend([h1[0], h2[0], *handles], ['$V_\mathrm{RF}$ [GV]', '$E_\mathrm{acc}$ [MV/m]', *hlabels])
        leg.get_frame().set_edgecolor('k')
        # plot decorations
        ax.set_xlabel('$E_\mathrm{acc}$ [MV/m]')
        ax.set_ylabel('$V_\mathrm{RF}$ [GV]')
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        # ax.set_xlim(5e3, 5e6)
        # ax.set_ylim(0, 1200)
        # ax.legend(loc='upper left')  #
        ax.minorticks_on()
        # ax.tick_params(which='minor', length=5, grid_color='k', grid_alpha=0.3)
        # ax.tick_params(which='major', length=10, grid_color='k', grid_alpha=0.3)
        ax.grid(which='both', alpha=0.2)

        fig.tight_layout()
        fig.show()

    def ql_pwp_mesh(self, labels, geometry, RF, QOI, Machine, p_data=None, Vrf=None, Eacc=None, hist=None):
        """
        Calculate the value of input power as a function of loaded quality factor

        Parameters
        ----------
        labels: list, array like
            Descriptive labels on matplotlib plot
        geometry: list, array like
            List of grouped geometric input parameters
        RF: list, array like
            List of grouped radio-frequency (RF) properties
        QOI:
            List of quantities of interest for cavities
        Machine:
            List of grouped machine related materials
        p_data:


        Returns
        -------

        """
        # check if entries are of same length

        it = iter(geometry)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        it = iter(RF)
        the_len = len(next(it))
        # if not all(len(l) == the_len for l in it):
        #     raise ValueError('not all lists have same length!')

        it = iter(QOI)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        it = iter(Machine)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all lists have same length!')

        n_cells, l_cells, G, b = [np.array(x) for x in geometry]
        E_acc, Vrf_ = [np.array(x) for x in RF]

        fig, ax = plt.subplots()

        # QOI
        f0, R_Q = [np.array(x) for x in QOI]

        # Machine
        I0, rho, E0 = [np.array(x) for x in Machine]

        l_active = 2 * n_cells * l_cells
        l_cavity = l_active + 8 * l_cells

        rc = PARTICLE[self.particle]['rc [m]']
        m = PARTICLE[self.particle]['m [kg]']
        m_GeV = m * c0 ** 2 * 6241506479.9632  # [GeV]

        C_y = (4 * np.pi / 3) * rc / m_GeV ** 3
        U_loss = C_y * E0 ** 4 / rho  # GeV # energy lost per turn per beam
        v_loss = U_loss / 1 * 1e9  # V # v loss per beam

        ic(Eacc)
        # filter Vrf, for values less than vloss, the computation does not make sense
        Vrf = Vrf[np.where(Vrf >= v_loss[0])]
        # add vloss to Vrf
        if v_loss[0] not in Vrf:
            Vrf = np.insert(Vrf, 0, v_loss[0])
        EACC, VRF = np.meshgrid(Eacc, Vrf)

        # CALCULATED
        v_cav = EACC * l_active
        NCAV = VRF / v_cav
        Q0 = 2.7e9
        RQ = 152.8
        PWP = 745 * NCAV * v_cav ** 2 / (RQ * Q0) + NCAV * l_active * 8 / (400 / 500) ** 2
        cvrf = ax.contour(EACC * 1e-6, VRF * 1e-9, PWP * 1e-6, levels=120, colors='r')
        ceacc = ax.contour(EACC * 1e-6, VRF * 1e-9, NCAV, levels=80, colors='k')
        ax.clabel(cvrf, inline=True, fontsize=10, colors='r', zorder=50, fmt='%.2f', use_clabeltext=True, manual=True)
        ax.clabel(ceacc, inline=True, fontsize=10, colors='k', zorder=50, use_clabeltext=True, manual=True)

        h1, _ = cvrf.legend_elements()
        h2, _ = ceacc.legend_elements()

        if hist:
            E_acc_hist, Vrf_hist = hist
            # selected points history
            v_cav_hist = E_acc_hist * l_active
            n_cav_hist = Vrf_hist / v_cav_hist
            print(n_cav_hist)

            hist_labels = ['2018', '2020', '2022']
            handles = []
            for i, (q, p) in enumerate(zip(E_acc_hist * 1e-6, Vrf_hist * 1e-9)):
                ax.scatter(q, p, s=80, marker='o', ec='k', lw=2, label=hist_labels[i], zorder=100)

        handles, hlabels = plt.gca().get_legend_handles_labels()
        leg = ax.legend([h1[0], h2[0], *handles], ['$P_\mathrm{wp}$ [MW]', '$N_\mathrm{cav}$ []', *hlabels],
                        bbox_to_anchor=(1.0, 1), loc='upper left')
        # leg.get_frame().set_edgecolor('k')
        # plot decorations
        ax.set_xlabel('$E_\mathrm{acc}$ [MV/m]')
        ax.set_ylabel('$V_\mathrm{RF}$ [GV]')
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        # ax.set_xlim(5e3, 5e6)
        # ax.set_ylim(0, 1200)
        # ax.legend(loc='upper left')  #
        ax.minorticks_on()
        # ax.tick_params(which='minor', length=5, grid_color='k', grid_alpha=0.3)
        # ax.tick_params(which='major', length=10, grid_color='k', grid_alpha=0.3)
        # ax.grid(which='both', alpha=0.2)

        fig.tight_layout()
        fig.show()

    def run_slans(self):
        for cav in self.cavities_list:
            cav.run_slans()

    def run_abci(self):
        for cav in self.cavities_list:
            cav.run_abci()

    def run_multipacting(self):
        for cav in self.cavities_list:
            cav.run_multipacting()

    @staticmethod
    def linspace(start, stop, step=1.):
        """
        Like np.linspace but uses step instead of num
        This is inclusive to stop, so if start=1, stop=3, step=0.5
        Output is: array([1., 1.5, 2., 2.5, 3.])
        """
        if start < stop:
            ll = np.linspace(start, stop, int((stop - start) / abs(step) + 1))
            if stop not in ll:
                ll = np.append(ll, stop)

            return ll
        else:
            ll = np.linspace(stop, start, int((start - stop) / abs(step) + 1))
            if start not in ll:
                ll = np.append(ll, start)
            return ll

    @staticmethod
    def lineTo(prevPt, nextPt, step):
        if prevPt[0] == nextPt[0]:
            # vertical line
            # chwxk id nextPt is greater
            if prevPt[1] < nextPt[1]:
                py = np.linspace(prevPt[1], nextPt[1], step)
            else:
                py = np.linspace(nextPt[1], prevPt[1], step)
                py = py[::-1]
            px = np.ones(len(py)) * prevPt[0]

        elif prevPt[1] == nextPt[1]:
            # horizontal line
            if prevPt[0] < nextPt[1]:
                px = np.linspace(prevPt[0], nextPt[0], step)
            else:
                px = np.linspace(nextPt[0], prevPt[0], step)

            py = np.ones(len(px)) * prevPt[1]
        else:
            # calculate angle to get appropriate step size for x and y
            ang = np.arctan((nextPt[1] - prevPt[1]) / (nextPt[0] - prevPt[0]))
            if prevPt[0] < nextPt[0] and prevPt[1] < nextPt[1]:
                px = np.arange(prevPt[0], nextPt[0], step * np.cos(ang))
                py = np.arange(prevPt[1], nextPt[1], step * np.sin(ang))
            elif prevPt[0] > nextPt[0] and prevPt[1] < nextPt[1]:
                px = np.arange(nextPt[0], prevPt[0], step * np.cos(ang))
                px = px[::-1]
                py = np.arange(prevPt[1], nextPt[1], step * np.sin(ang))
            elif prevPt[0] < nextPt[0] and prevPt[1] > nextPt[1]:
                px = np.arange(prevPt[0], nextPt[0], step * np.cos(ang))
                py = np.arange(nextPt[1], prevPt[1], step * np.sin(ang))
                py = py[::-1]
            else:
                px = np.arange(nextPt[0], prevPt[0], step * np.cos(ang))
                px = px[::-1]
                py = np.arange(nextPt[1], prevPt[1], step * np.sin(ang))
                py = py[::-1]

        # plt.plot(px, py)

    @staticmethod
    def arcTo2(x_center, y_center, a, b, step, start_angle, end_angle):
        u = x_center  # x-position of the center
        v = y_center  # y-position of the center
        a = a  # radius on the x-axis
        b = b  # radius on the y-axis
        sa = (start_angle / 360) * 2 * np.pi  # convert angle to radians
        ea = (end_angle / 360) * 2 * np.pi  # convert angle to radians

        if ea < sa:
            # end point of curve
            x_end, y_end = u + a * np.cos(sa), v + b * np.sin(sa)

            t = np.arange(ea, sa, np.pi / 100)
            # t = np.linspace(ea, sa, 100)
            # check if end angle is included, include if not
            if sa not in t:
                t = np.append(t, sa)
            t = t[::-1]
        else:
            # end point of curve
            x_end, y_end = u + a * np.cos(ea), v + b * np.sin(ea)

            t = np.arange(sa, ea, np.pi / 100)
            # t = np.linspace(ea, sa, 100)
            if ea not in t:
                t = np.append(t, ea)

        # print("t0 ", [(u + a * np.cos(t))[0], (v + b * np.sin(t))[0]])
        # ic([u + a * np.cos(t), v + b * np.sin(t)])
        # ic()

        # plt.plot(u + a * np.cos(t), v + b * np.sin(t))

        return [x_end, y_end]

    @staticmethod
    def arcTo(x_center, y_center, a, b, step, start, end):
        u = x_center  # x-position of the center
        v = y_center  # y-position of the center
        a = a  # radius on the x-axis
        b = b  # radius on the y-axis

        t = np.arange(0, 2 * np.pi, np.pi / 100)

        x = u + a * np.cos(t)
        y = v + b * np.sin(t)
        pts = np.column_stack((x, y))
        inidx = np.all(np.logical_and(np.array(start) < pts, pts < np.array(end)), axis=1)
        inbox = pts[inidx]
        inbox = inbox[inbox[:, 0].argsort()]

        # plt.plot(inbox[:, 0], inbox[:, 1])

        return inbox

    def make_latex_summary_tables(self):
        try:
            l1 = r"\begin{table}[!htb]"
            l2 = r"\centering"
            l3 = r"\caption{Geometric parameters and QoIs of cavities.}"
            l4 = r"\resizebox{\textwidth}{!}{\begin{tabular}{|l|" + f"{''.join(['c' for i in self.cavities_list])}" + "|}"
            toprule = r"\toprule"
            header = r" ".join([fr"& {cav.name} " for cav in self.cavities_list]) + r" \\"
            hline = r"\hline"
            hhline = r"\hline \hline"
            A = r"$A$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][0], 2)}/{round(cav.d_geom_params['OC'][0], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            B = r"$B$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][1], 2)}/{round(cav.d_geom_params['OC'][1], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            a = r"$a$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][2], 2)}/{round(cav.d_geom_params['OC'][2], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            b = r"$b$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][3], 2)}/{round(cav.d_geom_params['OC'][3], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            Ri = r"$R_\mathrm{i}$ " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][4], 2)}/{round(cav.d_geom_params['OC'][4], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            L = r"$L$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][5], 2)}/{round(cav.d_geom_params['OC'][5], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            Req = r"$R_\mathrm{eq}$ [mm] " + "".join(
                [fr"& {round(cav.d_geom_params['IC'][6], 2)}/{round(cav.d_geom_params['OC'][6], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            alpha = r"$ \alpha [^\circ]$" + "".join(
                [fr"& {round(cav.d_geom_params['IC'][7], 2)}/{round(cav.d_geom_params['OC'][7], 2)} " for cav in
                 self.cavities_list]) + r" \\"
            rf_freq = r"RF Freq. [MHz] " + "".join(
                [fr"& {round(cav.op_freq * 1e-6, 2)} " for cav in self.cavities_list]) + r" \\"
            Vrf = r"$V_\mathrm{RF}$ [V] " + "".join(
                [fr"& {round(cav.v_rf, 2):.2E} " for cav in self.cavities_list]) + r" \\"
            Eacc = r"$E_\mathrm{acc}$ [MV/m] " + "".join(
                [fr"& {round(cav.op_field * 1e-6, 2)} " for cav in self.cavities_list]) + r" \\"
            R_Q = r"$R/Q [\Omega$] " + "".join([fr"& {round(cav.R_Q, 2)} " for cav in self.cavities_list]) + r" \\"
            G = r"$G$ [$\Omega$] " + "".join([fr"& {round(cav.G, 2)} " for cav in self.cavities_list]) + r" \\"
            GR_Q = r"$G.R/Q [10^4\Omega^2]$ " + "".join(
                [fr"& {round(cav.GR_Q * 1e-4, 2)} " for cav in self.cavities_list]) + r" \\"
            epk = r"$E_{\mathrm{pk}}/E_{\mathrm{acc}}$ " + "".join(
                [fr"& {round(cav.e, 2)} " for cav in self.cavities_list]) + r" \\"
            bpk = r"$B_{\mathrm{pk}}/E_{\mathrm{acc}} \left[\mathrm{\frac{mT}{MV/m}}\right]$ " + "".join(
                [fr"& {round(cav.b, 2)} " for cav in self.cavities_list]) + r" \\"

            kfm = r"$|k_\mathrm{FM}|$ [V/pC] " + "".join(
                [fr"& {'/'.join([str(round(k_fm, 4)) for k_fm in cav.k_fm])} " for cav in self.cavities_list]) + r" \\"
            kloss = r"$|k_\mathrm{\parallel}|$ [V/pC]" + "".join(
                [fr"& {'/'.join([str(round(k_loss, 4)) for k_loss in cav.k_loss])} " for cav in
                 self.cavities_list]) + r" \\"
            kkick = r"$k_\mathrm{\perp}$ [V/pC/m]" + "".join(
                [fr"& {'/'.join([str(round(k_kick, 4)) for k_kick in cav.k_kick])} " for cav in
                 self.cavities_list]) + r" \\"

            Ncav = r"$N_\mathrm{cav}$ " + "".join(
                [fr"& {int(cav.n_cav_op_field)} " for cav in self.cavities_list]) + r" \\"
            Q0 = r"$Q_\mathrm{0}$ " + "".join(
                [fr"& {int(cav.Q0_desired):2E} " for cav in self.cavities_list]) + r" \\"
            Pin = r"$P_\mathrm{in}\mathrm{/cav} [\mathrm{kW}]$ " + "".join(
                [fr"& {round(cav.p_in * 1e-3, 2)} " for cav in self.cavities_list]) + r" \\"

            Pstat = r"$P_\mathrm{stat}$/cav [W] " + "".join(
                [fr"& {round(cav.pstat, 2)} " for cav in self.cavities_list]) + r" \\"

            Pdyn = r"$P_\mathrm{dyn}$/cav [W] " + "".join(
                [fr"& {round(cav.pdyn, 2)} " for cav in self.cavities_list]) + r" \\"

            Pwp = r"$P_\mathrm{wp}$/cav [kW] " + "".join(
                [fr"& {round(cav.p_wp * 1e-3, 2)} " for cav in self.cavities_list]) + r" \\"

            Phom = r"$P_\mathrm{HOM}$/cav [kW] " + "".join(
                [fr"& {'/'.join([str(round(phom, 4)) for phom in cav.phom])} " for cav in self.cavities_list]) + r" \\"

            phom_values = np.array([[phom for phom in cav.phom] for cav in self.cavities_list])
            ic(phom_values)
            if len(phom_values)%2 == 0:
                result_array = phom_values[1::2] - phom_values[::2]
                ic(result_array)
                Phom_diff = r"$\Delta P_\mathrm{HOM}$/cav [W] & \multicolumn{2}{c|}{" + " & \multicolumn{2}{c|}{".join(
                    [fr"{'/'.join([str(round(val*1e3, 2)) for val in arr])} " + '}' for arr in result_array]) + r" \\"

            PHOM = r"$P_\mathrm{HOM}$ [kW] " + "".join(
                [fr"& {'/'.join([str(round(phom * cav.n_cav_op_field, 2)) for phom in cav.phom])} " for cav in
                 self.cavities_list]) + r" \\"
            bottomrule = r"\bottomrule"
            l34 = r"\end{tabular}}"
            l35 = r"\label{tab: selected shape}"
            l36 = r"\end{table}"

            all_lines = (l1, l2, l3, l4,
                         toprule, hline,
                         header,
                         hhline,
                         A, B, a, b, Ri, L, Req, alpha,
                         hhline,
                         rf_freq, Vrf, Eacc, R_Q, G, GR_Q, epk, bpk,
                         hhline,
                         kfm, kloss, kkick, Phom,
                         hhline, Ncav, Q0, Pin, Pstat, Pdyn, Pwp, Phom, Phom_diff, PHOM,
                         hline,
                         bottomrule,
                         l34, l35, l36)

            # save plots
            fname = [cav.name for cav in self.cavities_list]
            fname = '_'.join(fname)
            ic(fname)
            with open(
                    fr"D:\Dropbox\Quick presentation files\{self.save_folder}\{fname}_latex_summary.txt",
                    'w') as f:
                for ll in all_lines:
                    f.write(ll + '\n')
        except KeyError as e:
            print("Either SLANS or ABCI results not available. Please use '<cav>.set_eigenmode_qois(<folder>)' "
                  "or '<cav>.set_abci_qois(<folder>)' to fix this. Error: ", e)

    def make_excel_summary(self):
        try:
            data = {'Name': [cav.name for cav in self.cavities_list],
                    'Project': [cav.project for cav in self.cavities_list],
                    'Type': [cav.type for cav in self.cavities_list],
                    'CW/Pulsed': [cav.cw_pulsed for cav in self.cavities_list],
                    'Material': [cav.material for cav in self.cavities_list],
                    'N_cells': [cav.n_cells for cav in self.cavities_list],
                    'Freq [MHz]': [cav.op_freq for cav in self.cavities_list],
                    'Beta': [cav.beta for cav in self.cavities_list],
                    'T_oper [K]': [cav.op_temp for cav in self.cavities_list],
                    'I0 [mA]': [cav.I0 for cav in self.cavities_list],
                    'sigma [mm]': [cav.sigma for cav in self.cavities_list],
                    'A_i [mm]': [round(cav.d_geom_params['IC'][0], 2) for cav in self.cavities_list],
                    'B_i [mm]': [round(cav.d_geom_params['IC'][1], 2) for cav in self.cavities_list],
                    'a_i [mm]': [round(cav.d_geom_params['IC'][2], 2) for cav in self.cavities_list],
                    'b_i [mm]': [round(cav.d_geom_params['IC'][3], 2) for cav in self.cavities_list],
                    'R_i [mm]': [round(cav.d_geom_params['IC'][4], 2) for cav in self.cavities_list],
                    'L_i [mm]': [round(cav.d_geom_params['IC'][5], 2) for cav in self.cavities_list],
                    'Req [mm]': [round(cav.d_geom_params['IC'][6], 2) for cav in self.cavities_list],
                    'alpha_i [deg]': [round(cav.d_geom_params['IC'][7], 2) for cav in self.cavities_list],
                    'A_el [mm]': [round(cav.d_geom_params['OC'][0], 2) for cav in self.cavities_list],
                    'B_el [mm]': [round(cav.d_geom_params['OC'][1], 2) for cav in self.cavities_list],
                    'a_el [mm]': [round(cav.d_geom_params['OC'][2], 2) for cav in self.cavities_list],
                    'b_el [mm]': [round(cav.d_geom_params['OC'][3], 2) for cav in self.cavities_list],
                    'R_el [mm]': [round(cav.d_geom_params['OC'][4], 2) for cav in self.cavities_list],
                    'L_el [mm]': [round(cav.d_geom_params['OC'][5], 2) for cav in self.cavities_list],
                    # 'Req [mm]': [round(cav.d_geom_params['OC'][6], 2) for cav in self.cavities_list],
                    'alpha_el [deg]': [round(cav.d_geom_params['OC'][7], 2) for cav in self.cavities_list],
                    'A_er [mm]': [round(cav.d_geom_params['OC'][0], 2) for cav in self.cavities_list],
                    'B_er [mm]': [round(cav.d_geom_params['OC'][1], 2) for cav in self.cavities_list],
                    'a_er [mm]': [round(cav.d_geom_params['OC'][2], 2) for cav in self.cavities_list],
                    'b_er [mm]': [round(cav.d_geom_params['OC'][3], 2) for cav in self.cavities_list],
                    'R_er [mm]': [round(cav.d_geom_params['OC'][4], 2) for cav in self.cavities_list],
                    'L_er [mm]': [round(cav.d_geom_params['OC'][5], 2) for cav in self.cavities_list],
                    # 'Req [mm]': [round(cav.d_geom_params['OC'][6], 2) for cav in self.cavities_list],
                    'alpha_er [deg]': [round(cav.d_geom_params['OC'][7], 2) for cav in self.cavities_list],
                    'R_shunt [Ohm]': ['' for cav in self.cavities_list],
                    'R/Q [Ohm]': [cav.R_Q for cav in self.cavities_list],
                    'k_cc [%]': [cav.k_cc for cav in self.cavities_list],
                    'field flatness [%]': [cav.ff for cav in self.cavities_list],
                    'L_active [m]': [cav.l_active for cav in self.cavities_list],
                    'Epk/Eacc []': [cav.e for cav in self.cavities_list],
                    'Bpk/Eacc [mT/MV/m]': [cav.b for cav in self.cavities_list],
                    'G [Ohm]': [cav.G for cav in self.cavities_list],
                    'R/Q.G [Ohm^2]': [cav.GR_Q for cav in self.cavities_list],
                    '|k_loss| [V/pC]': [cav.k_loss for cav in self.cavities_list],
                    '|k_kick| [V/pC/m]': [cav.k_kick for cav in self.cavities_list],
                    'P_HOM/cav [kW]': [cav.phom for cav in self.cavities_list],
                    'Reference': [cav.reference for cav in self.cavities_list]
                    }

            df = pd.DataFrame.from_dict(data)
            fname = [cav.name for cav in self.cavities_list]
            fname = '_'.join(fname)
            ic(self.save_folder)
            df.to_excel(
                fr"D:\Dropbox\Quick presentation files\{self.save_folder}\{fname}_excel_summary.xlsx",
                sheet_name='Cavities')
        except Exception as e:
            print("Either SLANS or ABCI results not available. Please use '<cav>.set_eigenmode_qois(<folder>)' "
                  "or '<cav>.set_abci_qois(<folder>)' to fix this.")
            print(e)

    def add_cavity(self, cav):
        """
        Adds cavity to cavities
        Parameters
        ----------
        cav: object
            Cavity object

        Returns
        -------

        """
        self.cavities_list.append(cav)

    def remove_cavity(self, cav):
        """
        Removes cavity from cavity list
        Parameters
        ----------
        cav: object
            Cavity object

        Returns
        -------

        """
        self.cavities_list.remove(cav)

    def save_all_plots(self, plot_name):
        """
        Save all plots
        Parameters
        ----------
        plot_name: str
            Name of saved plot

        Returns
        -------

        """
        if self.save_folder != '':
            # check if folder exists
            if os.path.exists(fr"D:\Dropbox\Quick presentation files\{self.save_folder}"):
                save_folder = fr"D:\Dropbox\Quick presentation files\{self.save_folder}"
                plt.savefig(f"{save_folder}/{plot_name}")
            else:
                save_folder = fr"D:\Dropbox\Quick presentation files\{self.save_folder}"
                os.mkdir(save_folder)
                plt.savefig(f"{save_folder}/{plot_name}")

    def __str__(self):
        return fr"{self.cavities_list}"


class Cavity:
    """
    Cavity class defines cavity parameters and operation settings for several cavity analysis
    """

    def __init__(self, slans_dir, abci_dir, vrf, inv_eta=219, name="Unnamed", plot_label=None,
                 op_field=1e6, wp='Z', sigma=[''], op_temp='2K', material='bulkNb', project='', Q0=None,
                 cell_type='normal'):
        """Constructs all the necessary attributes of the Cavity object

        Parameters
        ----------
        n_cells: int
            Number of cells of the cavity
        l_cell_mid: float
            Length of the mid cell of the cavity
        freq: float
            Fundamental mode frequency of the cavity
        vrf: float
            Radio-frequency voltage at which the cavity is designed to operate at
        R_Q: float
            R/Q of the cavity. R/Q is the energy stored to power loss ratio of the cavity.
        G: float
            Geometric factor of the cavity.
        Epk_Eacc: float, list, array like
            Peak electric field to accelerating electric field ratio for the fundamental mode
        Bpk_Eacc
            Peak magnetic field to accelerating electric field ratio for the fundamental mode
        inv_eta: float
            Inverse of eta which is the refrigerating efficiency of the cryomodule
        name: str
            Name given to the cavity. This is the name that will appear in plots.
        op_field:
            Electric field value at which the cavity is to operate.
        wp: str
            Working point of the cavity.
        op_temp: str
            Temperature at which the cavity is designed to operate
        material: str
            Cavity wall material
        """
        # geometric parameters
        # input
        if plot_label is None:
            self.plot_label = name
        else:
            self.plot_label = plot_label

        self.bunch_length = []
        # self.d_fm_all_results = None
        self.slans_dir = slans_dir
        self.abci_dir = abci_dir
        self.ff = ''
        self.k_cc = ''
        self.project = project
        self.reference = ''
        self.beta = 1
        self.cw_pulsed = ''
        self.op_temp = ''
        self.sigma = sigma
        self.op_temp = op_temp
        self.I0 = []
        self.material = material
        self.type = 'Elliptical'
        self.axis_field = None
        self.surface_field = None
        self.Rs = None
        self.d_geom_params = {}
        self.d_qois_fm = {}
        self.k_loss = []
        self.k_kick = []
        self.phom = []
        self.inv_eta = None
        self.name = name
        self.op_field = op_field
        self.p_wp = None
        self.pstat = None
        self.pdyn = None
        self.p_cryo = None
        self.p_in = None
        self.Q0_desired = Q0
        self.Q0 = Q0
        self.l_cavity = None
        self.n_cav = None
        self.n_cav_op_field = None
        self.E_acc = None
        self.n_cells = None
        self.l_cell_mid = None  # m
        self.v_rf = vrf
        self.R_Q = None
        self.k_fm = []
        self.GR_Q = None
        self.op_freq = None  # Hz
        self.e = None
        self.b = None
        self.G = None
        self.wp = wp  # working point
        self.cell_type = cell_type

        self.set_eigenmode_qois(solver='slans')
        self.set_abci_qois(wp)

        # calculated
        self.l_active = 2 * self.n_cells * self.l_cell_mid  # m
        self.n_cav_op_field = int(np.ceil(self.v_rf / (self.op_field * self.l_active)))

        # self.l_cavity = self.n_cav * self.n_cells * self.l_active + (self.n_cav - 1) * 6 * self.l_active

        self.E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m

        # cryo parameters
        # self.eta = 1/219  # 400 MHz
        self.eta = 1 / inv_eta  # 800 MHz

        # power parameters
        # calculated
        self.p_sr = 50e6  # W

    def set_Eacc(self, Eacc):
        """Sets the accelerating field range of analysis for the cavity

        Parameters
        ----------
        Eacc: list, array like
            Accelerating field range of interest

        Returns
        -------

        """
        self.E_acc = Eacc
        self.n_cav = self.v_rf / (self.E_acc * self.l_active)
        self.l_cavity = self.l_active + 8 * self.l_cell_mid  # L_cavity = L_active + 2 λ0

        Rs_dict = {
            "Rs_NbCu_2K_400.79Mhz": 0.57 * (Eacc * 1e-6 * self.b) + 28.4,  # nOhm
            "Rs_NbCu_4.5K_400.79Mhz": 39.5 * np.exp(0.014 * (Eacc * 1e-6 * self.b)) + 27,  # nOhm
            "Rs_bulkNb_2K_400.79Mhz": (2.33 / 1000) * (Eacc * 1e-6 * self.b) ** 2 + 26.24,  # nOhm
            "Rs_bulkNb_4.5K_400.79Mhz": 0.0123 * (Eacc * 1e-6 * self.b) ** 2 + 62.53,  # nOhm

            "Rs_NbCu_2K_801.58Mhz": 1.45 * (Eacc * 1e-6 * self.b) + 92,  # nOhm
            "Rs_NbCu_4.5K_801.58Mhz": 50 * np.exp(0.033 * (Eacc * 1e-6 * self.b)) + 154,  # nOhm
            "Rs_bulkNb_2K_801.58Mhz": (16.4 + Eacc * 1e-6 * self.b * 0.092) * (800 / 704) ** 2,  # nOhm
            "Rs_bulkNb_4.5K_801.58Mhz": 4 * (62.7 + (Eacc * 1e-6 * self.b) ** 2 * 0.012)  # nOhm
        }

        # if self.Q0 is None:
        # if np.isclose(self.op_freq, 801.58e6):
        print(f"Rs_{self.material}_{self.op_temp}_{round(self.op_freq * 1e-6, 2)}Mhz")
        try:
            self.Rs = Rs_dict[f"Rs_{self.material}_{self.op_temp}_{round(self.op_freq * 1e-6, 2)}Mhz"]
        except KeyError:
            self.Rs = Rs_dict[f"Rs_bulkNb_2K_801.58Mhz"]

        self.Q0 = self.G * 1e9 / self.Rs  # c
        # else:
        #     self.Q0 = np.ones(np.shape(Eacc))*self.Q0
        # ic("800 MHz")
        # elif np.isclose(self.op_freq, 400.79e6):
        #     self.Rs = Rs_bulkNb_2k_800Mhz
        #     self.Q0 = self.G * 1e9 / Rs_NbCu_4_5k_400Mhz  # c
        #     # ic("400 MHz")

        self.get_power()

    def get_Eacc(self):
        """Gets the accelerating field range of the Cavity object

        Returns
        -------
            Accelerating field range of cavity. To set the accelerating field range, use

        """
        return self.E_acc

    def set_op_field(self, op_field):
        """Sets the operating field of the cavity

        Parameters
        ----------
        op_field: float
            Electric field value at which the cavity is to operate.

        Returns
        -------

        """
        self.op_field = op_field

    def get_op_field(self):
        """Gets the operating field of the cavity

        Returns
        -------
            returns operating field value

        """
        return self.op_field

    def set_inv_eta(self, inv_eta):
        self.inv_eta = inv_eta
        self.eta = 1 / inv_eta

    def get_inv_eta(self):
        return self.inv_eta

    def set_vrf(self, vrf):
        self.v_rf = vrf

    def get_vrf(self):
        return self.v_rf

    def get_power(self):
        self.p_in = self.p_sr / self.n_cav_op_field  # maximum synchrotron radiation per beam

        # self.p_in = self.v_cav ** 2 / (4 * self.R_Q * Q) * ((1 + self.R_Q * Q * WP[self.wp]['I0'] / self.v_cav) ** 2
        #                                                     + (delta_f / f1_2 + self.R_Q * Q * WP[self.wp][
        #             'I0'] / self.v_cav * np.sin(phi)) ** 2)

        self.p_cryo = 8 / (np.sqrt(self.op_freq / 500e6))  # W/m

        # ic(self.v_rf, self.l_active, self.R_Q, self.l_cavity, self.p_in, self.n_cav_op_field, 1/self.eta, self.Q0)
        self.pdyn = self.v_rf * (self.op_field * self.l_active) / (
                    self.R_Q * self.Q0_desired * self.n_cav_op_field)  # per cavity
        self.pstat = (self.l_cavity * self.v_rf / (
                    self.l_active * self.op_field * self.n_cav_op_field)) * self.p_cryo  # per cavity
        self.p_wp = (1 / self.eta) * (self.pdyn + self.pstat)  # per cavity

        return self.p_in, self.p_wp

    def set_eigenmode_qois(self, solver='slans'):
        qois = 'qois.json'
        geom_params = 'geometric_parameters.json'

        with open(fr"{self.slans_dir}\monopole\{qois}") as json_file:
            self.d_qois_fm.update(json.load(json_file))

        with open(fr"{self.slans_dir}\monopole\{geom_params}") as json_file:
            self.d_geom_params.update(json.load(json_file))

        # if solver.lower() == 'slans':
        #     slans_svl = 'cavity_33.svl'
        #     self.d_fm_all_results = fr.svl_reader(fr"{folder_name}\{slans_svl}")
        # else:
        #     with open(fr"{folder_name}\{qois}") as json_file:
        #         self.d_fm_all_results.update(json.load(json_file))

        # ic(d_qois)
        self.l_cell_mid = self.d_geom_params['IC'][5] * 1e-3
        self.n_cells = self.d_qois_fm['N Cells']
        self.op_freq = self.d_qois_fm['freq [MHz]'] * 1e6
        self.k_cc = self.d_qois_fm['kcc [%]']
        self.ff = self.d_qois_fm['ff [%]']
        self.R_Q = self.d_qois_fm['R/Q [Ohm]']
        self.GR_Q = self.d_qois_fm['GR/Q [Ohm^2]']
        self.G = self.GR_Q / self.R_Q
        # print(G)
        # self.Q = d_qois['Q []']
        self.e = self.d_qois_fm['Epk/Eacc []']
        self.b = self.d_qois_fm['Bpk/Eacc [mT/MV/m]']

        # get axis field
        try:
            self.axis_field = fr.txt_reader(fr"{self.slans_dir}\cavity_33_{self.n_cells}.af", ' ')
        except FileNotFoundError:
            self.axis_field = []
        # print(self.axis_field)

        # get surface field
        try:
            self.surface_field = fr.txt_reader(fr"{self.slans_dir}\cavity_33_{self.n_cells}.sf", ' ')
        except FileNotFoundError:
            self.surface_field = []
        # print(self.surface_field)

    def set_abci_qois(self, working_point=''):

        qois = 'qois.json'
        geom_params = 'geometric_parameters.json'

        d_qois = {}
        with open(fr"{self.abci_dir}\{qois}") as json_file:
            d_qois.update(json.load(json_file))

        d_geom_params = {}
        with open(fr"{self.abci_dir}\{geom_params}") as json_file:
            d_geom_params.update(json.load(json_file))

        # ic(d_geom_params)
        ic(d_qois)
        for sigma in self.sigma:
            d_qoi = d_qois[f'{working_point}_{sigma}']

            self.bunch_length.append(d_qoi["sigma_z [mm]"])
            self.k_fm.append(d_qoi['k_FM [V/pC]'])
            self.k_loss.append(d_qoi['|k_loss| [V/pC]'])
            self.k_kick.append(d_qoi['|k_kick| [V/pC/m]'])
            self.phom.append(d_qoi['P_HOM [kW]'])
            self.I0.append(d_qoi['I0 [mA]'])

    def plot_multipac_triplot(self, folder, type='triplot'):

        if type == 'triplot':
            # create figure
            fig = plt.figure()
            gs = fig.add_gridspec(3, 1)
            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[1, 0])
            ax3 = fig.add_subplot(gs[2, 0])
            axs = [ax1, ax2, ax3]

        else:
            fig, axs = plt.subplots(1, 1)

        mpl.rcParams['figure.figsize'] = [6, 10]

        Eacc, Epk_Eacc, label = self.op_field, self.e, self.name

        # load_output_data
        # files
        fnames = ["Ccounter.mat", "Acounter.mat", "Atcounter.mat", "Efcounter.mat", "param",
                  "geodata.n", "secy1", "counter_flevels.mat", "counter_initials.mat"]
        data = {}
        # files_folder = "D:\Dropbox\multipacting\MPGUI21"
        for f in fnames:
            if ".mat" in f:
                data[f] = spio.loadmat(fr"{folder}\\{f}")
            else:
                data[f] = pd.read_csv(fr"{folder}\\{f}", sep='\s+', header=None)

        A = data["Acounter.mat"]["A"]
        At = data["Atcounter.mat"]["At"]
        C = data["Ccounter.mat"]["C"]
        Ef = data["Efcounter.mat"]["Ef"]
        flevel = data["counter_flevels.mat"]["flevel"]
        initials = data["counter_initials.mat"]["initials"]
        secy1 = data["secy1"].to_numpy()
        Pow = flevel
        n = len(initials[:, 0]) / 2  # number of initials in the bright set
        N = int(data["param"].to_numpy()[4])  # number of impacts
        U = flevel
        Efl = flevel
        q = 1.6021773e-19
        Efq = Ef / q

        e1 = np.min(np.where(secy1[:, 1] >= 1))  # lower threshold
        e2 = np.max(np.where(secy1[:, 1] >= 1))  # upper threshold
        val, e3 = np.max(secy1[:, 1]), np.argmax(secy1[:, 1])  # maximum secondary yield

        cl = 0
        ok, ok1, ok2 = 1, 1, 1
        if ok > 0:
            if n == 0:
                ic('Unable to plot the counters. No initial points.')
                return

            if ok1 * ok2 == 0:
                cl = ic('Counter functions or impact energy missing.')
            else:
                # if ss > 0:
                #     cl = ic(np.array(['Plotting the triplot (counter, enhanced ', 'counter and impact energy).']))

                if type == 'counter function' or type == 'triplot':
                    # fig, axs = plt.subplots(3)
                    axs[0].plot(Efl / 1e6, C / n, lw=2, label=label)
                    axs[0].set_ylabel("$c_" + "{" + f"{N}" + "}/ c_0 $")
                    axs[0].set_xlabel(r'$E_\mathrm{pk}$ [MV/m]')
                    # axs[0].set_title(r'$\mathbf{MultiPac 2.1~~~~~Counter function~~~~}$')
                    axs[0].set_xlim(np.amin(Efl) / 1e6, np.amax(Efl) / 1e6)
                    axs[0].set_ylim(0, np.max([0.1, axs[0].get_ylim()[1]]))

                    # plot peak operating field
                    axs[0].axvline(Eacc * Epk_Eacc, c='k', ls='--', lw=2)
                    axs[0].text(np.round(Eacc * Epk_Eacc, 2) - 1.5, 0.1,
                                f"{label[0]}: {np.round(Eacc * Epk_Eacc, 2)} MV/m",
                                size=12, rotation=90,
                                transform=axs[0].get_xaxis_transform())

                    axs[0].minorticks_on()

                if type == 'final impact energy' or type == 'triplot':
                    axs[1].semilogy(Efl / 1e6, Efq, lw=2)

                    # axs[1].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [secy1[e1, 0], secy1[e1, 0]], '-r')
                    e0 = sci.interp1d(secy1[0:e1 + 1, 1], secy1[0:e1 + 1, 0])(1)
                    axs[1].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [e0, e0], '-r')
                    axs[1].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [secy1[e2, 0], secy1[e2, 0]], '-r')
                    axs[1].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [secy1[e3, 0], secy1[e3, 0]], '--r')

                    axs[1].set_ylabel("$Ef_" + "{" + f"{N}" + "}$")
                    axs[1].set_xlabel(r'$E_\mathrm{pk}$ [MV/m]')
                    # axs[1].set_title('$\mathbf{Final~Impact~Energy~in~eV}$')
                    axs[1].set_xlim(np.min(Efl) / 1e6, np.max(Efl) / 1e6)
                    axs[1].set_ylim(0, axs[1].get_ylim()[1])

                    axs[1].axvline(Eacc * Epk_Eacc, c='k', ls='--', lw=2)
                    axs[1].text(np.round(Eacc * Epk_Eacc, 2) - 1.5, 0.1,
                                f"{label[0]}: {np.round(Eacc * Epk_Eacc, 2)} MV/m",
                                size=12, rotation=90,
                                transform=axs[1].get_xaxis_transform())

                    axs[1].minorticks_on()
                if type == 'enhanced counter function' or type == 'triplot':
                    axs[2].semilogy(Efl / 1e6, (A + 1) / n, lw=2)
                    axs[2].set_xlabel('$V$ [MV]')
                    axs[2].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [1, 1], '-r')
                    axs[2].set_xlim(np.min(Efl) / 1e6, np.max(Efl) / 1e6)
                    axs[2].set_ylim(np.min((A + 1) / n), axs[2].get_ylim()[1])
                    axs[2].set_ylabel("$e_" + "{" + f"{N}" + "}" + "/ c_0$")
                    axs[2].set_xlabel(r'$E_\mathrm{pk}$ [MV/m]')
                    # axs[2].set_title('$\mathbf{Enhanced~counter~function}$')

                    axs[2].axvline(Eacc * Epk_Eacc, c='k', ls='--', lw=2)
                    axs[2].text(np.round(Eacc * Epk_Eacc, 2) - 1, 0.1,
                                f"{label[0]}: {np.round(Eacc * Epk_Eacc, 2)} MV/m",
                                size=12, rotation=90,
                                transform=axs[2].get_xaxis_transform())

                    axs[2].minorticks_on()

        axs[0].legend(loc='upper left')

        fig.tight_layout()
        plt.show()

    def ql_pin(self, labels, geometry, RF, QOI, Machine, p_data=None):
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

        self.fig, ax = plt.subplots()

        self.cid_pick = self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        self.cid_motion = self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.cid_button_press = self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_button_release = self.fig.canvas.mpl_connect('button_release_event', self.on_release)

        # QOI
        f0, R_Q = [np.array(x) for x in QOI]

        # Machine
        I0, rho, E0 = [np.array(x) for x in Machine]

        l_active = 2 * n_cells * l_cells
        l_cavity = l_active + 8 * l_cells

        # CALCULATED
        v_cav = E_acc * l_active

        U_loss = 88.46 * E0 ** 4 / rho * 1e-6  # GeV # energy lost per turn per beam
        v_loss = U_loss * 1e9  # V # v loss per beam

        print(v_loss, Vrf, v_loss / Vrf)
        phi = np.arccos(v_loss / Vrf)
        delta_f = -R_Q * f0 * I0 * np.sin(phi) / (2 * v_cav)  # optimal df
        QL_0_x = v_cav / (R_Q * I0 * np.cos(phi))  # optimal Q loaded

        QL_0 = np.linspace(1e4, 1e9, 1000000)

        xy_list = [(0.7, 0.65), (0.7, 0.3), (0.6, 0.6), (0.77, 0.25)]
        for i in range(len(E_acc)):
            f1_2 = f0[i] / (2 * QL_0)  # 380.6
            ic(R_Q[i], v_cav[i], Vrf[i])
            pin = v_cav[i] ** 2 / (4 * R_Q[i] * QL_0) * \
                  ((1 + ((R_Q[i] * QL_0 * I0[i]) / v_cav[i]) * np.cos(phi[i])) ** 2 +
                   ((delta_f[i] / f1_2) + ((R_Q[i] * QL_0 * I0[i]) / v_cav[i]) * np.sin(phi[i])) ** 2)

            # p_cryo = 8 / (np.sqrt(f0[i] / 500e6))

            # material/ wall power
            e_acc = np.linspace(0.5, 25, 1000) * 1e6  # MV/m

            # Rs_NbCu_4_5k_400Mhz = 39.5 * np.exp(0.014 * (E_acc[i] * 1e-6 * b[i])) + 27

            # eta = 1 / 219  # c
            # Q0 = G[i] * 1e9 / Rs_NbCu_4_5k_400Mhz  # c
            # print(E_acc[i])
            # p_wp = (1 / eta) * Vrf[i] * (E_acc[i] * l_active[i] / (R_Q[i] * Q0)) + (1 / eta) * (
            #         l_cavity[i] * Vrf[i] / (l_active[i] * E_acc[i])) * p_cryo

            if "*" in labels[i]:
                l = ax.plot(QL_0, pin * 1e-3, label=f"${round(E_acc[i] * 1e-6, 2)}" + " ~[\mathrm{MV/m}]$", lw=4,
                            ls='--')
            else:
                l = ax.plot(QL_0, pin * 1e-3, label=f"${round(E_acc[i] * 1e-6, 2)}" + " ~[\mathrm{MV/m}]$", lw=4)

            # add annotations
            ic(l_active)
            txt = f"{labels[i]}, " \
                  f"{n_cells[i]}-Cell {int(f0[i] / 1e6)} MHz {int(np.ceil(Vrf[i] / (E_acc[i] * l_active[i])))} " \
                  f"cav" + " V$_\mathrm{RF}$ =" + f"{round(Vrf[i] * 1e-9, 2)} GV " \
                  + " V$_\mathrm{cav}$ =" + f"{round(v_cav[i] * 1e-6, 1)} MV " \
                                            "P$_{\mathrm{in}}$ = " + f"{round(min(pin) * 1e-3, 1)} kW" \
                                                                     "Q$_{\mathrm{L, 0}}^*$ = " + "{:.2e}".format(
                QL_0_x[i])

            annotext = ax.annotate(txt, xy=xy_list[i], xycoords='figure fraction', size=12, rotation=0,
                                   c=l[0].get_color(),
                                   weight='bold')

            dt = DraggableText(annotext)
            dt.connect()

            self.text_dict[id(annotext)] = dt

        if p_data:
            # plot QL with penetration
            ax_2 = ax.twinx()
            data = fr.excel_reader(p_data)
            data_ = data[list(data.keys())[0]]
            ax_2.plot(data_["QL"], data_["penetration"], lw=4)

        # plot decorations
        ax.set_xlabel(r"$Q_{L,0}$")
        ax.set_ylabel(r"$P_\mathrm{in}/\mathrm{cav} ~[\mathrm{kW}]$")
        ax.set_xscale('log')
        ax.set_xlim(5e3, 1e9)
        ax.set_ylim(100, 2000)
        ax.legend(loc='lower left', title=r"$E_\mathrm{acc}$")  #
        ax.minorticks_on()
        # ax.grid(which='both')
        plt.show()

    def on_motion(self, event):
        # vis = self.annot.get_visible()
        # if isinstance(self.ax, event.inaxes):  # not so good fix
        #     cont, ind = self.plot_object.contains(event)
        #     if cont:
        #         self.update_annot(
        #             ind["ind"][0])  # ind returns an array of close points. ind['ind'][0] returns just the first point
        #         self.annot.set_visible(True)
        #         self.fig.canvas.draw_idle()
        #     else:
        #         if vis:
        #             self.annot.set_visible(False)
        #             self.fig.canvas.draw_idle()
        self.MOTION = True

        # if self.PRESS and self.MOTION:

    def run_slans(self):
        pass

    def run_abci(self):
        pass

    def run_multipacting(self):
        pass

    def plot_ql_vs_pin(self):
        label = self.name

        # geometry
        n_cells = [self.n_cells]
        l_cell = [self.l_cell_mid]
        G = [self.G]
        b = [self.b]
        geometry = [n_cells, l_cell, G, b]

        # QOI
        R_Q = [self.R_Q]
        f0 = [self.op_freq]
        QOI = [f0, R_Q]

        # RF
        # Vrf = [2*0.1e9, 2*0.1e9, 2*0.75e9]  #   #2*0.75e9
        # Vrf = [2*0.12e9, 2*0.12e9, 2*1e9, 2*0.44e9]  #   #2*0.75e9 update
        Vrf = [self.v_rf]
        # Eacc = [20e6, 20e6, 20e6]
        Eacc = [self.op_field]
        RF = [Eacc, Vrf]

        # MACHINE
        # I0 = [1390e-3, 1390e-3, 147e-3, 147e-3]  # mA
        I0 = [WP[self.wp]['I0 [mA]'] * 1e-3]
        # rho = [10.76e3, 10.76e3, 10.76e3, 10.76e3]  # bending radius
        rho = [MACHINE['rho [m]']]  # bending radius
        E0 = [WP[self.wp]['E [GeV]']]  # Beam energy GeV
        machine = [I0, rho, E0]

        self.ql_pin(label, geometry, RF, QOI, machine)

    def view(self):
        # view cavity contour plot
        pass

    def write_cavity_for_multipac(self):
        pass

    def make_latex_summary_tables(self):
        try:
            l1m = r"\begin{table}[!htb]"
            l2m = r"\centering"
            l3m = r"\caption{Geometric parameters of optimized cavity" + f" {self.name}" + ".}"
            l4m = r"\resizebox{\textwidth}{!}{\begin{tabular}{cccccccc}"
            l5m = r"\toprule"
            l6m = r"\multicolumn{8}{c}{" + fr"{self.name}" + r" (Mid-cell) -  Geometric Properties} \\"
            l7m = r"\midrule"
            l8m = r"\midrule"
            l9m = r"$A_\mathrm{m}$ [mm] & $B_\mathrm{m}$ [mm] & $a_\mathrm{m}$ [mm] & $b_\mathrm{m}$ [mm] & " \
                  r"$R_\mathrm{i, m}$ [mm] & $L_\mathrm{m}$ [mm] & $R_{eq, m}$ [mm] & $\alpha_\mathrm{m} [^\circ]$ \\"
            l10m = r"\midrule"
            l11m = fr"{round(self.d_geom_params['IC'][0], 2)}" + \
                   r"".join([fr"& {round(x, 2)}" for i, x in enumerate(self.d_geom_params['IC']) if i > 0]) + r" \\"
            l12m = r"\bottomrule"
            l13m = r"\end{tabular}}"
            l14m = r"\label{tab: geometric cavity end-cell}"
            l15m = r"\end{table}"

            l_break = "\n" * 2

            l1e = r"\begin{table}[!htb]"
            l2e = r"\centering"
            l3e = r"\caption{Geometric parameters of optimized cavity" + f" {self.name}" + ".}"
            l4e = r"\resizebox{\textwidth}{!}{\begin{tabular}{cccccccc}"
            l5e = r"\toprule"
            l6e = r"\multicolumn{8}{c}{" + fr"{self.name}" + r" (End-cell) -  Geometric Properties} \\"
            l7e = r"\midrule"
            l8e = r"\midrule"
            l9e = r"$A_\mathrm{e}$ [mm] & $B_\mathrm{e}$ [mm] & $a_\mathrm{e}$ [mm] & $b_\mathrm{e}$ [mm] & " \
                  r"$R_\mathrm{i, e}$ [mm] & $L$ [mm] & $R_{eq, e}$ [mm] & $\alpha_\mathrm{e} [^\circ]$ \\"
            l10e = r"\midrule"
            l11e = fr"{round(self.d_geom_params['OC'][0], 2)}" + \
                   r"".join([fr"& {round(x, 2)}" for i, x in enumerate(self.d_geom_params['OC']) if i > 0]) + r" \\"
            l12e = r"\bottomrule"
            l13e = r"\end{tabular}}"
            l14e = r"\label{tab: geometric cavity end-cell}"
            l15e = r"\end{table}"

            l1_qois = r"\begin{table}[!htb]"
            l2_qois = r"\centering"
            l3_qois = r"\caption{QoIs of optimized cavity C$_{3795}$.}"
            l4_qois = r"\resizebox{\textwidth}{!}{\begin{tabular}{cccccccc}"
            l5_qois = r"\toprule"
            l6_qois = r"\multicolumn{8}{c}{C$_{3795}$ - QOIs} \\"
            l7_qois = r"\midrule"
            l8_qois = r"\midrule"
            l9_qois = r"No. of cells & $R/Q$  & $G$ & $G.R/Q$ & $E_{\mathrm{pk}}/E_{\mathrm{acc}}$  & " \
                      r"$B_{\mathrm{pk}}/E_{\mathrm{acc}}$  & $|k_\mathrm{\parallel}| \mathrm{[SR]}$ & " \
                      r"$k_\mathrm{\perp} \mathrm{[SR]}$ \\"
            l10_qois = r"& [$\Omega$] &[$\Omega$] & $[10^5\Omega^2]$ & [-] & $[\mathrm{\frac{mT}{MV/m}}]$ & " \
                       r"$[\mathrm{V/pC}]^$& $[\mathrm{V/pC/m}]^$ \\"
            l11_qois = r"\midrule"
            l12_qois = fr"{int(self.n_cells)} & {round(self.R_Q, 2)} & {round(self.G, 2)} & " \
                       fr"{round(self.GR_Q * 1e-5, 2)} & {round(self.e, 2)} & " \
                       fr"{round(self.b, 2)} & {round(self.k_loss, 2)} & {round(self.k_kick, 2)} \\"
            l13_qois = r"\bottomrule"
            l14_qois = r"\end{tabular}}"
            l15_qois = r"\label{tab: qois designed cavity}"
            l16_qois = r"\end{table}"

            all_lines = (l1m, l2m, l3m, l4m, l5m, l6m, l7m, l8m, l9m, l10m,
                         l11m, l12m, l13m, l14m, l15m, l_break,
                         l1e, l2e, l3e, l4e, l5e, l6e, l7e, l8e, l9e, l10e,
                         l11e, l12e, l13e, l14e, l15e, l_break,
                         l1_qois, l2_qois, l3_qois, l4_qois, l5_qois, l6_qois, l7_qois, l8_qois, l9_qois, l10_qois,
                         l11_qois, l12_qois, l13_qois, l14_qois, l15_qois, l16_qois
                         )

            with open(r"D:\Dropbox\Quick presentation files\latex_mid_end_table.txt", 'w') as f:
                for ll in all_lines:
                    f.write(ll + '\n')
        except KeyError as e:
            print("Either SLANS or ABCI results not available. Please use '<cav>.set_eigenmode_qois(<folder>)' "
                  "or '<cav>.set_abci_qois(<folder>)' to fix this.", e)

    def __repr__(self):
        return fr"{self.name}({self.n_cells} [], {self.op_freq * 1e-6} MHz, {self.e} [], {self.b} [mT/MV/m])"

    def __str__(self):
        return fr"{self.name}({self.n_cells} [], {self.op_freq * 1e-6} MHz, {self.e} [], {self.b} [mT/MV/m])"


class Machine:
    def __init__(self, name, rho):
        # machine properties
        self.name = name
        self.rho = rho

        self.E = []
        self.bunch_length = []
        self.I0 = []
        self.op_points = []

    def add_op_points(self, op_point):
        self.op_points.append(op_point)

    def remove_op_points(self, op_point):
        self.op_points.remove(op_point)

    def analyse_machine(self):
        pass

    def compare_cavities_for_op_point(self, op_point):
        pass


class OperatingPoint:
    def __init__(self, op_file_path):

        if isinstance(op_file_path, list):
            name, op_parameters = self.load_from_operating_points(op_file_path)
        else:
            name, op_parameters = self.load_operating_point(op_file_path)

        self.name = name
        self.E = op_parameters["E [GeV]"]
        self.I0 = op_parameters["I0 [mA]"]
        self.vrf = op_parameters["V [GV]"] * 1e9
        self.op_field = op_parameters["Eacc [MV/m]"] * 1e6
        self.op_temp = op_parameters["T [K]"]
        self.sigma_sr = op_parameters["sigma_SR [mm]"]
        self.sigma_bs = op_parameters["sigma_BS [mm]"]
        self.cavities = []

    def analyse(self):
        pass

    @staticmethod
    def load_from_operating_points(op_file_path):
        with open(fr"{op_file_path[0]}") as json_file:
            op_parameters = json.load(json_file)

        return op_file_path[1], op_parameters[op_file_path[1]]

    @staticmethod
    def load_operating_point(op_file_path):
        with open(fr"{op_file_path}") as json_file:
            op_parameters = json.load(json_file)

        return op_parameters.key(), op_parameters

    def set_sigma_sr(self, sigma):
        self.sigma_sr = sigma

    def get_sigma_sr(self):
        return self.sigma_sr

    def set_sigma_bs(self, sigma):
        self.sigma_bs = sigma

    def get_sigma_bs(self):
        return self.sigma_bs

    def set_I0(self, I0):
        self.I0 = I0

    def get_I0(self):
        return self.I0

    def set_energy(self, e):
        self.E = e

    def get_energy(self):
        return self.E

    def add_cavities(self, cavities):
        if isinstance(cavities, list):
            self.cavities.extend(cavities)
        else:
            self.cavities.append(cavities)

    def remove_cavities(self, cavities):
        if isinstance(cavities, list):
            for cav in cavities:
                self.cavities.remove(cav)
        else:
            self.cavities.remove(cavities)

    def set_beam_properties(self):
        pass


op_folder = fr"D:\Dropbox\CavityDesignHub\Cavity800\OperatingPoints\fcc.json"
Z_op_point = OperatingPoint([op_folder, "Z_2022"])
W_op_point = OperatingPoint([op_folder, "W_2022"])
H_op_point = OperatingPoint([op_folder, "H_2022"])
ttbar_op_point = OperatingPoint([op_folder, "ttbar_2022"])

H_op_point.add_cavities([3])


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
    # from labellines import labelLines

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

    sr_labels = ['Rs_Nb/Cu_2k_400MHz', 'Rs_Nb/Cu_4_5k_400MHz', 'Rs_bulkNb_2k_400MHz', 'Rs_bulkNb_4_5k_400MHz',
                 'Rs_Nb/Cu_2k_800MHz', 'Rs_Nb/Cu_4_5k_800MHz', 'Rs_bulkNb_2k_800MHz', 'Rs_bulkNb_4_5k_800MHz']

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
        plt.plot(Eacc, sr, label=sr_labels[i])

    # labelLines(plt.gca().get_lines(), zorder=2.5)
    # plt.legend()
    plt.xlabel("$E$ [MV/m]")
    plt.ylabel(r"$R_s [\mathrm{n\Omega}]$")
    plt.tight_layout()
    plt.yscale('log')
    plt.show()


# plot_surface_resistance()


def plot_brillouin(op_freq, p, cst_result_folder=None):
    fig, axs = plt.subplots(1, p)

    # get dictionary instead from cst run folder
    d = pd.read_csv(fr"D:\CST Studio\5. tt\Eigenmode\E3795_PBC.csv", sep=",", skipfooter=1, engine='python').to_dict(
        orient='list')
    # d = pd.read_csv(fr"D:\CST Studio\5. tt\Eigenmode\E3795_PBC_endcell.csv", sep=",", skipfooter=1,
    # engine='python').to_dict(orient='list')
    freq_dict = {key: val for key, val in d.items() if 'Mode' in key}
    phase = d['phase']

    nn = np.tile(np.array([0, 180]), p)
    nn = [nn[i:i + 2] for i in range(0, len(nn) - 1)]

    freq_list = np.arange(0, p + 3) * op_freq
    freq_list = [freq_list[i:i + 2] for i in range(0, len(freq_list))]

    # light line
    for i, x in enumerate(nn):
        for ax in axs:
            ax.plot(x, freq_list[i], c='k', label='light line')

    for key, mode_freqs in enumerate(freq_dict.values()):
        for i, ax in enumerate(axs):
            ax.plot(np.array(phase), mode_freqs, marker='o', mec='k', lw=3)
            ax.set_ylim(min(freq_list[i + 1]) - 50, max(freq_list[i + 1]))
            axs[i].set_xlim(-3, 183)

    plt.tight_layout()
    plt.legend()
    plt.show()


# plot_brillouin(801.55, 2)


def mucol_study():
    # MUCOL STUDY

    # RCS
    V = 20.87 * 1e9

    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\MuCol_Study\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\MuCol_Study\SimulationData\ABCI"

    ILC_LL = Cavity(vrf=V, inv_eta=745, name="ILC-LL", op_field=30e6, op_temp='2K', material='bulkNb',
                    wp=wp, sigma=sigma, plot_label="ILC-LL",
                    slans_dir=fr"{parent_dir_slans}\ILC-LL",
                    abci_dir=fr"{parent_dir_abci}\ILC-LL")

    ICHIRO = Cavity(vrf=V, inv_eta=745, name="ICHIRO", op_field=30e6, op_temp='2K', material='bulkNb',
                    wp=wp, sigma=sigma, plot_label="ICHIRO",
                    slans_dir=fr"{parent_dir_slans}\ICHIRO",
                    abci_dir=fr"{parent_dir_abci}\ICHIRO")

    NLSF = Cavity(vrf=V, inv_eta=745, name="NLSF", op_field=30e6, op_temp='2K', material='bulkNb',
                  wp=wp, sigma=sigma, plot_label="NLSF",
                  slans_dir=fr"{parent_dir_slans}\NLSF",
                  abci_dir=fr"{parent_dir_abci}\NLSF")

    NLSF_A = Cavity(vrf=V, inv_eta=745, name="NLSF-A", op_field=30e6, op_temp='2K', material='bulkNb',
                    wp=wp, sigma=sigma, plot_label="NLSF-A",
                    slans_dir=fr"{parent_dir_slans}\NLSF-A",
                    abci_dir=fr"{parent_dir_abci}\NLSF-A")

    # NLSF_RE = Cavity(9, l_cell_mid=57.7e-3, freq=1300e6, vrf=V, R_Q=1130.0682, G=279.9786,
    #                 Epk_Eacc=2.07, Bpk_Eacc=3.83, inv_eta=745, name="NLSF-B", op_field=30e6,
    #                 op_temp='2K', material='bulkNb')

    TESLA = Cavity(vrf=V, inv_eta=745, name="TESLA", op_field=30e6, op_temp='2K', material='bulkNb',
                   wp=wp, sigma=sigma, plot_label="TESLA",
                   slans_dir=fr"{parent_dir_slans}\TESLA",
                   abci_dir=fr"{parent_dir_abci}\TESLA")

    # C3795 = Cavity(vrf=V, inv_eta=745, name="C3795", op_field=30e6, op_temp='2K', material='bulkNb',
    #                 wp=wp, sigma=sigma, plot_label="C3795",
    #                 slans_dir=fr"{parent_dir_slans}\C3795",
    #                 abci_dir=fr"{parent_dir_abci}\C3795")

    slans_dirs = [
        fr"{parent_dir_slans}\ILC-LL",
        fr"{parent_dir_slans}\ICHIRO",
        fr"{parent_dir_slans}\NLSF",
        fr"{parent_dir_slans}\NLSF-A",
        # fr"{parent_dir_slans}\NLSF-B",
        fr"{parent_dir_slans}\TESLA",
        # fr"{parent_dir_slans}\C3795_1300MHz"
    ]

    abci_dirs = [
        fr"{parent_dir_abci}\ILC-LL",
        fr"{parent_dir_abci}\ICHIRO",
        fr"{parent_dir_abci}\NLSF",
        fr"{parent_dir_abci}\NLSF-A",
        # fr"{parent_dir_abci}\NLSF-B",
        fr"{parent_dir_abci}\TESLA",
        # fr"{parent_dir_slans}\C3795_1300MHz"
    ]

    cavities = Cavities([ILC_LL, ICHIRO, NLSF, NLSF_A, TESLA], 'TESLA_and_LL_Cavities')
    # cavities = Cavities([NLSF, TESLA], 'TESLA_ERLMA_NLSF_Cavities')

    cavities.set_cavities_slans(slans_dirs)
    cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 50, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()

    # print(cavities)
    # print(c3795_tt)
    cavities.plot_cavities_contour('mid')

    cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    cavities.plot_axis_fields()
    cavities.plot_surface_fields()

    # makr summaries
    cavities.make_latex_summary_tables()
    cavities.make_excel_summary()

    # multipacting
    multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]

    to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # for tp in to_plot:
    # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    cavities.plot_dispersion()
    # plt.show()


def mucol_study_2():
    # MUCOL STUDY
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\MuCol_Study\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\MuCol_Study\SimulationData\ABCI"

    abci_dirs = [
        fr"{parent_dir_abci}\NLSF_scale_1.621797_5",
        fr"{parent_dir_abci}\NLSF_scale_1.621797_7",
        fr"{parent_dir_abci}\NLSF_scale_1.621797_9",
        fr"{parent_dir_abci}\NLSF_scale_1.3_5",
        fr"{parent_dir_abci}\NLSF_scale_1.3_7",
        fr"{parent_dir_abci}\NLSF_scale_1.3_9",
        fr"{parent_dir_abci}\NLSF_5",
        fr"{parent_dir_abci}\NLSF_7",
        fr"{parent_dir_abci}\NLSF_9"
    ]

    # RCS
    V = 20.87 * 1e9
    wp = 'MuCol RCS Stage 1'  # working point
    sigma = 'SR_13.0mm'
    machine = "MuCol"

    # ILC_LL = Cavity(9, l_cell_mid=57.7e-3, freq=1300e6, vrf=V, R_Q=1201.2516, G=284.4373,
    #                 Epk_Eacc=2.31, Bpk_Eacc=3.62, inv_eta=745, name="ILC-LL", op_field=30e6,
    #                 op_temp='2K', material='bulkNb')

    # ICHIRO = Cavity(9, l_cell_mid=57.7e-3, freq=1300e6, vrf=V, R_Q=1204.0024, G=283.9091,
    #                 Epk_Eacc=2.32, Bpk_Eacc=3.61, inv_eta=745, name="ICHIRO", op_field=30e6,
    #                 op_temp='2K', material='bulkNb')

    NLSF5_0_8GHz = Cavity(vrf=V, inv_eta=745, name="NLSF5_0.8GHz", op_field=30e6, op_temp='2K', material='bulkNb',
                          wp=wp, sigma=sigma, plot_label="NLSF$_\mathrm{5-cell, 0.8~GHz}$",
                          slans_dir=fr"{parent_dir_slans}\NLSF_scale_1.621797_5",
                          abci_dir=fr"{parent_dir_abci}\NLSF_scale_1.621797_5", Q0=1e10)
    NLSF7_0_8GHz = Cavity(vrf=V, inv_eta=745, name="NLSF7_0.8GHz", op_field=30e6, op_temp='2K', material='bulkNb',
                          wp=wp, sigma=sigma, plot_label="NLSF$_\mathrm{7-cell, 0.8~GHz}$",
                          slans_dir=fr"{parent_dir_slans}\NLSF_scale_1.621797_7",
                          abci_dir=fr"{parent_dir_abci}\NLSF_scale_1.621797_7", Q0=1e10)
    NLSF9_0_8GHz = Cavity(vrf=V, inv_eta=745, name="NLSF9_0.8GHz", op_field=30e6, op_temp='2K', material='bulkNb',
                          wp=wp, sigma=sigma, plot_label="NLSF$_\mathrm{9-cell, 0.8~GHz}$",
                          slans_dir=fr"{parent_dir_slans}\NLSF_scale_1.621797_9",
                          abci_dir=fr"{parent_dir_abci}\NLSF_scale_1.621797_9", Q0=1e10)

    NLSF5_1GHz = Cavity(vrf=V, inv_eta=745, name="NLSF5_1.0GHz", op_field=30e6, op_temp='2K', material='bulkNb',
                        wp=wp, sigma=sigma, plot_label="NLSF$_\mathrm{5-cell, 1.0~GHz}$",
                        slans_dir=fr"{parent_dir_slans}\NLSF_scale_1.3_5",
                        abci_dir=fr"{parent_dir_abci}\NLSF_scale_1.3_5", Q0=1e10)
    NLSF7_1GHz = Cavity(vrf=V, inv_eta=745, name="NLSF7_1.0GHz", op_field=30e6, op_temp='2K', material='bulkNb',
                        wp=wp, sigma=sigma, plot_label="NLSF$_\mathrm{7-cell, 1.0~GHz}$",
                        slans_dir=fr"{parent_dir_slans}\NLSF_scale_1.3_7",
                        abci_dir=fr"{parent_dir_abci}\NLSF_scale_1.3_7", Q0=1e10)
    NLSF9_1GHz = Cavity(vrf=V, inv_eta=745, name="NLSF9_1.0GHz", op_field=30e6, op_temp='2K', material='bulkNb',
                        wp=wp, sigma=sigma, plot_label="NLSF$_\mathrm{9-cell, 1.0~GHz}$",
                        slans_dir=fr"{parent_dir_slans}\NLSF_scale_1.3_9",
                        abci_dir=fr"{parent_dir_abci}\NLSF_scale_1.3_9", Q0=1e10)

    NLSF5_1_3GHz = Cavity(vrf=V, inv_eta=745, name="NLSF5_1.3GHz", op_field=30e6, op_temp='2K', material='bulkNb',
                          wp=wp, sigma=sigma, plot_label="NLSF$_\mathrm{5-cell, 1.3~GHz}$",
                          slans_dir=fr"{parent_dir_slans}\NLSF_5",
                          abci_dir=fr"{parent_dir_abci}\NLSF_5", Q0=1e10)
    NLSF7_1_3GHz = Cavity(vrf=V, inv_eta=745, name="NLSF7_1.3GHz", op_field=30e6, op_temp='2K', material='bulkNb',
                          wp=wp, sigma=sigma, plot_label="NLSF$_\mathrm{7-cell, 1.3~GHz}$",
                          slans_dir=fr"{parent_dir_slans}\NLSF_7",
                          abci_dir=fr"{parent_dir_abci}\NLSF_7", Q0=1e10)
    NLSF9_1_3GHz = Cavity(vrf=V, inv_eta=745, name="NLSF9_1.3GHz", op_field=30e6, op_temp='2K', material='bulkNb',
                          wp=wp, sigma=sigma, plot_label="NLSF$_\mathrm{9-cell, 1.3~GHz}$",
                          slans_dir=fr"{parent_dir_slans}\NLSF_9",
                          abci_dir=fr"{parent_dir_abci}\NLSF_9", Q0=1e10)

    # NLSF_A = Cavity(9, l_cell_mid=57.7e-3, freq=1300e6, vrf=V, R_Q=1172.6666, G=277.9839,
    #                 Epk_Eacc=2.07, Bpk_Eacc=3.77, inv_eta=745, name="NLSF-A", op_field=30e6,
    #                 op_temp='2K', material='bulkNb')

    # NLSF_RE = Cavity(9, l_cell_mid=57.7e-3, freq=1300e6, vrf=V, R_Q=1130.0682, G=279.9786,
    #                 Epk_Eacc=2.07, Bpk_Eacc=3.83, inv_eta=745, name="NLSF-B", op_field=30e6,
    #                 op_temp='2K', material='bulkNb')

    # TESLA = Cavity(9, l_cell_mid=57.7e-3, freq=1300e6, vrf=V, R_Q=1022.8792, G=271.3334,
    #                Epk_Eacc=1.98, Bpk_Eacc=4.17, inv_eta=745, name="TESLA", op_field=30e6,
    #                op_temp='2K', material='bulkNb')

    cavities = Cavities([NLSF5_0_8GHz, NLSF7_0_8GHz, NLSF9_0_8GHz,
                         NLSF5_1GHz, NLSF7_1GHz, NLSF9_1GHz,
                         NLSF5_1_3GHz, NLSF7_1_3GHz, NLSF9_1_3GHz], save_folder='NLSF_Cavities_freq', machine=machine,
                        particle='muon')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 50, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()

    # print(cavities)
    # print(c3795_tt)
    cavities.plot_cavities_contour('mid')

    cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    # cavities.plot_axis_fields()
    # cavities.plot_surface_fields()

    # makr summaries
    cavities.make_latex_summary_tables()
    cavities.make_excel_summary()

    # multipacting
    multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]

    to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # for tp in to_plot:
    # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()


def mucol_study_3():
    # MUCOL STUDY
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\MuCol_Study\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\MuCol_Study\SimulationData\ABCI"

    abci_dirs = [
        fr"{parent_dir_abci}\NLSF_scale_1.621797_5",
        fr"{parent_dir_abci}\NLSF_scale_1.621797_7",
        fr"{parent_dir_abci}\NLSF_scale_1.621797_9",
        fr"{parent_dir_abci}\NLSF_scale_1.3_5",
        fr"{parent_dir_abci}\NLSF_scale_1.3_7",
        fr"{parent_dir_abci}\NLSF_scale_1.3_9",
        fr"{parent_dir_abci}\NLSF_5",
        fr"{parent_dir_abci}\NLSF_7",
        fr"{parent_dir_abci}\NLSF_9"
    ]

    # RCS
    V = 20.87 * 1e9
    wp = 'MuCol RCS Stage 1'  # working point
    sigma = 'SR_13.0mm'
    machine = "MuCol"

    NLSF = Cavity(vrf=V, inv_eta=745, name="NLSF", op_field=30e6, op_temp='2K', material='bulkNb',
                  wp=wp, sigma=sigma, plot_label="NLSF",
                  slans_dir=fr"{parent_dir_slans}\NLSF",
                  abci_dir=fr"{parent_dir_abci}\NLSF", Q0=1e10)
    ERL_MA = Cavity(vrf=V, inv_eta=745, name="ERL_MA", op_field=30e6, op_temp='2K', material='bulkNb',
                    wp=wp, sigma=sigma, plot_label="ERL-MA",
                    slans_dir=fr"{parent_dir_slans}\ERL_MA",
                    abci_dir=fr"{parent_dir_abci}\ERL_MA", Q0=1e10)
    TESLA = Cavity(vrf=V, inv_eta=745, name="TESLA", op_field=30e6, op_temp='2K', material='bulkNb',
                   wp=wp, sigma=sigma, plot_label="TESLA",
                   slans_dir=fr"{parent_dir_slans}\TESLA",
                   abci_dir=fr"{parent_dir_abci}\TESLA", Q0=1e10)

    cavities = Cavities([NLSF, ERL_MA, TESLA], save_folder='NLSF_ERL_TESLA', machine=machine,
                        particle='muon')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 50, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()

    # print(cavities)
    # print(c3795_tt)
    cavities.plot_cavities_contour('mid')

    cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    # cavities.plot_axis_fields()
    # cavities.plot_surface_fields()

    # makr summaries
    cavities.make_latex_summary_tables()
    cavities.make_excel_summary()

    # multipacting
    multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]

    to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # for tp in to_plot:
    # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()


def ttbar_study():
    # 5 cell cavities comparison for ttbar
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"
    wp = 'ttbar_2023'  # working point
    sigma = 'SR_1.67mm'

    c3795_tt = Cavity(vrf=9.2e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp=wp, sigma=sigma, plot_label="C3795",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    cFCCUROS5 = Cavity(vrf=9.2e9, inv_eta=745, name="FCCUROS5", op_field=20.12e6,
                       op_temp='2K', material='bulkNb',
                       wp=wp, sigma=sigma, plot_label="FCCUROS5",
                       slans_dir=fr"{parent_dir_slans}\FCC_UROS5_n5",
                       abci_dir=fr"{parent_dir_abci}\FCC_UROS5", Q0=3e10)

    # cTESLA = Cavity(vrf=9.2e9, inv_eta=745, name="TESLA", op_field=20.12e6,
    #                 op_temp='2K', material='bulkNb',
    #                 wp=wp, sigma=sigma, plot_label="TESLA",
    #                 slans_dir=fr"{parent_dir_slans}\TESLA_800MHz",
    #                 abci_dir=fr"{parent_dir_abci}\TESLA_800MHz", Q0=3e10)

    # cavities = Cavities([c3795_tt, cFCCUROS5, cTESLA], 'Cavities_C3795_FCCUROS5_TESLA')
    cavities = Cavities([c3795_tt, cFCCUROS5], 'Cavities_C3795_FCCUROS5')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)
    #
    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar()

    # print(cavities)
    # print(c3795_tt)
    cavities.make_latex_summary_tables()
    c3795_tt.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')
    #
    # cavities.plot_ql_vs_pin()
    # cavities.plot_cryomodule_comparison()
    # cavities.plot_axis_fields()
    # cavities.plot_surface_fields()
    # cavities.make_excel_summary()

    # multipacting
    multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795", r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
                         r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    # , r"D:\Dropbox\multipacting\MPGUI21\JLab"

    to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # for tp in to_plot:
    cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')

    cavities.plot_dispersion()
    # plt.show()
    # c3795_tt.plot_ql_vs_pin()


def tt_C3795_FCCUROS5():
    # 5 cell cavities comparison for ttbar
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    c3795_tt = Cavity(vrf=9.2e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="C3795 ($\mathrm{t \\bar t}$)",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    FCCUROS5_tt = Cavity(vrf=9.2e9, inv_eta=745, name="FCCUROS5", op_field=20.12e6,
                       op_temp='2K', material='bulkNb',
                       wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="FCCUROS5 ($\mathrm{t \\bar t}$)",
                       slans_dir=fr"{parent_dir_slans}\FCC_UROS5_n5",
                       abci_dir=fr"{parent_dir_abci}\FCC_UROS5", Q0=3e10)

    cavities = Cavities([c3795_tt, FCCUROS5_tt], 'Cavities_TT_C3795_FCCUROS5')
    #
    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    # cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar(ncols=2)

    # print(cavities)
    # print(c3795_tt)
    cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')
    #
    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    # cavities.plot_axis_fields()
    # cavities.plot_surface_fields()
    # cavities.make_excel_summary()

    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795", r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    #                      r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    # # , r"D:\Dropbox\multipacting\MPGUI21\JLab"
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    #
    # cavities.plot_dispersion()
    # plt.show()
    # c3795_tt.plot_ql_vs_pin()


def htt_C3795_FCCUROS5():
    # 5 cell cavities comparison for ttbar
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    c3795_H = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='H_2023', sigma=['SR_2.5mm', 'BS_4.45mm'], plot_label="C3795 ($\mathrm{t \\bar t}$)",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    FCCUROS5_H = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="FCCUROS5", op_field=20.12e6,
                       op_temp='2K', material='bulkNb',
                       wp='H_2023', sigma=['SR_2.5mm', 'BS_4.45mm'], plot_label="FCCUROS5 ($\mathrm{t \\bar t}$)",
                       slans_dir=fr"{parent_dir_slans}\FCC_UROS5_n5",
                       abci_dir=fr"{parent_dir_abci}\FCC_UROS5", Q0=3e10)

    c3795_tt = Cavity(vrf=9.2e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="C3795 ($\mathrm{t \\bar t}$)",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    FCCUROS5_tt = Cavity(vrf=9.2e9, inv_eta=745, name="FCCUROS5", op_field=20.12e6,
                       op_temp='2K', material='bulkNb',
                       wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="FCCUROS5 ($\mathrm{t \\bar t}$)",
                       slans_dir=fr"{parent_dir_slans}\FCC_UROS5_n5",
                       abci_dir=fr"{parent_dir_abci}\FCC_UROS5", Q0=3e10)

    cavities = Cavities([c3795_H, FCCUROS5_H, c3795_tt, FCCUROS5_tt], 'Cavities_HTT_C3795_FCCUROS5')
    #
    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    # cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar(ncols=2)

    # print(cavities)
    # print(c3795_tt)
    cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')
    #
    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    # cavities.plot_axis_fields()
    # cavities.plot_surface_fields()
    # cavities.make_excel_summary()

    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795", r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    #                      r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    # # , r"D:\Dropbox\multipacting\MPGUI21\JLab"
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    #
    # cavities.plot_dispersion()
    # plt.show()
    # c3795_tt.plot_ql_vs_pin()


def htt_C3795_C3795v2():
    # 5 cell cavities comparison for ttbar
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    c3795_H = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='H_2023', sigma=['SR_2.5mm', 'BS_4.45mm'], plot_label="C3795 ($\mathrm{t \\bar t}$)",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    c3795v2_H = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='H_2023', sigma=['SR_2.5mm', 'BS_4.45mm'], plot_label="C3795v2 ($\mathrm{t \\bar t}$)",
                      slans_dir=fr"{parent_dir_slans}\C3795v2_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795v2", Q0=3e10)

    FCCUROS5_H = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="FCCUROS5", op_field=20.12e6,
                       op_temp='2K', material='bulkNb',
                       wp='H_2023', sigma=['SR_2.5mm', 'BS_4.45mm'], plot_label="FCCUROS5 ($\mathrm{t \\bar t}$)",
                       slans_dir=fr"{parent_dir_slans}\FCC_UROS5_n5",
                       abci_dir=fr"{parent_dir_abci}\FCC_UROS5", Q0=3e10)

    c3795_tt = Cavity(vrf=9.2e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="C3795 ($\mathrm{t \\bar t}$)",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    c3795v2_tt = Cavity(vrf=9.2e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="C3795v2 ($\mathrm{t \\bar t}$)",
                      slans_dir=fr"{parent_dir_slans}\C3795v2_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795v2", Q0=3e10)

    FCCUROS5_tt = Cavity(vrf=9.2e9, inv_eta=745, name="FCCUROS5", op_field=20.12e6,
                       op_temp='2K', material='bulkNb',
                       wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="FCCUROS5 ($\mathrm{t \\bar t}$)",
                       slans_dir=fr"{parent_dir_slans}\FCC_UROS5_n5",
                       abci_dir=fr"{parent_dir_abci}\FCC_UROS5", Q0=3e10)

    cavities = Cavities([c3795_H, c3795v2_H, FCCUROS5_H, c3795_tt, c3795v2_tt, FCCUROS5_tt], 'Cavities_HTT_C3795_C3795v2_FCCUROS5')
    #
    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    # cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar(ncols=2)

    # print(cavities)
    # print(c3795_tt)
    cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')
    #
    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    # cavities.plot_axis_fields()
    # cavities.plot_surface_fields()
    # cavities.make_excel_summary()

    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795", r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    #                      r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    # # , r"D:\Dropbox\multipacting\MPGUI21\JLab"
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    #
    # cavities.plot_dispersion()
    # plt.show()
    # c3795_tt.plot_ql_vs_pin()


def htt_C3795_C3795v2_C3795v3():
    # 5 cell cavities comparison for ttbar
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    c3795_H = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='H_2023', sigma=['SR_2.5mm', 'BS_4.45mm'], plot_label="C3795 (H)",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    c3795_end_H = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="C3795_end_3", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='H_2023', sigma=['SR_2.5mm', 'BS_4.45mm'], plot_label="C3795_end (H)",
                      slans_dir=fr"{parent_dir_slans}\C3795_end_3_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795_end_3", Q0=3e10)

    FCCUROS5_H = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="FCCUROS5", op_field=20.12e6,
                       op_temp='2K', material='bulkNb',
                       wp='H_2023', sigma=['SR_2.5mm', 'BS_4.45mm'], plot_label="FCCUROS5 (H)",
                       slans_dir=fr"{parent_dir_slans}\FCC_UROS5_n5",
                       abci_dir=fr"{parent_dir_abci}\FCC_UROS5", Q0=3e10)

    c3795_tt = Cavity(vrf=9.2e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="C3795 ($\mathrm{t \\bar t}$)",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    c3795_end_tt = Cavity(vrf=9.2e9, inv_eta=745, name="C3795_end_3", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="C3795_end ($\mathrm{t \\bar t}$)",
                      slans_dir=fr"{parent_dir_slans}\C3795_end_3_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795_end_3", Q0=3e10)

    FCCUROS5_tt = Cavity(vrf=9.2e9, inv_eta=745, name="FCCUROS5", op_field=20.12e6,
                       op_temp='2K', material='bulkNb',
                       wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="FCCUROS5 ($\mathrm{t \\bar t}$)",
                       slans_dir=fr"{parent_dir_slans}\FCC_UROS5_n5",
                       abci_dir=fr"{parent_dir_abci}\FCC_UROS5", Q0=3e10)

    # cavities = Cavities([c3795_H, c3795_end_H, FCCUROS5_H, c3795_tt, c3795_end_tt, FCCUROS5_tt], 'Cavities_HTT_C3795_C3795v2_FCCUROS5')
    cavities = Cavities([c3795_H, c3795_end_H], 'Cavities_HTT_C3795_C3795v2_FCCUROS5')
    #
    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    # cavities.plot_power_comparison()
    cavities.plot_compare_bar(ncols=2)
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar(ncols=2)

    # print(cavities)
    # print(c3795_tt)
    cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')
    #
    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    # cavities.plot_axis_fields()
    # cavities.plot_surface_fields()
    # cavities.make_excel_summary()

    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795", r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    #                      r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    # # , r"D:\Dropbox\multipacting\MPGUI21\JLab"
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    #
    # cavities.plot_dispersion()
    # plt.show()
    # c3795_tt.plot_ql_vs_pin()



def tt_C3795_CEPCv5():
    # 5 cell cavities comparison for ttbar
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"
    parent_dir_ngsolve = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"

    c3795_tt = Cavity(vrf=9.2e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="C3795 ($\mathrm{t \\bar t}$)",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    CEPCv5_tt = Cavity(vrf=6.1e9, inv_eta=745, name="CEPCv5", op_field=27.6e6, op_temp='2K', material='bulkNb',
                       wp='CEPC_ttbar', sigma=['SR_2.2mm', 'BS_2.9mm'], plot_label="CEPCv5 ($\mathrm{t \\bar t}$)",
                       slans_dir=fr"{parent_dir_slans}\CEPCv5_n5",
                       abci_dir=fr"{parent_dir_abci}\CEPCv5", Q0=3e10)

    cavities = Cavities([c3795_tt, CEPCv5_tt], 'Cavities_TT_C3795_CEPCv5')
    #
    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    # cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar(ncols=2)

    # print(cavities)
    # print(c3795_tt)
    cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')
    #
    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    # cavities.plot_axis_fields()
    # cavities.plot_surface_fields()
    # cavities.make_excel_summary()

    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795", r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    #                      r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    # # , r"D:\Dropbox\multipacting\MPGUI21\JLab"
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    #
    # cavities.plot_dispersion()
    # plt.show()
    # c3795_tt.plot_ql_vs_pin()



def zwhtt_C3795_FCCUROS5_booster():
    # 5 cell cavities comparison for ttbar
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    c3795_Z_b = Cavity(vrf=0.0501e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='Z_b_2024', sigma=['SR_4.0mm'], plot_label="C3795 (Z)",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    FCCUROS5_Z_b = Cavity(vrf=0.0501e9, inv_eta=745, name="FCCUROS5", op_field=20.12e6,
                       op_temp='2K', material='bulkNb',
                       wp='Z_b_2024', sigma=['SR_4.0mm'], plot_label="FCCUROS5 (Z)",
                       slans_dir=fr"{parent_dir_slans}\FCC_UROS5_n5",
                       abci_dir=fr"{parent_dir_abci}\FCC_UROS5", Q0=3e10)

    c3795_W_b = Cavity(vrf=0.0501e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='W_b_2024', sigma=['SR_4.0mm'], plot_label="C3795 (W)",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    FCCUROS5_W_b = Cavity(vrf=0.0501e9, inv_eta=745, name="FCCUROS5", op_field=20.12e6,
                       op_temp='2K', material='bulkNb',
                       wp='W_b_2024', sigma=['SR_4.0mm'], plot_label="FCCUROS5 (W)",
                       slans_dir=fr"{parent_dir_slans}\FCC_UROS5_n5",
                       abci_dir=fr"{parent_dir_abci}\FCC_UROS5", Q0=3e10)

    c3795_H_b = Cavity(vrf=0.0501e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='H_b_2024', sigma=['SR_4.0mm'], plot_label="C3795 (H)",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    FCCUROS5_H_b = Cavity(vrf=0.0501e9, inv_eta=745, name="FCCUROS5", op_field=20.12e6,
                       op_temp='2K', material='bulkNb',
                       wp='H_b_2024', sigma=['SR_4.0mm'], plot_label="FCCUROS5 (H)",
                       slans_dir=fr"{parent_dir_slans}\FCC_UROS5_n5",
                       abci_dir=fr"{parent_dir_abci}\FCC_UROS5", Q0=3e10)

    c3795_tt_b = Cavity(vrf=0.0501e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='ttbar_b_2024', sigma=['SR_4.0mm'], plot_label="C3795 ($\mathrm{t \\bar t}$)",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    FCCUROS5_tt_b = Cavity(vrf=0.0501e9, inv_eta=745, name="FCCUROS5", op_field=20.12e6,
                       op_temp='2K', material='bulkNb',
                       wp='ttbar_b_2024', sigma=['SR_4.0mm'], plot_label="FCCUROS5 ($\mathrm{t \\bar t}$)",
                       slans_dir=fr"{parent_dir_slans}\FCC_UROS5_n5",
                       abci_dir=fr"{parent_dir_abci}\FCC_UROS5", Q0=3e10)

    cavities = Cavities([c3795_Z_b, FCCUROS5_Z_b, c3795_W_b, FCCUROS5_W_b, c3795_H_b, FCCUROS5_H_b, c3795_tt_b, FCCUROS5_tt_b], 'Cavities_ZWHTT_C3795_FCCUROS5')
    #
    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    # cavities.plot_power_comparison()
    cavities.plot_compare_bar(ncols=4)
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar(ncols=4)

    # print(cavities)
    # print(c3795_tt)
    cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')
    #
    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    # cavities.plot_axis_fields()
    # cavities.plot_surface_fields()
    # cavities.make_excel_summary()

    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795", r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    #                      r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    # # , r"D:\Dropbox\multipacting\MPGUI21\JLab"
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    #
    # cavities.plot_dispersion()
    # plt.show()
    # c3795_tt.plot_ql_vs_pin()


def h_study():
    c3794_H_400 = Cavity(2, l_cell_mid=187e-3, freq=400.79e6, vrf=2.1e9, R_Q=152.8, G=198.42,
                         Epk_Eacc=2.05, Bpk_Eacc=6.39, inv_eta=219, name="C3794_400", op_field=11.87e6,
                         op_temp='4.5K', material='NbCu')

    # c3794_H_400_bulkNb = Cavity(2, l_cell_mid=187e-3, freq=400.79e6, vrf=2.1e9, R_Q=152.8, G=198.42,
    #                      Epk_Eacc=2.05, Bpk_Eacc=6.39, inv_eta=219, name="C3794_400_bulkNb", op_field=11.87e6,
    #                      op_temp='2K', material='bulkNb')

    c3794_H_800 = Cavity(2, l_cell_mid=93.5e-3, freq=801.58e6, vrf=2.1e9, R_Q=152.8, G=198.42,
                         Epk_Eacc=2.05, Bpk_Eacc=6.39, inv_eta=219, name="C3794_800", op_field=11.87e6,
                         op_temp='2K', material='bulkNb')

    c3794_H_800_2 = Cavity(2, l_cell_mid=93.5e-3, freq=801.58e6, vrf=2.1e9, R_Q=152.8, G=198.42,
                           Epk_Eacc=2.05, Bpk_Eacc=6.39, inv_eta=219, name="C3794_800$_\mathrm{23.74 MV/m}$",
                           op_field=2 * 11.87e6,
                           op_temp='2K', material='bulkNb')

    # c3795_H_400 = Cavity(5, l_cell_mid=93.5e-3, freq=801.58e6, vrf=2.1e9 / 2, R_Q=448.12, G=261.63,
    #                      Epk_Eacc=2.43, Bpk_Eacc=4.88, inv_eta=745, name="C3795", op_field=24.72e6,
    #                      op_temp='4.5K', material='NbCu')

    c3795_H_800 = Cavity(5, l_cell_mid=93.5e-3, freq=801.58e6, vrf=2.1e9, R_Q=448.12, G=261.63,
                         Epk_Eacc=2.43, Bpk_Eacc=4.88, inv_eta=745, name="C3795_800", op_field=24.72e6,
                         op_temp='2K', material='bulkNb')

    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\ABCI"

    # 2 and 5 cell, 400MHz, 800MHz cavities comparison for H
    # wp = 'H'  # working point
    # sigma = 'SR_2.5mm'
    slans_dirs = [fr"{parent_dir_slans}\3794_400", fr"{parent_dir_slans}\3794_800",
                  fr"{parent_dir_slans}\C3795_800"]
    abci_dirs = [fr"{parent_dir_abci}\3794_400", fr"{parent_dir_abci}\3794_800",
                 fr"{parent_dir_abci}\C3795_800"]
    cavities = Cavities([c3794_H_400, c3794_H_800, c3795_H_800], 'Cavities_C3794_400_800_C3795_800')

    cavities.set_cavities_slans(slans_dirs)
    cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()

    # print(cavities)
    # print(c3795_tt)
    cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('end')

    cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    cavities.plot_axis_fields()
    cavities.plot_surface_fields()
    cavities.make_excel_summary()

    # multipacting
    multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]

    to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # for tp in to_plot:
    # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    cavities.plot_dispersion()
    # plt.show()


def w_c3794_cepcv2_study():
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_ngsolve = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    C3794 = Cavity(vrf=1 * 1e9, inv_eta=745, name="C3794", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                   wp='W_2023', sigma='SR_3.55mm', plot_label="C3794 (400 MHz)",
                   slans_dir=fr"{parent_dir_slans}\C3794_n2",
                   abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    CEPCv2 = Cavity(vrf=4.3e6, inv_eta=745, name="CEPCv2", op_field=9.3e6, op_temp='2K', material='bulkNb',
                    wp='CEPC_W', sigma='SR_3.4mm', plot_label="CEPCv2 (650 MHz)",
                    slans_dir=fr"{parent_dir_ngsolve}\CEPCv2_n2",
                    abci_dir=fr"{parent_dir_abci}\CEPCv2", Q0=1e10, cell_type='flat top')

    slans_dirs = [fr"{parent_dir_slans}\C3794_n2", fr"{parent_dir_ngsolve}\CEPCv2_n2"]
    abci_dirs = [fr"{parent_dir_abci}\3794", fr"{parent_dir_abci}\CEPCv2"]

    cavities = Cavities([C3794, CEPCv2], 'Cavities_W_C3794_CEPCv2')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar()
    #
    # # print(cavities)
    # # print(c3795_tt)
    # cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')

    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    cavities.plot_axis_fields()
    cavities.plot_surface_fields()
    cavities.make_excel_summary()
    #
    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()
    # plt.show()


def whtt_c3794_cepcv2_study():
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_ngsolve = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    C3794_W = Cavity(vrf=1e9, inv_eta=219, name="C3794", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                     wp='W_2023', sigma=['SR_3.55mm', 'BS_7.02mm'], plot_label="C3794 (W)",
                     slans_dir=fr"{parent_dir_slans}\C3794_n2",
                     abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    CEPCv2_W = Cavity(vrf=0.7e9, inv_eta=745, name="CEPCv2", op_field=9.1e6, op_temp='2K', material='bulkNb',
                      wp='CEPC_W', sigma=['SR_2.5mm', 'BS_4.9mm'], plot_label="CEPCv2 (W)",
                      slans_dir=fr"{parent_dir_ngsolve}\CEPCv2_n2",
                      abci_dir=fr"{parent_dir_abci}\CEPCv2", Q0=3e10, cell_type='flat top')

    C3794_H = Cavity(vrf=2.1 * 1e9, inv_eta=219, name="C3794", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                     wp='H_2023', sigma=['SR_2.5mm', 'BS_4.45mm'], plot_label="C3794 (H)",
                     slans_dir=fr"{parent_dir_slans}\C3794_n2",
                     abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    CEPCv2_H = Cavity(vrf=2.2e9, inv_eta=745, name="CEPCv2", op_field=14.2e6, op_temp='2K', material='bulkNb',
                      wp='CEPC_H', sigma=['SR_2.3mm', 'BS_4.1mm'], plot_label="CEPCv2 (H)",
                      slans_dir=fr"{parent_dir_ngsolve}\CEPCv2_n2",
                      abci_dir=fr"{parent_dir_abci}\CEPCv2", Q0=3e10, cell_type='flat top')

    C3794_tt = Cavity(vrf=2.1 * 1e9, inv_eta=219, name="C3794", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                      wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="C3794 ($\mathrm{t \\bar t}$)",
                      slans_dir=fr"{parent_dir_slans}\C3794_n2",
                      abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    CEPCv2_tt = Cavity(vrf=3.9e9, inv_eta=745, name="CEPCv2", op_field=25.2e6, op_temp='2K', material='bulkNb',
                       wp='CEPC_ttbar', sigma=['SR_2.2mm', 'BS_2.9mm'], plot_label="CEPCv2 ($\mathrm{t \\bar t}$)",
                       slans_dir=fr"{parent_dir_ngsolve}\CEPCv2_n2",
                       abci_dir=fr"{parent_dir_abci}\CEPCv2", Q0=3e10, cell_type='flat top')

    cavities = Cavities([C3794_W, CEPCv2_W, C3794_H, CEPCv2_H, C3794_tt, CEPCv2_tt], 'Cavities_WHTT_C3794_CEPCv2')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    # cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar()
    #
    # # print(cavities)
    # # print(c3795_tt)
    cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')

    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    cavities.plot_axis_fields()
    cavities.plot_surface_fields()
    cavities.make_excel_summary()
    #
    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()
    # plt.show()


def w_c3794_baseline_study():
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_ngsolve = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    C3794 = Cavity(vrf=1 * 1e9, inv_eta=745, name="C3794", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                   wp='W_2023', sigma=['SR_3.55mm', 'BS_7.02mm'], plot_label="C3794 (2-cell)",
                   slans_dir=fr"{parent_dir_slans}\C3794_n2",
                   abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    baseline = Cavity(vrf=1e9, inv_eta=745, name="FCCUROS", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                      wp='W_2023', sigma=['SR_3.55mm', 'BS_7.02mm'], plot_label="FCCUROS (4-cell)",
                      slans_dir=fr"{parent_dir_ngsolve}\Baseline_n4",
                      abci_dir=fr"{parent_dir_abci}\Baseline", Q0=1e10)

    cavities = Cavities([C3794, baseline], 'Cavities_W_C3794_Baseline')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    # cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar()
    #
    # # print(cavities)
    # # print(c3795_tt)
    # cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')

    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    cavities.plot_axis_fields()
    cavities.plot_surface_fields()
    cavities.make_excel_summary()
    #
    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()
    # plt.show()


def whtt_c3794_baseline_study():
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_ngsolve = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    C3794_W = Cavity(vrf=1 * 1e9, inv_eta=745, name="C3794", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                     wp='W_2023', sigma=['SR_3.55mm', 'BS_7.02mm'], plot_label="C3794 (W)",
                     slans_dir=fr"{parent_dir_slans}\C3794_n2",
                     abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    FCCUROS_W = Cavity(vrf=1e9, inv_eta=745, name="FCCUROS", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                       wp='W_2023', sigma=['SR_3.55mm', 'BS_7.02mm'], plot_label="FCCUROS (W)",
                       slans_dir=fr"{parent_dir_ngsolve}\Baseline_n4",
                       abci_dir=fr"{parent_dir_abci}\Baseline", Q0=1e10)

    C3794_H = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="C3794", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                     wp='H_2023', sigma=['SR_2.5mm', 'BS_4.45mm'], plot_label="C3794 (H)",
                     slans_dir=fr"{parent_dir_slans}\C3794_n2",
                     abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    FCCUROS_H = Cavity(vrf=2.1e9, inv_eta=745, name="FCCUROS", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                       wp='H_2023', sigma=['SR_2.5mm', 'BS_4.45mm'], plot_label="FCCUROS (H)",
                       slans_dir=fr"{parent_dir_ngsolve}\Baseline_n4",
                       abci_dir=fr"{parent_dir_abci}\Baseline", Q0=1e10)

    C3794_tt = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="C3794", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                      wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="C3794 ($\mathrm{t \\bar t}$)",
                      slans_dir=fr"{parent_dir_slans}\C3794_n2",
                      abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    FCCUROS_tt = Cavity(vrf=2.1e9, inv_eta=745, name="FCCUROS", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                        wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="FCCUROS ($\mathrm{t \\bar t}$)",
                        slans_dir=fr"{parent_dir_ngsolve}\Baseline_n4",
                        abci_dir=fr"{parent_dir_abci}\Baseline", Q0=1e10)

    cavities = Cavities([C3794_W, FCCUROS_W, C3794_H, FCCUROS_H, C3794_tt, FCCUROS_tt], 'Cavities_WHTT_C3794_FCCUROS')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    # cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar()
    #
    # # print(cavities)
    # # print(c3795_tt)
    cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')

    cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    cavities.plot_axis_fields()
    cavities.plot_surface_fields()
    cavities.make_excel_summary()
    #
    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()
    # plt.show()


def w_c3794_baseline_cepcv2_study():
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_ngsolve = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    C3794 = Cavity(vrf=1 * 1e9, inv_eta=745, name="C3794", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                   wp='W_2023', sigma='SR_3.55mm', plot_label="C3794 (2-cell)",
                   slans_dir=fr"{parent_dir_slans}\C3794_n2",
                   abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    baseline = Cavity(vrf=4.3e6, inv_eta=745, name="Baseline", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                      wp='W_2023', sigma='SR_3.55mm', plot_label="Baseline (4-cell)",
                      slans_dir=fr"{parent_dir_ngsolve}\Baseline_n4",
                      abci_dir=fr"{parent_dir_abci}\Baseline", Q0=1e10)

    CEPCv2 = Cavity(vrf=4.3e6, inv_eta=745, name="CEPCv2", op_field=9.3e6, op_temp='2K', material='bulkNb',
                    wp='CEPC_W', sigma='SR_3.4mm', plot_label="CEPCv2 (2_cell)",
                    slans_dir=fr"{parent_dir_ngsolve}\CEPCv2_n2",
                    abci_dir=fr"{parent_dir_abci}\CEPCv2", Q0=1e10, cell_type='flat top')

    slans_dirs = [fr"{parent_dir_slans}\C3794_n2", fr"{parent_dir_ngsolve}\CEPCv2_n2"]
    abci_dirs = [fr"{parent_dir_abci}\3794", fr"{parent_dir_abci}\CEPCv2"]

    cavities = Cavities([C3794, baseline, CEPCv2], 'Cavities_W_C3794_Baseline_CEPCv2')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar()
    #
    # # print(cavities)
    # # print(c3795_tt)
    # cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')

    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    cavities.plot_axis_fields()
    cavities.plot_surface_fields()
    cavities.make_excel_summary()
    #
    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()
    # plt.show()


def w_c3794_c3794v2_study():
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_ngsolve = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    C3794 = Cavity(vrf=1 * 1e9, inv_eta=745, name="C3794", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                   wp='W_2023', sigma='SR_3.55mm', plot_label="C3794 (400 MHz)",
                   slans_dir=fr"{parent_dir_slans}\C3794_n2",
                   abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    C3794v2 = Cavity(vrf=1 * 1e9, inv_eta=745, name="C3794v2", op_field=10.61e6, op_temp='2K', material='bulkNb',
                     wp='W_2023', sigma='SR_3.55mm', plot_label="C3794v2 (800 MHz)",
                     slans_dir=fr"{parent_dir_slans}\C3794_s0.5_n2",
                     abci_dir=fr"{parent_dir_abci}\C3794_s0.5", Q0=1e10)

    # slans_dirs = [fr"{parent_dir_slans}\C3794_n2", fr"{parent_dir_ngsolve}\CEPCv2_n2"]
    # abci_dirs = [fr"{parent_dir_abci}\3794", fr"{parent_dir_abci}\CEPCv2"]

    # cavities = Cavities([C3794, baseline, C3794v2, CEPCv2], 'Cavities_C3794_Baseline_C3794v2_CEPCv2')
    cavities = Cavities([C3794, C3794v2], 'Cavities_W_C3794_C3794v2')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    # cavities.plot_compare_fm_bar()
    # cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar()
    #
    # # print(cavities)
    # # print(c3795_tt)
    # cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')

    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    cavities.plot_axis_fields()
    cavities.plot_surface_fields()
    cavities.make_excel_summary()
    #
    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()
    # plt.show()


def wh_c3794_c3794v2_study():
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_ngsolve = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    C3794_W = Cavity(vrf=1 * 1e9, inv_eta=745, name="C3794 (W)", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                     wp='W_2023', sigma='SR_3.55mm', plot_label="C3794_W (400 MHz)",
                     slans_dir=fr"{parent_dir_slans}\C3794_n2",
                     abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    C3794v2_W = Cavity(vrf=1 * 1e9, inv_eta=745, name="C3794v2 (W)", op_field=21.22e6, op_temp='2K', material='bulkNb',
                       wp='W_2023', sigma='SR_3.55mm', plot_label="C3794v2_W (800 MHz)",
                       slans_dir=fr"{parent_dir_slans}\C3794_s0.5_n2",
                       abci_dir=fr"{parent_dir_abci}\C3794_s0.5", Q0=3e10)

    C3794_H = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="C3794_H", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                     wp='H_2023', sigma='SR_2.5mm', plot_label="C3794_H (400 MHz)",
                     slans_dir=fr"{parent_dir_slans}\C3794_n2",
                     abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    C3794v2_H = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="C3794v2_H", op_field=21.22e6, op_temp='2K', material='bulkNb',
                       wp='H_2023', sigma='SR_2.5mm', plot_label="C3794v2_H (800 MHz)",
                       slans_dir=fr"{parent_dir_slans}\C3794_s0.5_n2",
                       abci_dir=fr"{parent_dir_abci}\C3794_s0.5", Q0=3e10)

    # slans_dirs = [fr"{parent_dir_slans}\C3794_n2", fr"{parent_dir_ngsolve}\CEPCv2_n2"]
    # abci_dirs = [fr"{parent_dir_abci}\3794", fr"{parent_dir_abci}\CEPCv2"]

    # cavities = Cavities([C3794, baseline, C3794v2, CEPCv2], 'Cavities_C3794_Baseline_C3794v2_CEPCv2')
    cavities = Cavities([C3794_W, C3794v2_W, C3794_H, C3794v2_H], 'Cavities_WH_C3794_C3794v2')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar()
    #
    # # print(cavities)
    # # print(c3795_tt)
    # cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')

    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    cavities.plot_axis_fields()
    cavities.plot_surface_fields()
    cavities.make_excel_summary()
    #
    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()
    # plt.show()


def whtt_c3794_c3794v2_study():
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_ngsolve = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    C3794_W = Cavity(vrf=1 * 1e9, inv_eta=219, name="C3794 (W)", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                     wp='W_2023', sigma=['SR_3.55mm', 'BS_7.02mm'], plot_label="C3794 (W)",
                     slans_dir=fr"{parent_dir_slans}\C3794_n2",
                     abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    C3794v2_W = Cavity(vrf=1 * 1e9, inv_eta=745, name="C3794v2 (W)", op_field=21.22e6, op_temp='2K', material='bulkNb',
                       wp='W_2023', sigma=['SR_3.55mm', 'BS_7.02mm'], plot_label="C3794v2 (W)",
                       slans_dir=fr"{parent_dir_slans}\C3794_s0.5_n2",
                       abci_dir=fr"{parent_dir_abci}\C3794_s0.5", Q0=3e10)

    C3794_H = Cavity(vrf=2.1 * 1e9, inv_eta=219, name="C3794_H", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                     wp='H_2023', sigma=['SR_2.5mm', 'BS_4.45mm'], plot_label="C3794 (H)",
                     slans_dir=fr"{parent_dir_slans}\C3794_n2",
                     abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    C3794v2_H = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="C3794v2_H", op_field=21.22e6, op_temp='2K', material='bulkNb',
                       wp='H_2023', sigma=['SR_2.5mm', 'BS_4.45mm'], plot_label="C3794v2 (H)",
                       slans_dir=fr"{parent_dir_slans}\C3794_s0.5_n2",
                       abci_dir=fr"{parent_dir_abci}\C3794_s0.5", Q0=3e10)

    C3794_ttbar = Cavity(vrf=2.1e9, inv_eta=219, name="C3794_ttabr", op_field=10.61e6,
                         op_temp='4.5K', material='NbCu',
                         wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'], plot_label="C3794 ($\mathrm{t \\bar t}$)",
                         slans_dir=fr"{parent_dir_slans}\C3794_n2",
                         abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    C3794v2_ttbar = Cavity(vrf=2.1e9, inv_eta=745, name="C3794v2_ttbar", op_field=21.22e6,
                           op_temp='2K', material='bulkNb',
                           wp='ttbar_2023', sigma=['SR_1.67mm', 'BS_2.54mm'],
                           plot_label="C3794v2 ($\mathrm{t \\bar t}$)",
                           slans_dir=fr"{parent_dir_slans}\C3794_s0.5_n2",
                           abci_dir=fr"{parent_dir_abci}\C3794_s0.5", Q0=3e10)

    # slans_dirs = [fr"{parent_dir_slans}\C3794_n2", fr"{parent_dir_ngsolve}\CEPCv2_n2"]
    # abci_dirs = [fr"{parent_dir_abci}\3794", fr"{parent_dir_abci}\CEPCv2"]

    # cavities = Cavities([C3794, baseline, C3794v2, CEPCv2], 'Cavities_C3794_Baseline_C3794v2_CEPCv2')
    cavities = Cavities([C3794_W, C3794v2_W, C3794_H, C3794v2_H, C3794_ttbar, C3794v2_ttbar],
                        'Cavities_WHttbar_C3794_C3794v2')

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    # cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar(ncols=3)
    #
    # # print(cavities)
    # # print(c3795_tt)
    cavities.make_latex_summary_tables()
    cavities.make_excel_summary()
    # cavities.plot_cavities_contour('mid')
    # cavities.plot_cavities_contour('end')

    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    # cavities.plot_axis_fields()
    # cavities.plot_surface_fields()
    #
    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()
    # plt.show()


def w_c3794_baseline_c3794v2_cepcv2_study():
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_ngsolve = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    C3794 = Cavity(vrf=1 * 1e9, inv_eta=745, name="C3794", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                   wp='W_2023', sigma='SR_3.55mm', plot_label="C3794 (2-cell,400 MHz)",
                   slans_dir=fr"{parent_dir_slans}\C3794_n2",
                   abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    baseline = Cavity(vrf=4.3e6, inv_eta=745, name="Baseline", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                      wp='W_2023', sigma='SR_3.55mm', plot_label="Baseline (4-cell,400 MHz)",
                      slans_dir=fr"{parent_dir_ngsolve}\Baseline_n4",
                      abci_dir=fr"{parent_dir_abci}\Baseline", Q0=1e10)

    C3794v2 = Cavity(vrf=1 * 1e9, inv_eta=745, name="C3794v2", op_field=10.61e6, op_temp='2K', material='bulkNb',
                     wp='W_2023', sigma='SR_3.55mm', plot_label="C3794v2 (2-cell,800 MHz)",
                     slans_dir=fr"{parent_dir_slans}\C3794_s0.5_n2",
                     abci_dir=fr"{parent_dir_abci}\C3794_s0.5", Q0=1e10)

    CEPCv2 = Cavity(vrf=1e9, inv_eta=745, name="CEPCv2", op_field=9.3e6, op_temp='2K', material='bulkNb',
                    wp='CEPC_W', sigma='SR_3.4mm', plot_label="CEPCv2 (2-cell,650 MHz)",
                    slans_dir=fr"{parent_dir_ngsolve}\CEPCv2_n2",
                    abci_dir=fr"{parent_dir_abci}\CEPCv2", Q0=1e10, cell_type='flat top')

    # slans_dirs = [fr"{parent_dir_slans}\C3794_n2", fr"{parent_dir_ngsolve}\CEPCv2_n2"]
    # abci_dirs = [fr"{parent_dir_abci}\3794", fr"{parent_dir_abci}\CEPCv2"]

    # cavities = Cavities([C3794, baseline, C3794v2, CEPCv2], 'Cavities_C3794_Baseline_C3794v2_CEPCv2')
    cavities = Cavities([C3794, baseline, C3794v2, CEPCv2], 'Cavities_W_C3794_Baseline_C3794v2_CEPCv2')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar()
    #
    # # print(cavities)
    # # print(c3795_tt)
    # cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')

    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    cavities.plot_axis_fields()
    cavities.plot_surface_fields()
    cavities.make_excel_summary()
    #
    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()
    # plt.show()


def h_c3794_c3794v2_study():
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    parent_dir_ngsolve = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"

    C3794 = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="C3794", op_field=20.12e6, op_temp='4.5K', material='NbCu',
                   wp='H_2023', sigma='SR_2.5mm', plot_label="C3794 (400 MHz)",
                   slans_dir=fr"{parent_dir_slans}\C3794_n2",
                   abci_dir=fr"{parent_dir_abci}\C3794", Q0=3e10)

    C3794v2 = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="C3794v2", op_field=20.12e6, op_temp='2K', material='bulkNb',
                     wp='H_2023', sigma='SR_2.5mm', plot_label="C3794 (800 MHz)",
                     slans_dir=fr"{parent_dir_slans}\C3794_s0.5_n2",
                     abci_dir=fr"{parent_dir_abci}\C3794_s0.5", Q0=3e10)

    cavities = Cavities([C3794, C3794v2], 'Cavities_H_C3794_C3794v2')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar()
    #
    # # print(cavities)
    # # print(c3795_tt)
    # cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')

    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    cavities.plot_axis_fields()
    cavities.plot_surface_fields()
    cavities.make_excel_summary()
    #
    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()
    # plt.show()


def h_c3794_cepcv2_study():
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    parent_dir_ngsolve = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"

    C3794 = Cavity(vrf=2.1 * 1e9, inv_eta=745, name="C3794", op_field=20.12e6, op_temp='4.5K', material='NbCu',
                   wp='H_2023', sigma='SR_2.5mm', plot_label="C3794 (400 MHz)",
                   slans_dir=fr"{parent_dir_slans}\C3794_n2",
                   abci_dir=fr"{parent_dir_abci}\C3794", Q0=3e10)

    CEPCv2 = Cavity(vrf=2.1e9, inv_eta=745, name="CEPCv2", op_field=20.12e6, op_temp='2K', material='bulkNb',
                    wp='CEPC_H', sigma='SR_2.9mm', plot_label="CEPCv2",
                    slans_dir=fr"{parent_dir_ngsolve}\CEPCv2_n2",
                    abci_dir=fr"{parent_dir_abci}\CEPCv2", Q0=3e10, cell_type='flat top')

    cavities = Cavities([C3794, CEPCv2], 'Cavities_H_C3794_CEPCv2')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    cavities.plot_compare_all_bar()
    #
    # # print(cavities)
    # # print(c3795_tt)
    # cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')

    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    cavities.plot_axis_fields()
    cavities.plot_surface_fields()
    cavities.make_excel_summary()
    #
    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()
    # plt.show()


def h_c3794_C3795_cepcv2_study():
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    parent_dir_ngsolve = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"

    c3795_tt = Cavity(vrf=9.2e9, inv_eta=745, name="C3795", op_field=20.12e6,
                      op_temp='2K', material='bulkNb',
                      wp='ttbar_2023', sigma='SR_1.67mm', plot_label="C3795",
                      slans_dir=fr"{parent_dir_slans}\C3795_n5",
                      abci_dir=fr"{parent_dir_abci}\C3795", Q0=3e10)

    CEPCv2 = Cavity(vrf=2.1e9, inv_eta=745, name="CEPCv2", op_field=9.3e6, op_temp='2K', material='bulkNb',
                    wp='CEPC_W', sigma='SR_3.4mm', plot_label="CEPCv2",
                    slans_dir=fr"{parent_dir_ngsolve}\CEPCv2_n2",
                    abci_dir=fr"{parent_dir_abci}\CEPCv2", Q0=1e10, cell_type='flat top')

    slans_dirs = [fr"{parent_dir_slans}\C3794_n2", fr"{parent_dir_ngsolve}\CEPCv2_n2"]
    abci_dirs = [fr"{parent_dir_abci}\3794", fr"{parent_dir_abci}\CEPCv2"]

    cavities = Cavities([c3795_tt, CEPCv2], 'Cavities_H_C3794_400_800_C3795_800')

    # cavities.set_cavities_slans(slans_dirs)
    # cavities.set_cavities_abci(abci_dirs)

    E_acc = np.linspace(0.5, 30, 100) * 1e6  # V/m
    cavities.compare_power(E_acc=E_acc)
    cavities.plot_power_comparison()
    cavities.plot_compare_bar()
    cavities.plot_compare_fm_bar()
    cavities.plot_compare_hom_bar()
    #
    # # print(cavities)
    # # print(c3795_tt)
    # cavities.make_latex_summary_tables()
    cavities.plot_cavities_contour('mid')
    cavities.plot_cavities_contour('end')

    # cavities.plot_ql_vs_pin()
    cavities.plot_cryomodule_comparison()
    cavities.plot_axis_fields()
    cavities.plot_surface_fields()
    cavities.make_excel_summary()
    #
    # # multipacting
    # multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3795"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]
    #
    # to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # # for tp in to_plot:
    # # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()
    # plt.show()


def c3794_study():
    # MUCOL STUDY
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    abci_dirs = [
        fr"{parent_dir_abci}\C3794"
    ]

    # RCS
    V = 1 * 1e9  # <- V
    wp = 'W_2023'  # working point
    sigma = 'SR_3.55mm'
    machine = "FCC-ee"

    C3794 = Cavity(vrf=V, inv_eta=745, name="C3794", op_field=10.61e6, op_temp='4.5K', material='NbCu',
                   wp=wp, sigma=sigma, plot_label="C3794",
                   slans_dir=fr"{parent_dir_slans}\C3794_n2",
                   abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    cavities = Cavities([C3794], save_folder='C3794', machine=machine, particle='electron')

    E_acc = np.arange(5, 51, 1) * 1e6  # V/m
    cavities.set_cavities_field(E_acc)
    # cavities.compare_power(E_acc=E_acc)
    # cavities.plot_power_comparison()
    # cavities.plot_compare_bar()
    # cavities.plot_compare_fm_bar()
    # cavities.plot_compare_hom_bar()

    # print(cavities)
    # print(c3795_tt)
    # cavities.plot_cavities_contour('mid')

    # cavities.plot_ql_vs_pin(np.arange(5, 16, 0.5)*1e6)

    Eacc_ = np.arange(4, 20 + 2, 1) * 1e6
    Vrf = np.arange(0.2, 2 + 0.2, 0.1) * 1e9
    history = [np.array([9.6, 11.91, 10.61]) * 1e6, np.array([0.75, 0.44, 1.05]) * 1e9]
    cavities.plot_ql_vs_pin(Eacc_, Vrfin=Vrf, mesh=True, hist=history)
    plt.tight_layout()
    plt.show()
    # cavities.plot_ql_vs_pin(np.arange(5, 16, 5)*1e6)
    # cavities.plot_cryomodule_comparison()
    # cavities.plot_axis_fields()
    # cavities.plot_surface_fields()

    # # makr summaries
    # cavities.make_latex_summary_tables()
    # cavities.make_excel_summary()

    # multipacting
    multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3794"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]

    to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # for tp in to_plot:
    # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()


def c3794_study_Z():
    # MUCOL STUDY
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    abci_dirs = [
        fr"{parent_dir_abci}\C3794"
    ]

    # RCS
    V = 1 * 1e9  # <- V
    wp = 'Z_2023'  # working point
    sigma = 'SR_4.32mm'
    machine = "FCC-ee"

    C3794 = Cavity(vrf=V, inv_eta=745, name="C3794", op_field=3.81e6, op_temp='4.5K', material='NbCu',
                   wp=wp, sigma=sigma, plot_label="C3794",
                   slans_dir=fr"{parent_dir_slans}\C3794_n2",
                   abci_dir=fr"{parent_dir_abci}\C3794", Q0=1e10)

    cavities = Cavities([C3794], save_folder='C3794', machine=machine, particle='electron')

    E_acc = np.arange(1, 51, 1) * 1e6  # V/m
    cavities.set_cavities_field(E_acc)
    # cavities.compare_power(E_acc=E_acc)
    # cavities.plot_power_comparison()
    # cavities.plot_compare_bar()
    # cavities.plot_compare_fm_bar()
    # cavities.plot_compare_hom_bar()

    # print(cavities)
    # print(c3795_tt)
    # cavities.plot_cavities_contour('mid')

    # cavities.plot_ql_vs_pin(np.arange(5, 16, 0.5)*1e6)

    Eacc_ = np.arange(1, 10 + 2, 1) * 1e6
    Vrf = np.arange(0.01, 2 + 0.01, 0.01) * 1e9
    # history = [np.array([9.6, 11.91, 10.61])*1e6, np.array([0.75, 0.44, 1.05])*1e9]
    cavities.plot_ql_vs_pin(Eacc_, Vrfin=Vrf, mesh=True)
    plt.tight_layout()
    plt.show()
    # cavities.plot_ql_vs_pin(np.arange(5, 16, 5)*1e6)
    # cavities.plot_cryomodule_comparison()
    # cavities.plot_axis_fields()
    # cavities.plot_surface_fields()

    # # makr summaries
    # cavities.make_latex_summary_tables()
    # cavities.make_excel_summary()

    # multipacting
    multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3794"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]

    to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # for tp in to_plot:
    # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()


def fcc_uros_1_1_study_Z():
    # MUCOL STUDY
    parent_dir_slans = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"
    parent_dir_abci = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\ABCI"

    abci_dirs = [
        fr"{parent_dir_abci}\FCC_UROS1_1"
    ]

    # RCS
    V = 1 * 1e9  # <- V
    wp = 'Z_2023'  # working point
    sigma = 'SR_4.32mm'
    machine = "FCC-ee"

    C3794 = Cavity(vrf=V, inv_eta=745, name="FCC_UROS1_1", op_field=3.81e6, op_temp='4.5K', material='NbCu',
                   wp=wp, sigma=sigma, plot_label="FCC_UROS1_1",
                   slans_dir=fr"{parent_dir_slans}\FCC_UROS1_1_n1",
                   abci_dir=fr"{parent_dir_abci}\FCC_UROS1_1", Q0=2.7e9)

    cavities = Cavities([C3794], save_folder='FCC_UROS1_1', machine=machine, particle='electron')

    E_acc = np.arange(1, 51, 1) * 1e6  # V/m
    cavities.set_cavities_field(E_acc)
    # cavities.compare_power(E_acc=E_acc)
    # cavities.plot_power_comparison()
    # cavities.plot_compare_bar()
    # cavities.plot_compare_fm_bar()
    # cavities.plot_compare_hom_bar()

    # print(cavities)
    # print(c3795_tt)
    # cavities.plot_cavities_contour('mid')

    # cavities.plot_ql_vs_pin(np.arange(5, 16, 0.5)*1e6)

    Eacc_ = np.arange(1, 10 + 2, 1) * 1e6
    Vrf = np.arange(0.01, 0.1 + 0.01, 0.01) * 1e9
    # history = [np.array([9.6, 11.91, 10.61])*1e6, np.array([0.75, 0.44, 1.05])*1e9]
    cavities.plot_ql_vs_pin(Eacc_, Vrfin=Vrf, mesh=True)
    plt.tight_layout()
    plt.show()
    # cavities.plot_ql_vs_pin(np.arange(5, 16, 5)*1e6)
    # cavities.plot_cryomodule_comparison()
    # cavities.plot_axis_fields()
    # cavities.plot_surface_fields()

    # # makr summaries
    # cavities.make_latex_summary_tables()
    # cavities.make_excel_summary()

    # multipacting
    multipact_folders = [r"D:\Dropbox\multipacting\MPGUI21\C3794"]  # , r"D:\Dropbox\multipacting\MPGUI21\FCCUROS5",
    # r"D:\Dropbox\multipacting\MPGUI21\TESLA"]

    to_plot = ['counter function', 'final impact energy', 'enhanced counter function']
    # for tp in to_plot:
    # cavities.plot_multipac_triplot(multipact_folders, 'enhanced counter function')
    # cavities.plot_dispersion()


def write_thresholds(filepath, fm=1, filter=False, order=None):
    with open(filepath, 'r') as f:
        op = json.load(f)

    df = pd.DataFrame.from_dict(op)
    if filter:
        filt = [col for col in df.columns if filter in col]
        df = df[filt]
    df = df.T
    print(df)

    # calculate thresholds
    df['long. threshold'] = (2 * df['E [GeV]'] * 1e9 * df['nu_s []'])/(df['N_c []'] * df['I0 [mA]'] * 1e-3 * df['alpha_p [1e-5]'] * 1e-5 * df['tau_z [ms]'] * 1e-3 * fm) * 1e-3
    df['trans. threshold'] = (2 * df['E [GeV]']) * 1e9 / (df['N_c []'] * df['I0 [mA]'] * 1e-3 * df['beta_xy [m]'] * df['tau_xy [ms]'] * 1e-3 * df['f_rev [kHz]'] * 1e3) * 1e-3
    if order:
        print(df.T[order].to_latex(float_format="{:.2f}".format))
    else:
        print(df.T.to_latex(float_format="{:.2f}".format))



if __name__ == '__main__':
    # mucol_study_2()
    # mucol_study_3()

    # c3794_study()
    # c3794_study_Z()
    # fcc_uros_1_1_study_Z()

    # w_c3794_cepcv2_study()
    # whtt_c3794_cepcv2_study()
    # w_c3794_baseline_study()
    # whtt_c3794_baseline_study()
    # w_c3794_c3794v2_study()
    # wh_c3794_c3794v2_study()
    # whtt_c3794_c3794v2_study()
    # w_c3794_baseline_cepcv2_study()
    # w_c3794_baseline_c3794v2_cepcv2_study()
    # h_c3794_cepcv2_study()
    # h_c3794_c3794v2_study()
    # h_c3794_C3795_cepcv2_study()
    # h_study()
    # ttbar_study()
    # htt_C3795_FCCUROS5()
    # tt_C3795_CEPCv5()
    # htt_C3795_C3795v2()
    htt_C3795_C3795v2_C3795v3()
    # tt_C3795_FCCUROS5()
    # zwhtt_C3795_FCCUROS5_booster()

    # write_thresholds(r'D:\Dropbox\CavityDesignHub\PhD_Thesis\OperatingPoints\cepc.json', fm=650e6)
    # print()
    # write_thresholds(r'D:\Dropbox\CavityDesignHub\PhD_Thesis\OperatingPoints\fcc.json', fm=400.79e6, filter='2023')

    # order = ['Z_b_2024', 'Z_b_2024_FB', 'W_b_2024', 'W_b_2024_FB', 'H_b_2024', 'H_b_2024_FB', 'ttbar_b_2024', 'ttbar_b_2024_FB']
    # write_thresholds(r'D:\Dropbox\CavityDesignHub\PhD_Thesis\OperatingPoints\fcc.json', fm=801.58e6, filter='b_2024', order=None)
    # write_thresholds(r'D:\Dropbox\CavityDesignHub\PhD_Thesis\OperatingPoints\fcc.json', fm=400.79e6, filter='b_2023')
