import os
import shutil

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from icecream import ic
import scipy.io as spio
import scipy.interpolate as sci


def plot_settings(mode):
    # plt.style.use(['science', 'no-latex'])
    if mode == 'presentation':
        # plt.rcParams["figure.figsize"] = (11, 6)
        #
        # mpl.rcParams['lines.markersize'] = 10
        #
        # # axes
        # mpl.rcParams['axes.labelsize'] = 20
        # mpl.rcParams['axes.titlesize'] = 25
        #
        # mpl.rcParams['xtick.labelsize'] = 15
        # mpl.rcParams['ytick.labelsize'] = 15
        #

        # plot settings
        mpl.rcParams['xtick.labelsize'] = 20
        mpl.rcParams['ytick.labelsize'] = 20

        mpl.rcParams['axes.labelsize'] = 20
        mpl.rcParams['axes.titlesize'] = 20
        mpl.rcParams['legend.fontsize'] = 20

        mpl.rcParams['figure.figsize'] = [7, 4.5]
        mpl.rcParams['figure.dpi'] = 100

    if mode.lower() == 'square':
        # plot settings
        mpl.rcParams['xtick.labelsize'] = 20
        mpl.rcParams['ytick.labelsize'] = 20

        mpl.rcParams['axes.labelsize'] = 20
        mpl.rcParams['axes.titlesize'] = 20
        mpl.rcParams['legend.fontsize'] = 20

        mpl.rcParams['figure.figsize'] = [6, 6]
        mpl.rcParams['figure.dpi'] = 100


def plot_electron_evolution_spark3d():
    fig, ax = plt.subplots()
    # x_label = "Time [ns]"
    x_label = "Peak electric field [MV/m]"
    y_label = "Multipactor order"
    # ax.set_title("$\mathbf{Spark3D}$")
    Power = np.array([[1, 0], [2, 0], [5, 0], [7, 0], [9, 0], [10, 0], [11, 1.01],
                      [12, 1], [13, 1], [14, 1], [15, 1], [16, 1], [17, 0.998], [18, 0.996],
                      [19, 0.994], [20, 0.992], [25, 0.978], [30, 0.951], [31, 0], [32, 0],
                      [33, 0], [34, 0], [35, 0], [36, 0], [37, 0], [38, 0], [39, 0], [40, 0],
                      [45, 0], [50, 0], [55, 0], [60, 0]])
    Power_Cu = np.array([[1, 0], [2, 0], [5, 0], [7, 0], [9, 1.01], [10, 1], [11, 0.996],
                      [12, 0.995], [13, 0.991], [14, 0.986], [15, 0.984], [16, 0.979], [17, 0.972], [18, 0.965],
                      [19, 0.958], [20, 0.948], [25, 0.891], [30, 0.822], [31, 0.809], [32, 0.795],
                      [33, 0.781], [34, 0.769], [35, 0.776], [36, 0.739], [37, 0.72], [38, 0.711], [39, 0.697], [40, 0.684],
                      [45, 0.641], [50, 0.603], [55, 0.578], [60, 0.556], [70, 0.524], [80, 0.499], [100, 0.447],
                      [150, 0.367], [2000, 0]])

    Power_HOM_Coupler = np.array([[25000, 0], [30000, 0], [35000, 0], [40000, 0.841], [50000, 0.783],
                                  [60000, 0.728], [70000, 0.67], [80000, 0.624], [90000, 0], [100000, 0]])

    Power_HOM_Coupler_Cu = np.array([[10000, 0], [15000, 1.02], [16000, 1.01], [18000, 1], [20000, 1],
                                     [30000, 0.998], [35000, 0.995], [40000, 0.991], [45000, 0.988],
                                     [50000, 0.984], [60000, 0.975], [70000, 0.965], [80000, 0.955], [90000, 0.937], [100000, 0]])

    Power = np.array([[0.01, 0], [0.05, 0], [0.1, 0], [0.15, 0], [0.2, 0], [0.25, 0], [0.3, 1.01],
                      [0.5, 1], [1, 1], [1.5, 1], [2, 1], [2.5, 1], [3, 0.998], [3.5, 0.996],
                      [4, 0.994], [4.5, 0.992], [5, 0.978], [5.5, 0.951], [6, 0], [6.5, 0],
                      [7, 0], [7.5, 0], [8, 0], [10, 0], [11, 0], [12, 0], [13, 0], [14, 0],
                      [15, 0], [16, 0], [20, 0], [22, 0], [25, 0], [30, 0], [35, 0], [40, 0], [50, 0],
                      [60, 0], [70, 0], [80, 0]])

    Power_list = [Power, Power_Cu]
    Power_list_coupler = [Power_HOM_Coupler, Power_HOM_Coupler_Cu]
    labels = ['Nb', "Cu"]
    E = 1.7e6/2.05  # peak electric field of the imported field for the analysis

    # for i, p in enumerate(Power_list):
    #     ax.plot(np.sqrt(2) * np.sqrt(p[:, 0]) * E * 1e-6, p[:, 1], label=labels[i])

    # data = pd.read_csv("D:\CST Studio\Multipacting\Results\C3794_1-60W.csv", delimiter=' ', header=None)
    data = pd.read_csv("D:\CST Studio\Multipacting\E_C3794_HC\R1.csv", delimiter=' ', header=None)
    for i, p in enumerate(Power):

        # if p[1] > 0:
        #     ax.plot(data[0][data[i+1]>0]*1e9, data[i+1][data[i+1]>0], label=f"{np.sqrt(2)*np.sqrt(p[0])*E*1e-6:.2f} MV/m", lw=2)
        # else:
        #     ax.plot(data[0][data[i+1]>0]*1e9, data[i+1][data[i+1]>0], label=f"{np.sqrt(2)*np.sqrt(p[0])*E*1e-6:.2f} MV/m", lw=2, ls='--', alpha=0.5)

        ax.plot(data[0][data[i + 1] > 0] * 1e9, data[i + 1][data[i + 1] > 0], label=f"{np.sqrt(2) * np.sqrt(p[0]) * E * 1e-6:.2f}", lw=2)

    # ax.set_xlabel(x_label)
    # ax.set_ylabel(y_label)
    ax.set_ylabel("Electron Evolution")
    ax.set_xlabel("Time [ns]")
    # ax.set_yscale("log")
    ax.legend(title="$E_\mathrm{acc} ~\mathrm{[MV/m]}$", ncol=4, loc="upper right")

    ax.grid(True, which="both", ls=":")
    ax.minorticks_on()
    # ax.set_xlim(0, ax.get_xlim()[-1])
    ax.set_xlim(0, 1000)
    ax.set_ylim(0, ax.get_ylim()[-1])
    # plt.legend()
    plt.tight_layout()


def plot_sey():
    fig, ax = plt.subplots()
    x_label = "Incident Energy [eV]"
    y_label = "SEY"
    data = pd.read_csv("D:\CST Studio\Multipacting\SEY\secy1.txt", sep='\s+', header=None)
    data2 = pd.read_csv("D:\CST Studio\Multipacting\SEY\secy2.txt", sep='\s+', header=None)
    sey_list = [data, data2]
    label = ["Nb", "Cu"]

    for i, sey in enumerate(sey_list):
        ax.plot(sey[0], sey[1], lw=2, label=label[i])

    ax.axhline(1, ls='--', c='r')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(0, 1000)
    ax.set_ylim(0, 2.2)
    plt.legend()
    ax.grid(True, which="both", ls=":")
    ax.minorticks_on()
    plt.tight_layout()


def plot_multipac_triplot(Eacc, Epk_Eacc):
    mpl.rcParams['figure.figsize'] = [6, 7.5]
    # load_output_data
    # files
    fnames = ["Ccounter.mat", "Acounter.mat", "Atcounter.mat", "Efcounter.mat", "param",
              "geodata.n", "secy1", "counter_flevels.mat", "counter_initials.mat"]
    data = {}
    files_folder = "D:\Dropbox\multipacting\MPGUI21"
    for f in fnames:
        if ".mat" in f:
            data[f] = spio.loadmat(fr"{files_folder}\\{f}")
        else:
            data[f] = pd.read_csv(fr"{files_folder}\\{f}", sep='\s+', header=None)

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

            fig, axs = plt.subplots(3)
            axs[0].plot(Efl / 1e6, C / n)
            axs[0].set_ylabel("$c_" + "{" + f"{N}" + "}/ c_0 $")
            axs[0].set_xlabel('Peak electric field  [MV/m]')
            axs[0].set_title(r'$\mathbf{MultiPac 2.1~~~~~~~~~~~Counter function~~~~~~~~~~~~~~~~~}$')
            axs[0].set_xlim(np.amin(Efl) / 1e6, np.amax(Efl) / 1e6)
            axs[0].set_ylim(0, np.max([0.1, axs[0].get_ylim()[1]]))

            axs[1].semilogy(Efl / 1e6, Efq)

            # axs[1].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [secy1[e1, 0], secy1[e1, 0]], '-r')
            e0 = sci.interp1d(secy1[0:e1, 1], secy1[0:e1, 0])(0)
            axs[1].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [e0, e0], '-r')
            axs[1].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [secy1[e2, 0], secy1[e2, 0]], '-r')
            axs[1].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [secy1[e3, 0], secy1[e3, 0]], '--r')

            axs[1].set_ylabel("$Ef_" + "{" + f"{N}" + "}$")
            axs[1].set_xlabel('Peak electric field  [MV/m]')
            axs[1].set_title('$\mathbf{Final~Impact~Energy~in~eV}$')
            axs[1].set_xlim(np.min(Efl) / 1e6, np.max(Efl) / 1e6)
            axs[1].set_ylim(0, axs[1].get_ylim()[1])

            axs[2].semilogy(Efl / 1e6, (A + 1) / n)
            axs[2].set_xlabel('Voltage   [MV]')
            axs[2].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [1, 1], '-r')
            axs[2].set_xlim(np.min(Efl) / 1e6, np.max(Efl) / 1e6)
            axs[2].set_ylim(np.min((A + 1) / n), axs[2].get_ylim()[1])
            axs[2].set_ylabel("$e_" + "{" + f"{N}" + "}" + "/ c_0$")
            axs[2].set_xlabel('Peak electric field   [MV/m]')
            axs[2].set_title('$\mathbf{Enhanced~counter~function}$')

            # plot peak operating field
            axs[0].axvline(Eacc*Epk_Eacc, c='k', ls='--', lw=3)
            axs[1].axvline(Eacc*Epk_Eacc, c='k', ls='--', lw=3)
            axs[2].axvline(Eacc*Epk_Eacc, c='k', ls='--', lw=3)

            axs[0].grid(True, which="both", ls=":")
            axs[1].grid(True, which="both", ls=":")
            axs[2].grid(True, which="both", ls=":")
            axs[0].minorticks_on()
            axs[1].minorticks_on()
            axs[2].minorticks_on()


def plot_trajectory():
    files_folder = "D:\Dropbox\multipacting\MPGUI21"
    fieldparams = pd.read_csv(fr"{files_folder}\\fieldparam", sep='\s+', header=None).to_numpy()
    geodata = pd.read_csv(fr"{files_folder}\\geodata.n", sep='\s+', header=None).to_numpy()
    param = pd.read_csv(fr"{files_folder}\\param", sep='\s+', header=None).to_numpy()
    elecpath = pd.read_csv(fr"{files_folder}\\elecpath", sep='\s+', header=None).to_numpy()

    gtype = fieldparams[0]

    ng = len(geodata[:, 0])
    bo = geodata[4:ng, 0: 1].T
    wr = []
    wz = []

    eps = np.spacing(1.0)
    par = param
    
    n = np.shape(elecpath)[0]
    if n == 1:
        pat = []
        ic(['No electron emission. Please, define a new initial point.'])
    else:
        pat = elecpath[1:n, [0, 2, 3, 5, 6, 7]]
        N = par[5]
        hit = np.where(pat[:, 5] != 0)
        hit = hit[np.array(0, len(hit))]
        speed = np.sqrt(pat[hit, 2]**2 + pat[hit, 3]**2)
        c = 2.9979e8
        M = 9.1093879e-31
        q = 1.6021773e-19
        energy = np.multiply(np.divide(1., np.sqrt(1.0 - np.divide(speed**2, c**2))) - 1), M*c**2./q
        avegene = np.mean(energy)
        finaene = energy[len(energy)]
        maxiene = np.max(energy)

        fig, axs = plt.subplots(3)

        axs[0].plot(bo[1, :], bo[0, :], '-b')
        axs[0].plot(pat[:, 1], pat[:, 0], '-r')
        axs[0].set_title(f'MultiPac 2.1       Electron Trajectory,   N = {N},     ')
        dt = abs(np.min(pat[:, 4]) - np.max(pat[:, 4])) *par[1]
        axs[0].set_xlabel(f'z-axis [m],  flight time {dt} periods')
        axs[0].set_ylabel('r-axis [m]')

        axs[1].plot(bo[1,:], bo[0,:], '-b')
        axs[1].plot(pat[:, 1], pat[:, 0], '-r', pat(hit, 2), pat(hit, 1), 'ro')

        min1 = 0.9 * np.min(pat[:, 0])-eps
        max1 = 1.1 * np.max(pat[:, 0])+eps
        if np.min(pat[: ,1]) < 0:
            min2 = 1.1 * np.min(pat[:, 1])-eps
        else:
            min2 = 0.9 * np.min(pat[:, 2])-eps
        
    if np.max(pat[: ,1]) < 0:
        max2 = 0.9 * np.max(pat[:, 1]) + eps
    else:
        max2 = 1.1 * max(pat[:, 1])+eps
        # axis([min2, max2, min1, max1])
        axs[1].set_xlabel('z-axis [m]')
        axs[1].set_ylabel('r-axis [m]')

        axs[2].plot(pat[:, 4]*par[0], pat[:, 0], 'r', pat[hit, 4] * par[0], pat[hit, 0], 'ro'),
        # midax = axis
        # midax(3) = max(min(bo(1,:)), min(pat(:, 1) - 1e-3))
        # midax(4) = min(max(bo(1,:)), max(pat(:, 1) + 1e-3))
        # axis(midax)
        axs[2].set_xlabel(f"time in [1/f], average energy {avegene} final energy {finaene}")
        axs[2].set_ylabel('r-axis [m]')


def sensitivity():
    fig, ax = plt.subplots()
    x = [1, 2, 3, 4, 5, 6, 7]
    s = [-0.987984504, 0.21315, -0.21626, 0.086534, 0.45319, 0.275492, -1.25225]
    labels = ["$A$", "$B$", "$a$", "$b$", "$R_\mathrm{i}$", "$L$", "$R_\mathrm{eq}$" ]
    ax.bar(x, s, align='center', width=1, color=['#1f77b4' if v < 0 else '#ff7f0e' for v in s])
    ax.set_xticks(x, labels)
    ax.set_ylabel(r"$\mathrm{\Delta}f/\mathrm{\Delta}p_i$ [MHz/mm]")
    ax.set_title(r"Sensi½tivity of $f_{\mathrm{FM}}$ [MHz] to geometric variables $p_i$ [mm]")
    for bars in ax.containers:
        ax.bar_label(bars, fontsize=18, label_type='center')


def plot_cavity():
    # data = pd.read_csv(r"D:\Dropbox\CEMCodesHub\C1092V\PostprocessingData\Data\3794_geom.txt", sep='\s+', header=None)
    # data1 = pd.read_csv(r"D:\Dropbox\CEMCodesHub\C1092V\PostprocessingData\Data\2183_geom.txt", sep='\s+', header=None)
    # data3 = pd.read_csv(r"D:\Dropbox\CEMCodesHub\C1092V\PostprocessingData\Data\650_geom.txt", sep='\s+', header=None)
    # data4 = pd.read_csv(r"D:\Dropbox\CEMCodesHub\C1092V\PostprocessingData\Data\770_geom.txt", sep='\s+', header=None)

    # ll = [650, 770, 2183, 3345, 3794, 4123, 4250, 4618]
    ll = ['C40866_geom', 'C3794_800MHz_geom', "G6_C170_M_geom"]
    laf = ['C40866', 'C3794_800MHz', "G6_C170_M"]
    for i, x in enumerate(ll):
        data = pd.read_csv(fr"D:\Dropbox\CEMCodesHub\C800MHz\PostprocessingData\Data\{x}.txt", sep='\s+', header=None)
        plt.rcParams["figure.figsize"] = (5, 5)
        plt.plot(data[1]*1000, data[0]*1000, lw=5, label=laf[i], ls='--')
        # plt.plot(data1[1]*1e3, data1[0]*1e3, lw=6, label="C2183", ls='--')
        # plt.plot(data3[1]*1e3, data3[0]*1e3, lw=6, label="C650", ls='--')
        # plt.plot(data4[1]*1e3, data4[0]*1e3, lw=6, label="C770", ls='--')
        # plt.plot(data1[1]*1e3, data1[0]*1e3, lw=3, label="$\mathrm{FCC_{UROS1.0}}$")
        # plt.plot(data2[1]*1e3, data2[0]*1e3, lw=3, label="$\mathrm{FCC_{UROS1.1}}$")
        plt.legend(loc='lower left')

        x_label = "z [mm]"
        y_label = "r [mm]"
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.xlim(-0.1, 95)
        plt.ylim(-0.1, 200)
        # plt.savefig(fr'D:\Dropbox\Quick presentation files\{x}.png', format='png', transparent=True)
        # plt.cla()


# plt.style.use('fivethirtyeight')
plot_settings("presentation")

# plot_electron_evolution_spark3d()
# plot_sey()
plot_multipac_triplot(25, 2.38)
# plot_trajectory()
# sensitivity()
# plot_cavity()


plt.tight_layout()
plt.show()

# folder = fr"D:\Dropbox\CEMCodesHub\Cavity800\SimulationData"
# folders = os.listdir(folder)
# for d in folders:
#     # delete SLANS_EXE folder
#     if os.path.exists(fr'{folder}/{d}/SLANS_exe'):
#         shutil.rmtree(fr'{folder}/{d}/SLANS_exe')
#         print(fr"Removed from {d}")