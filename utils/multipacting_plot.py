import os
import random
import shutil
from math import atan2, degrees
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy
from icecream import ic
import scipy.io as spio
import scipy.interpolate as sci
import matplotlib.patheffects as patheffects
import vtk

from scipy.signal import butter, filtfilt

# from labellines import labelLines


# Label line with line2D label data
from matplotlib import cm
from matplotlib.colorbar import ColorbarBase


def labelLine(line, x, label=None, align=True, **kwargs):
    ax = line.axes
    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if (x < xdata[0]) or (x > xdata[-1]):
        print('x label location is outside data range!')
        return

    # Find corresponding y co-ordinate and angle of the line
    ip = 1
    for i in range(len(xdata)):
        if x < xdata[i]:
            ip = i
            break

    y = ydata[ip - 1] + (ydata[ip] - ydata[ip - 1]) * (x - xdata[ip - 1]) / (xdata[ip] - xdata[ip - 1])

    if not label:
        label = line.get_label()

    if align:
        # Compute the slope
        if ip > 20:
            dx = xdata[ip + 4] - xdata[ip - 20]
            dy = ydata[ip + 4] - ydata[ip - 20]
        else:
            dx = xdata[ip] - xdata[ip - 1]
            dy = ydata[ip] - ydata[ip - 1]
        ang = degrees(atan2(dy, dx))

        # Transform to screen co-ordinates
        pt = np.array([x, y]).reshape((1, 2))
        trans_angle = ax.transData.transform_angles(np.array((ang,)), pt)[0]

    else:
        trans_angle = 0

    # Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = line.get_color()

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'center'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
        kwargs['backgroundcolor'] = ax.get_facecolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5

    txt = ax.text(x, y, label, rotation=trans_angle, **kwargs)
    txt.set_size(8)
    txt.set_path_effects([patheffects.withStroke(linewidth=0.4, foreground='k')])
    txt.set_bbox(dict(facecolor='white', alpha=0.0, edgecolor='white'))


def labelLines(lines, align=True, xvals=None, **kwargs):
    ax = lines[0].axes
    labLines = []
    labels = []

    # Take only the lines which have labels other than the default ones
    for line in lines:
        label = line.get_label()
        if "_line" not in label:
            labLines.append(line)
            labels.append(label)

    if xvals is None:
        xmin, xmax = ax.get_xlim()
        xvals = np.linspace(xmin, xmax, len(labLines) + 2)[1:-1]

    for line, x, label in zip(labLines, xvals, labels):
        labelLine(line, x, label, align, **kwargs)


def plot_settings():
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12

    mpl.rcParams['axes.labelsize'] = 12
    mpl.rcParams['axes.titlesize'] = 12
    mpl.rcParams['legend.fontsize'] = 12
    mpl.rcParams['legend.title_fontsize'] = 12

    mpl.rcParams['figure.dpi'] = 100

    # Set the desired colormap
    plt.rcParams['axes.prop_cycle'] = plt.cycler('color', plt.cm.Set2.colors)


def get_power_from_directories(root_folder, keyword):
    power = []
    for dirpath, dirnames, filenames in os.walk(root_folder):
        for filename in filenames:
            if keyword in filename:
                power.append(float(filename.split('_')[-2]))
    return np.sort(power)


def get_files_from_directories(root_folder, keyword):
    filename_list = []
    for dirpath, dirnames, filenames in os.walk(root_folder):
        for filename in filenames:
            if keyword in filename:
                filename_list.append(os.path.join(dirpath, filename))
    return np.sort(filename_list)


def plot_electron_evolution_spark3d(filename, power_folder, Eacc, sep='\s+', layout=None):
    if layout is None:
        layout = [[0, 1]]
    kwargs = {"width": 0.1, "alpha": 1, "ec": 'k', "lw": 0.5}

    figsize = (12, 3 * len(layout))

    fig, axs = plt.subplot_mosaic(layout, figsize=figsize)

    df = pd.read_csv(filename, sep=sep, header=None)
    t = df.iloc[:, 0]
    evolution = df.iloc[:, 1:]

    growth_factor = evolution.T[evolution.index[-1]] / evolution.T[evolution.index[0]]

    # calculate Eacc from input power
    # get input powers from folder
    power = get_power_from_directories(power_folder, "Surface3D_stats")

    axs[0].semilogy(t * 1e9, evolution, label=[f"{np.sqrt(2) * np.sqrt(p) * Eacc * 1e-6:.2f}" for p in power])
    axs[0].set_xlabel('Time [ns]')
    axs[0].set_ylabel('Electron evolution')
    # axs[0].legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
    #               mode="expand", borderaxespad=0, ncol=6)

    Eacc_list = np.sqrt(2) * np.sqrt(power) * Eacc * 1e-6
    # axs[1].semilogy(Eacc_list, growth_factor, label='Growth factor', marker='o', mec='k')
    axs[1].bar(Eacc_list, growth_factor, label='Growth factor', **kwargs)

    axs[1].set_ylabel("Growth factor")
    axs[1].set_xlabel("$E_\mathrm{acc} ~\mathrm{[MV/m]}$")
    axs[1].legend(ncol=6, loc="upper right")

    axs[0].minorticks_on()
    axs[0].set_xlim(0, 50)
    axs[1].minorticks_on()
    axs[1].set_xlim(1, 12)

    kwargs = {}

    # get label locations
    xvals_index = []
    for column_name in evolution.columns:
        column = evolution[column_name]
        xvals_index.append((column != 0).cumsum().idxmax())

    xvals_index = [xval - random.randint(5, 10) for xval in xvals_index if xval > 10]
    labelLines(axs[0].get_lines(), xvals=t[xvals_index] * 1e9, zorder=2.5, align=False, **kwargs)

    plt.tight_layout()
    plt.show()


def plot_sey():
    fig, ax = plt.subplots()
    x_label = "Incident Energy [eV]"
    y_label = "SEY"
    # data = pd.read_csv("D:\CST Studio\Multipacting\SEY\secy1.txt", sep='\s+', header=None)
    # data2 = pd.read_csv("D:\CST Studio\Multipacting\SEY\secy2.txt", sep='\s+', header=None)
    data = pd.read_csv("D:\Dropbox\multipacting\MPGUI21\secy1", sep='\s+', header=None)
    sey_list = [data]  # , data2]
    label = ["Nb", "Cu"]

    for i, sey in enumerate(sey_list):
        ax.plot(sey[0], sey[1], lw=1.5, label=label[i])

    ax.axhline(1, ls='--', c='r')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(0, 1000)
    # ax.set_ylim(0, 2.2)
    plt.legend()
    # ax.grid(True, which="both", ls=":")
    ax.minorticks_on()
    plt.tight_layout()
    plt.show()


def plot_multipac_triplot(Eacc_list, Epk_Eacc_list, folders, labels, kind='triplot', layout=None, min_max=None):
    use_input_min_max = False
    if layout is None:
        layout = [[0], [1], [2]]

    # mpl.rcParams['figure.figsize'] = [6, 7]

    # if kind == 'triplot':
    figsize = (10, 3 * len(layout))
    fig, axs = plt.subplot_mosaic(layout, figsize=figsize)
    kwargs = {"width": 0.1, "alpha": 1, "ec": 'k', "lw": 0.5}

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

        A = data["Acounter.mat"]["A"][:, 0]
        At = data["Atcounter.mat"]["At"]
        C = data["Ccounter.mat"]["C"][:, 0]
        Ef = data["Efcounter.mat"]["Ef"][:, 0]
        flevel = data["counter_flevels.mat"]["flevel"]
        initials = data["counter_initials.mat"]["initials"]

        secy1 = data["secy1"].to_numpy()
        Pow = flevel
        n = len(initials[:, 0]) / 2  # number of initials in the bright set
        N = int(data["param"].to_numpy()[4])  # number of impacts
        U = flevel
        Efl = flevel[:, 0]
        q = 1.6021773e-19
        Efq = Ef / q

        e1 = np.min(np.where(secy1[:, 1] >= 1))  # lower threshold
        e2 = np.max(np.where(secy1[:, 1] >= 1))  # upper threshold
        val, e3 = np.max(secy1[:, 1]), np.argmax(secy1[:, 1])  # maximum secondary yield

        if kind == 'counter_function' or kind == 'triplot':
            if kind == 'counter_function':
                ax = axs[0]
            else:
                ax = axs[0]
            # ax.semilogy(Efl / (Epk_Eacc * 1e6), C / n, lw=1.5, label=label, marker='o', mec='k')
            ax.bar(Efl / (Epk_Eacc * 1e6), C / n, label=label, **kwargs)
            # ax.set_yscale('log')
            ax.set_ylabel("$c_" + "{" + f"{N}" + "}/ c_0 $")
            ax.set_xlabel(r'$E_\mathrm{acc}$ [MV/m]')
            # ax.set_title(r'$\mathbf{MultiPac 2.1~~~~~Counter function~~~~}$')

            if min_max:
                if len(min_max) == 2:
                    use_input_min_max = True

            if use_input_min_max:
                ax.set_xlim(min_max[0], min_max[1])
            else:
                ax.set_xlim(np.amin(Efl) / (Epk_Eacc * 1e6), np.amax(Efl) / (Epk_Eacc * 1e6))

            ax.set_ylim(0, np.max([0.1, axs[0].get_ylim()[1]]))

            # plot peak operating field
            ax.axvline(Eacc, c='k', ls='--', lw=1.5)
            ax.text(Eacc, 0.5, f"{np.round(Eacc, 2)} MV/m",
                    size=12, rotation=90,
                    transform=ax.get_xaxis_transform(), ha='right', va='center')

            ax.minorticks_on()

        if kind == 'final_impact_energy' or kind == 'triplot':
            if kind == 'final_impact_energy':
                ax = axs[0]
            else:
                ax = axs[1]
            # ax.plot(Efl / (Epk_Eacc * 1e6), Efq, lw=1.5, label=label, marker='o', mec='k')
            ax.bar(Efl / (Epk_Eacc * 1e6), Efq, label=label, **kwargs)
            ax.set_yscale('log')

            # axs[1].plot([np.min(Efl) / 1e6, np.max(Efl) / 1e6], [secy1[e1, 0], secy1[e1, 0]], '-r')
            e0 = sci.interp1d(secy1[0:e1 + 1, 1], secy1[0:e1 + 1, 0])(1)
            ax.plot([np.min(Efl) / (Epk_Eacc * 1e6), np.max(Efl) / (Epk_Eacc * 1e6)], [e0, e0], '-r')
            ax.plot([np.min(Efl) / (Epk_Eacc * 1e6), np.max(Efl) / (Epk_Eacc * 1e6)],
                    [secy1[e2, 0], secy1[e2, 0]], '-r')
            ax.plot([np.min(Efl) / (Epk_Eacc * 1e6), np.max(Efl) / (Epk_Eacc * 1e6)],
                    [secy1[e3, 0], secy1[e3, 0]], '--r')

            ax.set_ylabel("$Ef_" + "{" + f"{N}" + "} [eV]$")
            ax.set_xlabel(r'$E_\mathrm{acc}$ [MV/m]')
            # ax.set_title('$\mathbf{Final~Impact~Energy~in~eV}$')

            if min_max:
                if len(min_max) == 2:
                    use_input_min_max = True

            if use_input_min_max:
                ax.set_xlim(min_max[0], min_max[1])
            else:
                ax.set_xlim(np.min(Efl) / (Epk_Eacc * 1e6), np.max(Efl) / (Epk_Eacc * 1e6))
            ax.set_ylim(0, axs[0].get_ylim()[1])

            ax.axvline(Eacc, c='k', ls='--', lw=1.5)
            ax.text(Eacc, 0.5, f"{np.round(Eacc, 2)} MV/m",
                    size=12, rotation=90,
                    transform=ax.get_xaxis_transform(), ha='right', va='center')

            ax.minorticks_on()

        if kind == 'enhanced_counter_function' or kind == 'triplot':
            if kind == 'enhanced_counter_function':
                ax = axs[0]
            else:
                ax = axs[2]

            # ax.plot(Efl / (Epk_Eacc * 1e6), (A + 1) / n, lw=1.5, label=label, marker='o', mec='k')
            ax.bar(Efl / (Epk_Eacc * 1e6), (A + 1) / n, label=label, **kwargs)
            ax.set_yscale('log')
            ax.set_xlabel('$V$ [MV]')
            ax.plot([np.min(Efl) / (Epk_Eacc * 1e6), np.max(Efl) / (Epk_Eacc * 1e6)], [1, 1], '-r')

            if min_max:
                if len(min_max) == 2:
                    use_input_min_max = True

            if use_input_min_max:
                ax.set_xlim(min_max[0], min_max[1])
            else:
                ax.set_xlim(np.min(Efl) / (Epk_Eacc * 1e6), np.max(Efl) / (Epk_Eacc * 1e6))
            ax.set_ylim(np.min((A + 1) / n), ax.get_ylim()[1])
            ax.set_ylabel("$e_" + "{" + f"{N}" + "}" + "/ c_0$")
            ax.set_xlabel(r'$E_\mathrm{acc}$ [MV/m]')
            # ax.set_title('$\mathbf{Enhanced~counter~function}$')

            ax.axvline(Eacc, c='k', ls='--', lw=1.5)
            ax.text(Eacc, 0.5, f"{np.round(Eacc, 2)} MV/m",
                    size=12, rotation=90,
                    transform=ax.get_xaxis_transform(), ha='right', va='center')

            ax.minorticks_on()

    axs[0].legend(loc='upper right')

    fig.tight_layout()
    fig.savefig(fr'D:\Dropbox\Quick presentation files\Multipacting_{label}_multipac_{kind}.png')
    plt.show()
    plt.close(fig)


def plot_trajectory(files_folder, loc='center'):
    fieldparams = pd.read_csv(fr"{files_folder}\\fieldparam", sep='\s+', header=None).to_numpy()
    geodata = pd.read_csv(fr"{files_folder}\\geodata.n", sep='\s+', header=None).to_numpy()
    param = pd.read_csv(fr"{files_folder}\\param", sep='\s+', header=None).to_numpy()
    elecpath = pd.read_csv(fr"{files_folder}\\elecpath", sep='\s+', header=None).to_numpy()

    gtype = fieldparams[0]

    ng = len(geodata[:, 0])
    bo = geodata[3:ng, 0:2].T
    wr = []
    wz = []

    # convert boundary and electron path to mm
    bo = bo*1e3

    eps = np.spacing(1.0)
    par = param

    n = np.shape(elecpath)[0]
    if n == 1:
        pat = []
        ic(['No electron emission. Please, define a new initial point.'])
    else:
        pat = elecpath[1:n, [0, 2, 3, 5, 6, 7]]

        N = par[4]
        hit = np.array(np.where(pat[:, 5] != 0))
        hit = hit[:, np.arange(1, len(hit[0]), 2)]
        speed = np.sqrt(pat[hit, 2] ** 2 + pat[hit, 3] ** 2)
        c = 2.9979e8
        M = 9.1093879e-31
        q = 1.6021773e-19
        energy = (1 / np.sqrt(1.0 - (speed ** 2 / c ** 2)) - 1) * M * c ** 2. / q
        avegene = np.mean(energy)
        finaene = energy[:, len(energy)]
        maxiene = np.max(energy)
        dt = abs(np.min(pat[:, 4]) - np.max(pat[:, 4])) * par[0]

        fig, axs = plt.subplot_mosaic([[0, 0, 2, 2, 2]], figsize=(10, 4))
        pat = pat * 1e3  # convert to mm for plotting

        axs[0].plot(bo[1, :], bo[0, :], lw=2)
        axs[0].plot(pat[:, 1], pat[:, 0], '-r', lw=1)
        axs[0].margins(x=0.01, y=0.01)
        # fig.suptitle(f'MultiPac 2.1       Electron Trajectory,   N = {N},     ')
        axs[0].set_xlabel(f'z-axis [mm],  \nFlight time {dt[0]:.2f} periods')
        axs[0].set_ylabel('r-axis [mm]')

        # axs[0].set_xlim(right=2*(axs[0].get_xlim()[1] - axs[0].get_xlim()[0]))

        if loc == 'center':
            axin = axs[0].inset_axes([0.275, 0.1, 0.45, 0.45])
        elif loc == 'left':
            axin = axs[0].inset_axes([0.05, 0.5, 0.45, 0.45])
        else:
            axin = axs[0].inset_axes([0.5, 0.5, 0.45, 0.45])

        axin.set_xticklabels([])
        axin.set_yticklabels([])

        #
        axin.plot(bo[1, :], bo[0, :], lw=2)
        axin.plot(pat[:, 1], pat[:, 0], '-r', pat[hit, 1], pat[hit, 0], 'ro', lw=1)

        min1 = np.min(pat[:, 0]) - eps*1e3
        max1 = np.max(pat[:, 0]) + eps*1e3

        if np.min(pat[:, 1]) < 0:
            min2 = 1.1 * np.min(pat[:, 1]) - eps*1e3
        else:
            min2 = 0.9 * np.min(pat[:, 1]) - eps*1e3

        if np.max(pat[:, 1]) < 0:
            max2 = 0.9 * np.max(pat[:, 1]) + eps*1e3
        else:
            max2 = 1.1 * max(pat[:, 1]) + eps*1e3

        axin.set_xlim([min2, max2])
        axin.set_ylim([min1-0.05, max1+0.05])
        axs[0].indicate_inset_zoom(axin)
        # axin.set_xlim([max(min(bo[1, :]), min(pat[:, 1] - 1)), min(max(bo[1, :]), max(pat[:, 1] + 1))])

        # axin.set_xlabel('z-axis [mm]')
        # axin.set_ylabel('r-axis [mm]')

        axs[2].plot(pat[:, 4] * par[0]*1e-3, pat[:, 0], 'r', pat[hit, 4] * par[0]*1e-3, pat[hit, 0], 'ro', markerfacecolor="None")

        # axs[2].set_ylim(top=max(pat[:, 0]))

        axs[2].set_xlabel(f"Time in [1/f], \nAverage energy {avegene:.2f} eV \nFinal energy " + "$E_f=$" + fr"{finaene[0]:.2f} eV")
        axs[2].set_ylabel('r-axis [m]')
        axs[2].margins(x=0.01, y=0.01)

        fig.tight_layout(pad=2.5)
        plt.show()


def sensitivity():
    fig, ax = plt.subplots()
    x = [1, 2, 3, 4, 5, 6, 7]
    s = [-0.987984504, 0.21315, -0.21626, 0.086534, 0.45319, 0.275492, -1.25225]
    labels = ["$A$", "$B$", "$a$", "$b$", "$R_\mathrm{i}$", "$L$", "$R_\mathrm{eq}$"]
    ax.bar(x, s, align='center', width=1, color=['#1f77b4' if v < 0 else '#ff7f0e' for v in s])
    ax.set_xticks(x, labels)
    ax.set_ylabel(r"$\mathrm{\Delta}f/\mathrm{\Delta}p_i$ [MHz/mm]")
    ax.set_title(r"Sensi½tivity of $f_{\mathrm{FM}}$ [MHz] to geometric variables $p_i$ [mm]")
    for bars in ax.containers:
        ax.bar_label(bars, fontsize=18, label_type='center')


def plot_cavity():
    # data = pd.read_csv(r"D:\Dropbox\CavityDesignHub\C1092V\PostprocessingData\Data\3794_geom.txt", sep='\s+', header=None)
    # data1 = pd.read_csv(r"D:\Dropbox\CavityDesignHub\C1092V\PostprocessingData\Data\2183_geom.txt", sep='\s+', header=None)
    # data3 = pd.read_csv(r"D:\Dropbox\CavityDesignHub\C1092V\PostprocessingData\Data\650_geom.txt", sep='\s+', header=None)
    # data4 = pd.read_csv(r"D:\Dropbox\CavityDesignHub\C1092V\PostprocessingData\Data\770_geom.txt", sep='\s+', header=None)

    # ll = [650, 770, 2183, 3345, 3794, 4123, 4250, 4618]
    ll = ['C40866_geom', 'C3794_800MHz_geom', "G6_C170_M_geom"]
    laf = ['C40866', 'C3794_800MHz', "G6_C170_M"]
    for i, x in enumerate(ll):
        data = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\C800MHz\PostprocessingData\Data\{x}.txt", sep='\s+',
                           header=None)
        plt.rcParams["figure.figsize"] = (5, 5)
        plt.plot(data[1] * 1000, data[0] * 1000, lw=5, label=laf[i], ls='--')
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


def get_cell_data_array_from_vtk(file_path, name):
    cell_data_values = []

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_path)
    reader.Update()

    data = reader.GetOutput()
    cell_data = data.GetCellData()
    cell_data_array = cell_data.GetArray(name)

    for i in range(cell_data_array.GetNumberOfTuples()):
        value = cell_data_array.GetValue(i)
        cell_data_values.append(value)

    return np.array(cell_data_values)


def get_color(val, min_, max_):
    # Normalize the desired value to the range [0, 1] to map it to the colormap
    normalized_value = (np.log(val) - np.log(min_)) / (np.log(max_) - np.log(min_))

    # Get the colors from the "jet" colormap
    color_map = cm.get_cmap('coolwarm')
    desired_color = color_map(normalized_value)

    return desired_color


def spark_3D_statistics(Eacc_cst):
    plot_settings()
    # The source file
    files_folder = fr"D:\CST Studio\3. W\Multipacting\Spark3D\Project 1\Results\@Mod1\@ConfGr1\@EMConfGr1\@MuConf1\Mesh2\Field2"
    file_folders = get_files_from_directories(files_folder, 'Surface3D')

    # read input SEY
    input_sey = pd.read_csv(fr"D:\CST Studio\3. W\Multipacting\Spark3D\user_defined_sey.csv", sep='\s+', header=None)

    # Define the range of values
    min_value = 0.72
    max_value = 1300

    fig, ax = plt.subplots()
    E_ave_list, sey_ave_list, Eacc_list, power_list = [], [], [], []
    for file in file_folders:
        power = float(file.split('_')[-2])
        power_list.append(power)
        Eacc = np.sqrt(2) * np.sqrt(power) * Eacc_cst * 1e-6
        Eacc_list.append(Eacc)
        impact_density = get_cell_data_array_from_vtk(file, "Impact_Density")
        avg_sey = get_cell_data_array_from_vtk(file, "Avg_SEY")
        avg_impact_energy = get_cell_data_array_from_vtk(file, "Avg_Impact_Energy")
        ic(len(avg_impact_energy))
        emission_density = get_cell_data_array_from_vtk(file, "Emission_Density")

        desired_color = get_color(power, min_value, max_value)

        ax.scatter(avg_impact_energy, avg_sey, color=desired_color, edgecolor='k', zorder=power)
        ax.set_xlabel('Electron impact energy [eV]')
        ax.set_ylabel('SEY')

        E_ave = np.sum(np.dot(impact_density, avg_impact_energy)) / np.sum(impact_density)
        sey_ave = np.sum(np.dot(impact_density, avg_sey)) / np.sum(impact_density)

        E_ave_list.append(E_ave)
        sey_ave_list.append(sey_ave)

    ax.plot(input_sey[0], input_sey[1], lw=2,
            path_effects=[patheffects.Stroke(linewidth=4, foreground='k'), patheffects.Normal()], zorder=10000)
    plt.xlim(0, 1000)
    plt.show()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Eacc_list_sorted = np.sort(Eacc_list)
    ind = np.array(Eacc_list).argsort()
    ax.plot(np.sort(power_list), np.array(sey_ave_list)[ind], label='Avg. SEY', marker='o', mec='k')

    ax2 = ax.twinx()
    ax2.plot(np.sort(power_list), np.array(E_ave_list)[ind], label='Avg. Impact Energy', marker='o', mec='k')
    fig.legend()
    ax2.set_ylim(0, 1000)
    plt.show()

    plt.scatter(np.array(E_ave_list)[ind], np.array(sey_ave_list)[ind])
    plt.xlim(0, 1000)
    plt.show()


# spark_3D_statistics(829326.1902)
def butter_lowpass_filter(data, cutoff, fs, order):
    nyq = 0.5 * fs  # Nyquist Frequency
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y


def plot_cst_pic_results(folders, label, Eaccs_eigen, Eacc, xlim=None, ylim=None):
    ########################### CST PIC ####################################
    # c3894
    figsize = (10, 3)
    fig, ax = plt.subplots(figsize=figsize)
    pffs = [fr"{folder}\Export\Particle vs. Time_0D_yAtX.txt" for folder in folders]

    if xlim:
        kwargs = {"width": 0.007143 * (xlim[1] - xlim[0]), "alpha": 1, "ec": 'k', "lw": 0.5}
    else:
        kwargs = {"width": 0.1, "alpha": 1, "ec": 'k', "lw": 0.5}

    for Eacc_eigen, pff in zip(Eaccs_eigen, pffs):
        particle_vs_ff = pd.read_csv(pff, sep='\s+', names=['field_factor', 'particles'], header=None)
        particle_vs_ff = particle_vs_ff.sort_values('field_factor')
        # plt.plot(particle_vs_ff['field_factor'] * Eacc_eigen * 1e-6, particle_vs_ff[result_name], marker='o',
        #              mec='k')
        ax.bar(particle_vs_ff['field_factor'] * Eacc_eigen * 1e-6, particle_vs_ff['particles'], **kwargs)

    ax.set_xlabel(r'$E_\mathrm{acc}$ [MV/m]')
    ax.set_ylabel('Particles')
    if xlim:
        ax.set_xlim(xlim[0], xlim[1])

    ax.axvline(Eacc, c='k', ls='--', lw=1.5)
    ax.text(Eacc, 0.5, f"{np.round(Eacc, 2)} MV/m",
            size=12, rotation=90,
            transform=ax.get_xaxis_transform(), ha='right', va='center')

    plt.minorticks_on()

    fig.tight_layout()
    fig.savefig(fr'D:\Dropbox\Quick presentation files\Multipacting_{label}_PIC_no_of_particles_9.9s_va_eacc.png')
    plt.show()
    plt.close(fig)

    figsize = (10, 3)
    fig, ax = plt.subplots(figsize=figsize)
    for Eacc_eigen, folder in zip(Eaccs_eigen, folders):
        particle_vs_ff = pd.read_csv(fr"{folder}\Export\Particle vs. Time_0D_yAtX.txt", sep='\s+',
                                     names=['field_factor', 'particles'], header=None)
        # get secondary emisison and low pass filter

        # Filter requirements.
        T = 1  # Sample Period
        fs = 50000  # sample rate, Hz
        cutoff = 2 * np.pi / T  # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz

        order = 2  # sin wave can be approx represented as quadratic
        n = int(T * fs)  # total number of samples

        # multipacting intervals
        interval = [6, 10]
        b = []
        for ff in particle_vs_ff['field_factor']:
            if ff.is_integer():
                ff = int(ff)

            pemission = fr"{folder}\\Export\\Emitted Secondaries_(field_factor={ff}).txt"
            sec_particles_vs_time = pd.read_csv(pemission, names=['time', 'particles'], header=None, sep='\s+')

            # low pass filter
            sec_low_pass_filter = butter_lowpass_filter(sec_particles_vs_time['particles'], cutoff, fs, order)
            sec_particles_vs_time['particles_filtered'] = sec_low_pass_filter

            # fitting
            sec_particles_vs_time_int = sec_particles_vs_time.loc[
                (sec_particles_vs_time['time'] > interval[0]) & (sec_particles_vs_time['time'] < interval[1])]

            x, y = sec_particles_vs_time_int['time'], sec_particles_vs_time_int['particles_filtered']
            sec_particles_vs_time_int['log_y'] = np.log(y)
            sec_particles_vs_time_int = sec_particles_vs_time_int.dropna(axis=0)

            x, y = sec_particles_vs_time_int['time'], sec_particles_vs_time_int['particles_filtered']

            AB = np.polyfit(x, sec_particles_vs_time_int['log_y'], 1)
            if np.isnan(AB[0]) or np.isnan(AB[1]):
                AB = [0, 0]
            else:
                AB = [np.exp(AB[0]), np.exp(AB[1])]

            b.append(AB[0])
            # ic(ff, AB[0])

        ax.bar(particle_vs_ff['field_factor'] * Eacc_eigen * 1e-6, b, label='Exponent base b', **kwargs)
        if ylim:
            ax.set_ylim(ylim[0], ylim[1])
        else:
            ax.set_ylim(0.75, max(b))

    if xlim:
        ax.set_xlim(xlim[0], xlim[1])

    ax.axvline(Eacc, c='k', ls='--', lw=1.5)
    ax.text(Eacc, 0.5, f"{np.round(Eacc, 2)} MV/m",
            size=12, rotation=90,
            transform=ax.get_xaxis_transform(), ha='right', va='center')
    ax.axhline(1, c='r')
    ax.set_xlabel(r'$E_\mathrm{acc}$ [MV/m]')
    ax.set_ylabel('Exponential growth factor')
    fig.tight_layout()
    plt.minorticks_on()
    fig.savefig(fr'D:\Dropbox\Quick presentation files\Multipacting_{label}_PIC_exponent_base_vs_eacc.png')
    plt.show()


if __name__ == '__main__':
    # plt.style.use('fivethirtyeight')
    plot_settings()

    # ########################### Multipac ####################################
    # # C3794
    # folders = [r"D:\Dropbox\multipacting\MPGUI21\C3794"]
    # labels = ['C3794']
    # Epk_Eacc_list = [2.05]
    # Eacc_list = [10.61]
    # min_max = [1, 15]
    #
    # for k in ['counter_function', 'final_impact_energy', 'enhanced_counter_function']:
    #     plot_multipac_triplot(Eacc_list, Epk_Eacc_list, folders, labels, kind=k, layout=[[0]], min_max=min_max)
    #
    # # C3795
    # folders = [r"D:\Dropbox\multipacting\MPGUI21\Mid_cell_C3795", r"D:\Dropbox\multipacting\MPGUI21\End_cell_C3795"]
    #
    # labels = ['C3795 (mid cell)', 'C3795 (end cell)']
    # Epk_Eacc_list = [2.31, 2.56]
    # Eacc_list = [20.12, 20.12]
    # min_max = [1, 26]
    # # fig, axs = plt.subplots(3)
    #
    # for k in ['counter_function', 'final_impact_energy', 'enhanced_counter_function']:
    #     plot_multipac_triplot(Eacc_list, Epk_Eacc_list, folders, labels, kind=k, layout=[[0]], min_max=min_max)

    files_folder = r"D:\Dropbox\multipacting\MPGUI21\C3794"
    plot_trajectory(files_folder, loc='left')

    files_folder = r"D:\Dropbox\multipacting\MPGUI21\Mid_cell_C3795"
    plot_trajectory(files_folder)

    ######################### CST PIC #####################################
    # c3794
    folders = [fr'D:\CST Studio\3. W\Multipacting\MP_C3794_end_cell_hex']
    Eaccs_eigen = [1.208299394e6]
    plot_cst_pic_results(folders, label='C3794', Eaccs_eigen=Eaccs_eigen, Eacc=10.61, xlim=[1, 15], ylim=[0.75, 1.3])

    # c3795
    folders = [fr'D:\CST Studio\5. tt\Multipacting\Sturmvogel\c3795_mid_cell',
               fr'D:\CST Studio\5. tt\Multipacting\Sturmvogel\c3795_end_cell']
    Eaccs_eigen = [3.650251e6, 3.672625e6]
    plot_cst_pic_results(folders, label='C3795', Eaccs_eigen=Eaccs_eigen, Eacc=20.12, xlim=[1, 26])

    ############################### Spark3D ########################################################
    filename = r'D:\CST Studio\3. W\Multipacting\Spark3D\MultipactorConfig 1.csv'
    power_folder = fr'D:\CST Studio\3. W\Multipacting\Spark3D\Project 1\Results\@Mod1\@ConfGr1\@EMConfGr1\@MuConf1\Mesh2\Field2'
    plot_electron_evolution_spark3d(filename, power_folder, 829326.1902)
