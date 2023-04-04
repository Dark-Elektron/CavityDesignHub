import os.path
import matplotlib as mpl
import numpy as np
import pandas as pd
from icecream import ic
from matplotlib import pyplot as plt


def eigenmode_analysis(sim_folder, folders, requests):
    df = pd.DataFrame()
    for request in requests:
        d = pd.DataFrame()
        for folder in folders:
            d = pd.concat([d, pd.read_csv(fr"{sim_folder}\{folder}\Export\{request}.txt", sep="\t", header=None)])
            # print(d)

        df[request] = d.loc[:, 1]

    # other
    f = df["Frequency (Multiple Modes)"][df["RQT_kOhm_m"] == max(df["RQT_kOhm_m"])]
    f = f.values[0]
    print(f, max(df["RQT_kOhm_m"]))
    write_threshold(fr"{sim_folder}\Lossy_Eigenmode_results")

    # df = df.sort_values(by=['Frequency (Multiple Modes)'])
    recursive_save(df, fr"{sim_folder}\Lossy_Eigenmode_results")


def split_axis_x(data, bounds, scale='linear', divisions=None, labels=None):
    if labels is None:
        labels = ['']
    if divisions is None:
        divisions = [100]
    f, axs = plt.subplots(1, len(bounds), sharey=True)

    # plot the same data on both axes
    for n, xy in enumerate(data):

        x, y = xy

        for i, ax in enumerate(axs):
            if len(divisions) > 1:
                try:
                    width = (bounds[i][1] - bounds[i][0])/divisions[n]
                except:
                    width = (bounds[i][1] - bounds[i][0])/divisions[-1]
            else:
                width = (bounds[i][1] - bounds[i][0])/divisions[0]

            ax.bar(x, y, width=width, label=labels[n])
            # zoom-in / limit the view to different portions of the data
            ax.set_xlim(bounds[i][0], bounds[i][1])  # outliers only

            d = .015  # how big to make the diagonal lines in axes coordinates
            # arguments to pass to plot, just so we don't keep repeating them
            kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
            # hide the spines between ax and ax2
            if i == 0:
                ax.spines['right'].set_visible(False)

                ax.tick_params(labelright=False)  # don't put tick labels at the right

                ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # top-right diagonal
                ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # bottom-right diagonal
                # ax.axhline(100, c='k')
                ax.axhline(1e10/670, c='k')

            if 0 < i < len(bounds) - 1:
                ax.spines['left'].set_visible(False)
                ax.spines['right'].set_visible(False)

                ax.set_yscale(scale)
                ax.tick_params(axis='y', which='minor', left=False)
                ax.tick_params(axis='y', which='minor', right=False)
                ax.tick_params(top=False, bottom=True, left=False, right=False,
                               labelleft=False, labelbottom=True)

                ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # top-right diagonal
                ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # bottom-right diagonal

                ax.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
                ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
                # ax.axhline(100, c='k')
                ax.axhline(1e10/670, c='k')

            if i == len(bounds) - 1:
                ax.spines['left'].set_visible(False)

                ax.set_yscale(scale)
                # plt.minorticks_off()
                ax.tick_params(axis='y', which='minor', left=False)
                ax.tick_params(top=False, bottom=True, left=False, right=False,
                               labelleft=False, labelbottom=True)

                d = .015  # how big to make the diagonal lines in axes coordinates
                # arguments to pass to plot, just so we don't keep repeating them
                kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
                ax.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
                ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # top-right diagonal
                # ax.axhline(100, c='k')
                ax.axhline(1e10/670, c='k')

    f.supxlabel("$f$ [MHz]")
    f.supylabel(r'$Z_\perp \mathrm{/cav ~[k \Omega/m]}$')
    # f.supylabel(r'$R/Q_\perp \cdot n_\mathrm{cavs} \cdot f^2 ~\mathrm{[M \Omega/m \cdot GHz^2]}$')


def plot_threshold():
    plt.rcParams["figure.figsize"] = (10, 4.5)
    df = pd.read_excel(fr"D:\CST Studio\MuCol_Study\Assembly\Eigenmode\NLSF_Lossy_Eigenmode_results_.xlsx", "Sheet1")
    intervals = [[1600, 1990], [2450, 2590]]

    df = pd.read_excel(fr"D:\CST Studio\MuCol_Study\Assembly\Eigenmode\ERL_MA_Lossy_Eigenmode_results.xlsx", "Sheet1")
    intervals = [[1540, 1920], [2470, 2515]]

    df = pd.read_excel(fr"D:\CST Studio\MuCol_Study\Assembly\Eigenmode\TESLA_Lossy_Eigenmode_results.xlsx", "Sheet1")
    intervals = [[1600, 1905], [2470, 2590]]

    # get only dipole modes
    df = df[df['Pol'] == 'D']
    df = df.sort_values(by=['Frequency (Multiple Modes)'])
    f = df["Frequency (Multiple Modes)"]  # MHz
    ic(f)
    Q = df["Q-Factor (lossy E) (Multiple Modes)"]
    f_Q = f * 1e6 / Q

    # Zth = (100 * 1e3 * Q / (f * 1e-3) ** 2) / 670
    # Zth[f_Q < 1e5] = 1e10 / 670
    df['Single Turn'] = 100  # [MOhm/m x GHz^2]
    df['Multi Turn'] = 1e7 * (f*1e-3)**2 / Q  # [MOhm/m x GHz^2]
    Zth = df[['Single Turn', 'Multi Turn']].min(axis=1)

    divisions = [30, 80]  # to ensure the bars have equal width
    labels = ['$R/Q_\mathrm{\perp, threshold} \cdot f^2$', '$R/Q_\perp \cdot f^2$ (Lossy eig.)']

    Z = (df['RQT_kOhm_m']*1e-3*670) * (f*1e-3)**2
    split_axis_x([[f, Zth], [f, Z]], intervals, scale='log', divisions=divisions, labels=labels)
    plt.ylim(0, 1e3)
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.show()


def plot_threshold_Z():
    plt.rcParams["figure.figsize"] = (10, 5)
    df = pd.read_excel(fr"D:\CST Studio\MuCol_Study\Assembly\Eigenmode\NLSF_Lossy_Eigenmode_results_.xlsx", "Sheet1")
    intervals = [[1600, 1990], [2450, 2590]]

    df = pd.read_excel(fr"D:\CST Studio\MuCol_Study\Assembly\Eigenmode\ERL_MA_Lossy_Eigenmode_results.xlsx", "Sheet1")
    intervals = [[1540, 1920], [2470, 2515]]

    df = pd.read_excel(fr"D:\CST Studio\MuCol_Study\Assembly\Eigenmode\TESLA_Lossy_Eigenmode_results.xlsx", "Sheet1")
    intervals = [[1600, 1905], [2470, 2590]]

    # get only dipole modes
    df = df[df['Pol'] == 'D']
    df = df.sort_values(by=['Frequency (Multiple Modes)'])
    f = df["Frequency (Multiple Modes)"]  # MHz
    ic(f)
    Q = df["Q-Factor (lossy E) (Multiple Modes)"]
    f_Q = f * 1e6 / Q

    # Zth = (100 * 1e3 * Q / (f * 1e-3) ** 2) / 670
    # Zth[f_Q < 1e5] = 1e10 / 670
    df['Single Turn'] = 100 * 1e3 * Q / (f*1e-3)**2 / 670  # [MOhm/m]
    df['Multi Turn'] = 1e10 / 670  # [MOhm/m]
    Zth = df[['Single Turn', 'Multi Turn']].min(axis=1)

    divisions = [30, 80]  # to ensure the bars have equal width
    labels = ['$Z_\mathrm{\perp, threshold}$', '$Z_\perp$ (Lossy eig.)']

    split_axis_x([[f, Zth], [f, df['Z_T_kOhm_m']]], intervals, scale='log', divisions=divisions, labels=labels)
    plt.ylim(0, 1e11 / 670)
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.show()


def write_threshold(filename):
    # get frequency of max R/Q [kOhm/m] mode
    df_threshold = pd.DataFrame()
    Q = 1e8
    df_threshold["f [MHz]"] = np.linspace(0, 5000, 1000)

    df_threshold['f/Q'] = df_threshold["f [MHz]"] * 1e6 / Q
    df_threshold['ZTth [kOhm/m]'] = (100 * 1e3 * Q / (
            df_threshold["f [MHz]"] * 1e-3) ** 2) / 670  # normalize by number of cavities
    df_threshold['ZTth [kOhm/m]'][(df_threshold['f/Q'] < 1e5) & (
            1e10 / 670 < df_threshold['ZTth [kOhm/m]'])] = 1e10 / 670  # normalize by number of cavities
    df_threshold.to_excel(f"{filename}_threshold.xlsx")


def s_parameters(sim_folder, folders, requests):
    for folder in folders:
        df = pd.DataFrame()
        for request in requests:
            if 'f [MHz]' in df.keys():
                d = pd.read_csv(fr"{sim_folder}\{folder}\Export\{request}.txt", sep="\t", header=None)
                re = d.loc[:, 1]
                im = d.loc[:, 2]
                mag = np.sqrt(re ** 2 + im ** 2)
                magdB = 20 * np.log10(mag)
                df[request] = magdB
            else:
                d = pd.read_csv(fr"{sim_folder}\{folder}\Export\{request}.txt", sep="\t", header=None)
                df['f [MHz]'] = d.loc[:, 0]
                re = d.loc[:, 1]
                im = d.loc[:, 2]
                mag = np.sqrt(re ** 2 + im ** 2)
                magdB = 20 * np.log10(mag)
                df[request] = magdB

        df.to_excel(fr"{sim_folder}\{folder}.xlsx")


def recursive_save(df, filename):
    if os.path.exists(f"{filename}.xlsx"):
        filename = filename + '(1)'
        recursive_save(df, filename)
    else:
        df.to_excel(f"{filename}.xlsx")


if __name__ == '__main__':
    sim_folder = r"D:\CST Studio\MuCol_Study\Assembly\Eigenmode"
    folders = ["NLSF_9_cell_2DQW_D1", "NLSF_9_cell_2DQW_D2", "NLSF_9_cell_2DQW_D3", "NLSF_9_cell_2DQW_M1",
               "NLSF_9_cell_2DQW_M2"]
    folders = ["ERL_MA_9_cell_2DQW_D1", "ERL_MA_9_cell_2DQW_D2", "ERL_MA_9_cell_2DQW_D3", "ERL_MA_9_cell_2DQW_M1", "ERL_MA_9_cell_2DQW_M2"]
    folders = ["TESLA_9_cell_2DQW_D1", "TESLA_9_cell_2DQW_D2", "TESLA_9_cell_2DQW_M1", "TESLA_9_cell_2DQW_D3"]
    req = ["Frequency (Multiple Modes)", "Q-Factor (lossy E) (Multiple Modes)", "R over Q beta=1 (Multiple Modes)",
           "RQT", "RQT_kOhm_m", "Z_kOhm", "Z_T_kOhm_m"]
    eigenmode_analysis(sim_folder, folders, req)

    # sim_folder = r"D:\CST Studio\MuCol_Study\Couplers"
    # folders = ["DQW_1300MHz_Ri38mm"]
    # req = ["S-Parameters_S1(1),2(1)", "S-Parameters_S1(2),2(1)", "S-Parameters_S1(3),2(1)",
    #        "S-Parameters_S2(1),2(1)", "S-Parameters_S3(1),2(1)", "S-Parameters_S3(2),2(1)", "S-Parameters_S3(3),2(1)"]
    # s_parameters(sim_folder, folders, req)

    # plot_threshold()
    # plot_threshold_Z()
