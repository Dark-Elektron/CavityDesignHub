import math
import os.path
import matplotlib as mpl
import numpy as np
import pandas as pd
import scipy.integrate
from icecream import ic
from matplotlib import pyplot as plt
from numpy import loadtxt


def compile_monopole(folder, request, cavity_half_len, n_cells):
    """
    This function filters monopole modes from eigenmode analysis
    Returns
    -------

    """

    monopole_mode_num, beampipe_modes = [], []

    for path in os.listdir(fr'{folder}\Export'):
        if request in path:
            # get axis field profile

            mode_number = int(path.split(' ')[-1].split('.')[0])

            mode = pd.read_csv(fr"{folder}\Export\{path}", sep="\t", header=None)

            z, Ez = mode[0], mode[1]

            # filter axis field if maximum is less than 1e4
            if max(abs(Ez)) < 1e4:
                # ic('Low axis field value: ', mode_number, max(abs(Ez)))
                continue

            # # detect noisy field values
            # peaks, _ = scipy.signal.find_peaks(abs(Ez), height=2000)
            # if len(peaks) > 20:
            #     # ic('Noisy: Possibly Quadrupole: ', mode_number, "Number of peaks: ", len(peaks))
            #     continue

            # check if maximum field is in beampipe
            Ez_cav = np.array(Ez)[np.where((z > -cavity_half_len*n_cells) & (z < cavity_half_len*n_cells))]
            Ez_bp = np.array(Ez)[np.where((z < -cavity_half_len*n_cells) | (z > cavity_half_len*n_cells))]
            # plt.plot(np.array(z)[np.where((z < -cavity_half_len*n_cells) | (z > cavity_half_len*n_cells))], Ez_bp)
            # plt.plot(np.array(z)[np.where((z > -cavity_half_len*n_cells) & (z < cavity_half_len*n_cells))], Ez_cav)
            # plt.show()

            Ez_cav_mx, Ez_bp_mx = max(abs(Ez_cav)), max(abs(Ez_bp))
            # ic(mode_number, Ez_cav_mx, Ez_bp_mx)
            field_ratio = Ez_cav_mx/Ez_bp_mx
            # ic(mode_number, field_ratio)
            if field_ratio < 0.2:
                beampipe_modes.append(mode_number)
                # ic('Beampipe mode: ', mode_number)
                continue

            monopole_mode_num.append(mode_number)
    sorted_mode_number = np.sort(monopole_mode_num)
    ic(np.sort(beampipe_modes))
    return sorted_mode_number


def monopole_mode(z, Ez, cavity_half_len, n_cells):
    # filter axis field if maximum is less than 1e4
    if max(abs(Ez)) < 1e5:
        return False

    # # detect noisy field values
    # peaks, _ = scipy.signal.find_peaks(abs(Ez))
    # if len(peaks) > 20:
    #     return False

    # check if maximum field is in beampipe
    Ez_cav = np.array(Ez)[np.where((z > -cavity_half_len * n_cells) & (z < cavity_half_len * n_cells))]
    Ez_bp = np.array(Ez)[np.where((z < -cavity_half_len * n_cells) | (z > cavity_half_len * n_cells))]

    Ez_cav_mx, Ez_bp_mx = max(abs(Ez_cav)), max(abs(Ez_bp))
    field_ratio = Ez_cav_mx / Ez_bp_mx
    if field_ratio < 0.2:
        return False

    return True


def dipole_mode(z, Ez, cavity_half_len, n_cells):
    Ez = np.abs(Ez)
    # filter axis field if maximum is less than 1e4
    if max(abs(Ez)) < 1e5:
        ic('Low axis field value: ', max(abs(Ez)))
        return False

    # detect noisy field values
    peaks, _ = scipy.signal.find_peaks(abs(Ez), height=2000)
    # plt.plot(z, Ez)
    # plt.scatter(z[peaks], Ez[peaks])
    # plt.show()
    if len(peaks) > 100:
        ic('Noisy: Possibly Quadrupole: ', "Number of peaks: ", len(peaks))
        return False

    # check if maximum field is in beampipe
    Ez_cav = np.array(Ez)[np.where((z > -cavity_half_len * n_cells) & (z < cavity_half_len * n_cells))]
    Ez_bp = np.array(Ez)[np.where((z < -cavity_half_len * n_cells) | (z > cavity_half_len * n_cells))]
    # plt.plot(np.array(z)[np.where((z < -cavity_half_len*n_cells) | (z > cavity_half_len*n_cells))], Ez_bp)
    # plt.plot(np.array(z)[np.where((z > -cavity_half_len*n_cells) & (z < cavity_half_len*n_cells))], Ez_cav)
    # plt.show()
    Ez_cav_mx, Ez_bp_mx = max(abs(Ez_cav)), max(abs(Ez_bp))
    field_ratio = Ez_cav_mx / Ez_bp_mx
    if field_ratio < 0.2:
        ic('Beampipe mode: ')
        return False

    return True


def eigenmode_analysis(sim_folder, folders, requests, name):
    df = pd.DataFrame()
    for request in requests:
        d = pd.DataFrame()
        for folder in folders:
            new_d = pd.read_csv(fr"{sim_folder}\{folder}\Export\{request}.txt", sep="\t", header=None)
            new_d.index += 1
            d = pd.concat([d, new_d])

        df[request] = d.loc[:, 1]

    # # get monopole mode inde and filter dataframe
    # monopole_index = compile_monopole(fr"{sim_folder}\{name}", 'e_Z (Z)_', 93.5, 5) - 1
    # ic(monopole_index)
    # df = df.iloc[monopole_index, :]

    # # other
    # f = df["Frequency (Multiple Modes)"][df["RQT_kOhm_m"] == max(df["RQT_kOhm_m"])]
    # f = f.values[0]
    # print(f, max(df["RQT_kOhm_m"]))
    # write_threshold(fr"{sim_folder}\Lossy_Eigenmode_results")

    # df = df.sort_values(by=['Frequency (Multiple Modes)'])
    recursive_save(df, fr"{sim_folder}\{name}")


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
                # ax.axhline(1e10/670, c='k')

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
                # ax.axhline(1e10/670, c='k')

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
                # ax.axhline(1e10/670, c='k')

    f.supxlabel("$f$ [MHz]")
    f.supylabel(r'$Z_\perp \mathrm{/cav ~[k \Omega/m]}$')
    # f.supylabel(r'$R/Q_\perp \cdot n_\mathrm{cavs} \cdot f^2 ~\mathrm{[M \Omega/m \cdot GHz^2]}$')


def plot_threshold():
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
    plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
              mode="expand", borderaxespad=0, ncol=2)
    plt.tight_layout()
    plt.show()


def plot_threshold_Z():
    df = pd.read_excel(fr"D:\CST Studio\MuCol_Study\Assembly\Eigenmode\NLSF_Lossy_Eigenmode_results_.xlsx", "Sheet1")
    intervals = [[1600, 1990], [2450, 2590]]

    df = pd.read_excel(fr"D:\CST Studio\MuCol_Study\Assembly\Eigenmode\ERL_MA_Lossy_Eigenmode_results.xlsx", "Sheet1")
    intervals = [[1540, 1920], [2470, 2515]]

    # df = pd.read_excel(fr"D:\CST Studio\MuCol_Study\Assembly\Eigenmode\TESLA_Lossy_Eigenmode_results.xlsx", "Sheet1")
    # intervals = [[1600, 1905], [2470, 2590]]

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
    plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
              mode="expand", borderaxespad=0, ncol=2)
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


def plot_particle_vs_time():

    df = pd.read_csv(fr"D:\CST Studio\Multipacting_CST\tesla_mid_cell.csv", sep=",", skipfooter=1, engine='python')
    df1 = pd.read_csv(fr"D:\CST Studio\Multipacting_CST\tesla_mid_cell_ref.csv", sep=",", skipfooter=1, engine='python')
    df = pd.concat([df, df1])
    ic(df)

    plt.stem(df['field_factor'], df[r'Tables\0D Results\Particle vs. Time_0D_yAtX'])
    plt.show()

# plot_particle_vs_time()


def integrate_wake_potential():
    sigma = 100e-3  # m
    WL = pd.read_csv(fr"D:\CST Studio\TESLA\TESLA\Export\Particle Beams_ParticleBeam1_Wake potential_Z.txt", sep='\s+', engine='python', header=None)
    WTx = pd.read_csv(fr"D:\CST Studio\TESLA\TESLA\Export\Particle Beams_ParticleBeam1_Wake potential_X.txt", sep='\s+', engine='python', header=None)
    WTy = pd.read_csv(fr"D:\CST Studio\TESLA\TESLA\Export\Particle Beams_ParticleBeam1_Wake potential_X.txt", sep='\s+', engine='python', header=None)

    WT = (WTx + WTy)/2
    s = WL[0]*1e-3  # convert to m

    # scale with bunch length
    # for gaussian bunch
    lambda_s = 1/(np.sqrt(2*np.pi)*sigma) * np.exp(-s**2/(2*sigma**2))

    WsL = np.multiply(lambda_s, WL[1])
    plt.plot(s, WsL)
    plt.show()
    WsT = np.multiply(lambda_s, WT[1])

    k_loss = -np.trapz(WsL, s)
    k_kick = np.trapz(WsT, s)

    ds = np.array(s[1:]) - np.array(s[:-1])
    left_area = WsL[:-1] * ds
    right_area = WsL[1:] * ds
    left_int = math.fsum(left_area)
    right_int = math.fsum(right_area)

    ic(left_int, right_int, np.abs(left_int-right_int))

    ic(k_loss, k_kick)


# Define a function to find approximate matches
def find_approx_matches(df1, df2, col_name, tolerance):
    val_1, val_2 = df1[col_name], df2[col_name]

    close_freq_list = []
    for v1 in val_1:
        subtr = np.abs(val_2 - v1)
        ind, min_diff = np.argmin(subtr), np.min(subtr)

        if min_diff <= tolerance:
            close_freq_list.append(v1)

    ic(close_freq_list)
    return df1[df1[col_name].isin(close_freq_list)]


def plot_frequency_line(ee_file, mm_file):
    fig, axs = plt.subplot_mosaic([[0, 1, 2], [3, 4, 5]], figsize=(14, 4))
    df_ee = pd.read_excel(ee_file, 'Sheet1')[['Frequency (Multiple Modes)', 'Pol']]
    df_mm = pd.read_excel(mm_file, 'Sheet1')[['Frequency (Multiple Modes)', 'Pol']]
    f_ee = df_ee['Frequency (Multiple Modes)']
    f_mm = df_mm['Frequency (Multiple Modes)']

    # get modes with same frqquency
    # Define your tolerance for approximate equality
    tolerance = 0.1
    # Find approximate matches using Pandas merge
    merged_df = pd.merge_asof(df_ee, df_mm, on='Frequency (Multiple Modes)', direction='nearest', tolerance=tolerance)
    # Filter the merged DataFrame to get approximate matches
    approx_matches_pd = merged_df.dropna()
    approx_matches = find_approx_matches(df_ee, df_mm, 'Frequency (Multiple Modes)', tolerance=tolerance)
    f_tp = approx_matches['Frequency (Multiple Modes)']

    # Print the matches
    if not approx_matches_pd.empty:
        f_tp_pd = approx_matches_pd['Frequency (Multiple Modes)']
    else:
        ic("No approximate matches found.")

    intervals = [[370, 600], [600, 1000], [1000, 1300], [1300, 1500], [1500, 1800], [1800, 2000]]

    for i, interval in enumerate(intervals):
        axs[i].scatter(f_ee[(f_ee > interval[0]) & (f_ee < interval[1])], f_ee[(f_ee > interval[0]) & (f_ee < interval[1])], marker='o', edgecolors='k', facecolors='none', s=25)
        axs[i].scatter(f_mm[(f_mm > interval[0]) & (f_mm < interval[1])], f_mm[(f_mm > interval[0]) & (f_mm < interval[1])], marker='+', facecolors='r', s=25)
        axs[i].scatter(f_tp[(f_tp > interval[0]) & (f_tp < interval[1])], f_tp[(f_tp > interval[0]) & (f_tp < interval[1])], marker='o', edgecolors='b', facecolors='none', s=100)
        axs[i].scatter(f_tp_pd[(f_tp_pd > interval[0]) & (f_tp_pd < interval[1])], f_tp_pd[(f_tp > interval[0]) & (f_tp_pd < interval[1])], marker='o', edgecolors='pink', facecolors='pink', s=50)
        ylim = axs[i].get_ylim()
        # plot vertical lines
        for f in f_tp[(f_tp > interval[0]) & (f_tp < interval[1])]:
            axs[i].plot([f, f], [0, f], ls='--', c='k')
        axs[i].set_ylim(ylim)

    plt.tight_layout()
    plt.show()


# integrate_wake_potential()

if __name__ == '__main__':
    plt.rcParams["figure.figsize"] = (10, 3)
    # sim_folder = r"D:\CST Studio\3. W\Eigenmode"
    sim_folder = r"D:\CST Studio\5. tt\Assembly"
    # folders = ["E_C3794_390_665_EE", "E_C3794_665_1000_EE", "E_C3794_1000_1400_EE", "E_C3794_1400_1800_EE", "E_C3794_1800_2200_EE"]
    # folders = ["E_C3794_390_665", "E_C3794_665_1000", "E_C3794_1000_1400_MM", "E_C3794_1400_1800_MM", "E_C3794_1800_2200_MM"]
    # folders = ["E_C3795_700_1250", "E_C3795_1250_1450", "E_C3795_1450_1650", "E_C3795_1650_1850", "E_C3795_1850_2050", "E_C3795_2050_2250"]
    folders = ["E_C3795_2DQW_700_1250_ref_mesh", "E_C3795_2DQW_1250_1450_ref_mesh", "E_C3795_2DQW_1450_1650_ref_mesh",
               "E_C3795_2DQW_1650_1850_ref_mesh", "E_C3795_2DQW_1850_2050_ref_mesh", "E_C3795_2DQW_2050_2250_ref_mesh"]
    req = ["Frequency (Multiple Modes)", "Q-Factor (lossy E) (Multiple Modes)", "R over Q beta=1 (Multiple Modes)",
           "RQT", "Z_kOhm", "Z_T_kOhm_m"]
    # req = ["Frequency (Multiple Modes)", "R over Q beta=1 (Multiple Modes)", "RQT"]
    name = "E_C3795"

    eigenmode_analysis(sim_folder, folders, req, name)

    # sim_folder = r"D:\CST Studio\MuCol_Study\Couplers"
    # folders = ["DQW_1300MHz_Ri38mm"]
    # req = ["S-Parameters_S1(1),2(1)", "S-Parameters_S1(2),2(1)", "S-Parameters_S1(3),2(1)",
    #        "S-Parameters_S2(1),2(1)", "S-Parameters_S3(1),2(1)", "S-Parameters_S3(2),2(1)", "S-Parameters_S3(3),2(1)"]
    # s_parameters(sim_folder, folders, req)

    # plot_threshold()
    # plot_threshold_Z()

    # sim_folder = r"D:\CST Studio\3. W\Eigenmode\E_C3794"
    # monopole = compile_monopole(sim_folder, 'e_Z (Z)_x=0', 187, 2)
    # ic(monopole, len(monopole))
    # monopole = compile_monopole(sim_folder, 'e_Z (Z)_x=Ri', 187, 2)
    # ic(monopole, len(monopole))


    # ee_file = fr'D:\CST Studio\3. W\Eigenmode\E_C3794_EE.xlsx'
    # mm_file = fr'D:\CST Studio\3. W\Eigenmode\E_C3794_MM.xlsx'
    # plot_frequency_line(ee_file, mm_file)

