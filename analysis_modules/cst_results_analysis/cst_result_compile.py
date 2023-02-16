import os.path

import numpy as np
import pandas as pd


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

    df = df.sort_values(by=['Frequency (Multiple Modes)'])
    recursive_save(df, fr"{sim_folder}\Lossy_Eigenmode_results")


def write_threshold(filename):
    # get frequency of max R/Q [kOhm/m] mode
    df_threshold = pd.DataFrame()
    Q = 1e4
    df_threshold["f [MHz]"] = np.linspace(0, 5000, 1000)

    df_threshold['f/Q'] = df_threshold["f [MHz]"] * 1e6 / Q
    df_threshold['ZTth [kOhm/m]'] = (100 * 1e3 * Q / (df_threshold["f [MHz]"] * 1e-3) ** 2)/670  # normalize by number of cavities
    df_threshold['ZTth [kOhm/m]'][(df_threshold['f/Q'] < 1e5) & (1e10 / 670 < df_threshold['ZTth [kOhm/m]'])] = 1e10 / 670  # normalize by number of cavities
    df_threshold.to_excel(f"{filename}_threshold.xlsx")


def s_parameters(sim_folder, folders, requests):
    for folder in folders:
        df = pd.DataFrame()
        for request in requests:
            if 'f [MHz]' in df.keys():
                d = pd.read_csv(fr"{sim_folder}\{folder}\Export\{request}.txt", sep="\t", header=None)
                re = d.loc[:, 1]
                im = d.loc[:, 2]
                mag = np.sqrt(re**2+im**2)
                magdB = 20*np.log10(mag)
                df[request] = magdB
            else:
                d = pd.read_csv(fr"{sim_folder}\{folder}\Export\{request}.txt", sep="\t", header=None)
                df['f [MHz]'] = d.loc[:, 0]
                re = d.loc[:, 1]
                im = d.loc[:, 2]
                mag = np.sqrt(re**2+im**2)
                magdB = 20*np.log10(mag)
                df[request] = magdB

        df.to_excel(fr"{sim_folder}\{folder}.xlsx")


def recursive_save(df, filename):
    if os.path.exists(f"{filename}.xlsx"):
        filename = filename + '(1)'
        recursive_save(df, filename)
    else:
        df.to_excel(f"{filename}.xlsx")


if __name__ == '__main__':
    # sim_folder = r"D:\CST Studio\MuCol_Study\Assembly\Eigenmode"
    # folders = ["NLSF_9_cell_2DQW_D1", "NLSF_9_cell_2DQW_D2", "NLSF_9_cell_2DQW_D3", "NLSF_9_cell_2DQW_M1",
    #            "NLSF_9_cell_2DQW_M2"]
    # req = ["Frequency (Multiple Modes)", "Q-Factor (lossy E) (Multiple Modes)", "R over Q beta=1 (Multiple Modes)",
    #        "RQT", "RQT_kOhm_m", "Z_kOhm", "Z_T_kOhm_m"]
    # eigenmode_analysis(sim_folder, folders, req)

    sim_folder = r"D:\CST Studio\MuCol_Study\Couplers"
    folders = ["DQW_1300MHz_Ri38mm"]
    req = ["S-Parameters_S1(1),2(1)", "S-Parameters_S1(2),2(1)", "S-Parameters_S1(3),2(1)",
           "S-Parameters_S2(1),2(1)", "S-Parameters_S3(1),2(1)", "S-Parameters_S3(2),2(1)", "S-Parameters_S3(3),2(1)"]
    s_parameters(sim_folder, folders, req)
