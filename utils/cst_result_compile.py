import pandas as pd


def eigenmode_analysis(folder, requests):
    df = pd.DataFrame()
    for request in requests:
        d = pd.read_csv(fr"{folder}\{request}.txt", sep="\t", header=None)
        df[request] = d.loc[:, 1]

    df.to_excel(fr"{folder}\Lossy_Eigenmode_results.xlsx")


req = ["Frequency (Multiple Modes)", "Q-Factor (lossy E) (Multiple Modes)", "R over Q beta=1 (Multiple Modes)",
       "RQT", "Z_kOhm", "Z_T_kOhm_m"]
eigenmode_analysis(r"D:\CST Studio\5. tt\Eigenmode\E_C3795_from_1700MHz\Export", req)
