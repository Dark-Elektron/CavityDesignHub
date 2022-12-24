import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from icecream import ic

if __name__ == '__main__':

    n_modes = 10
    fig, axs = plt.subplots(2, n_modes, sharex=True, sharey=True, figsize=(15, 4))
    for nm in range(1, n_modes+1):
        result_folder_reference = fr'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\E_3794_4HC_1FPC_Optimize_Reference\Export'
        df_u = pd.read_csv(fr"{result_folder_reference}\3d\Mode {nm}_e.txt", sep=' \s+', header=0, skiprows=[1], engine='python')
        freq_reference = pd.read_csv(fr"{result_folder_reference}\Frequency (Multiple Modes).txt", sep='\t', header=None)

        # ic(df_u.shape)
        df_u['r'] = np.sqrt(df_u['x [mm]']**2 + df_u['y [mm]']**2)
        df_u['Emag [V/m]'] = np.sqrt(df_u['ExRe [V/m]']**2 + df_u['ExIm [V/m]']**2 +
                                     df_u['EyRe [V/m]']**2 + df_u['EyIm [V/m]']**2 +
                                     df_u['EzRe [V/m]']**2 + df_u['EzIm [V/m]']**2)
        df_u['|Emag|'] = df_u['Emag [V/m]']/df_u['Emag [V/m]'].max()  # normalize

        F_u = df_u.loc[df_u['r'] <= 150].reset_index(drop=True)
        # ic(F_u.shape)

        result_folder = fr'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\E_3794_4HC_1FPC_Optimize\Export\3d'
        df = pd.read_csv(fr"{result_folder}\Mode {nm}_e.txt", sep=' \s+', header=0, skiprows=[1], engine='python')
        freq = pd.read_csv(fr"{result_folder_reference}\Frequency (Multiple Modes).txt", sep='\t', header=None)
        df['Emag [V/m]'] = np.sqrt(df['ExRe [V/m]']**2 + df['ExIm [V/m]']**2 +
                                   df['EyRe [V/m]']**2 + df['EyIm [V/m]']**2 +
                                   df['EzRe [V/m]']**2 + df['EzIm [V/m]']**2)
        df['|Emag|'] = df['Emag [V/m]']/df['Emag [V/m]'].max()  # normalize
        # ic(df.shape)

        df['r'] = np.sqrt(df['x [mm]']**2 + df['y [mm]']**2)
        F = df.loc[df['r'] <= 150].reset_index(drop=True)
        # ic(F.shape)

        cart = ['x', 'y', 'z']
        angle_list = []
        # fig, axs = plt.subplots(3, 2)
        # for i, a in enumerate(cart):
        # E_u = F_u[f'E{a}Re [V/m]'] + 1j*F_u[f'E{a}Im [V/m]']
        # E = F[f'E{a}Re [V/m]'] + 1j*F[f'E{a}Im [V/m]']
        E_u = df_u[f'|Emag|']
        E = df[f'|Emag|']
        dp = np.arccos(np.abs(np.vdot(E, E_u))/(np.linalg.norm(E, ord=2)*np.linalg.norm(E_u, ord=2)))
        angle_list.append(dp)

        min_x, max_x = min(F_u['x [mm]']), max(F_u['x [mm]'])
        step_size = 1  # mm from CST studio export options
        sh = int((max_x-min_x)/step_size + 1)
        # ic(sh)
        # ic(df_u['x [mm]'].tolist())

        dd = np.reshape(df_u['|Emag|'].tolist(), (sh, sh))
        Z_u = np.ma.array(dd, mask=dd == 0)

        dd = np.reshape(df['|Emag|'].tolist(), (sh, sh))
        Z = np.ma.array(dd, mask=dd == 0)
        ic(freq_reference)
        axs[0, nm - 1].contourf(np.reshape(df_u['x [mm]'].tolist(), (sh, sh)), np.reshape(df_u['y [mm]'].tolist(), (sh, sh)), Z_u, cmap='jet')
        axs[1, nm - 1].contourf(np.reshape(df['x [mm]'].tolist(), (sh, sh)), np.reshape(df['y [mm]'].tolist(), (sh, sh)), Z, cmap='jet')

        ic(np.abs(np.vdot(E, E_u))/(np.linalg.norm(E, ord=2)*np.linalg.norm(E_u, ord=2)))
        ic(angle_list)
        ic('')


    plt.tight_layout()
    plt.show()
