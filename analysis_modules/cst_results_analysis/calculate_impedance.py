import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mplcursors
from icecream import ic

c0 = 299792458
mu0 = 4 * np.pi * 1e-7
imp_cal = 1


def rrc_filter(signal):
    import numpy as np
    import matplotlib.pyplot as plt

    # Create a time signal (for example, a simple sine wave)
    t = np.linspace(0, 1, 1000)  # Time from 0 to 1 second
    original_signal = np.sin(2 * np.pi * 5 * t)  # A sine wave with frequency 5 Hz
    t = signal[0]*1e-3/299792458
    original_signal = signal[1]

    # Define the raised cosine filter
    T_symbol = 1e-10 # Symbol duration in seconds
    roll_off = 1  # Roll-off factor (adjust as needed)
    t_filter = np.linspace(-5 * T_symbol, 5 * T_symbol, 1000)  # Time for the filter
    t_filter = np.arange(1000) - (1000 - 1) // 2
    raised_cosine_filter = 1 / T_symbol * np.sinc(t_filter / T_symbol) * (np.cos(np.pi * roll_off * t_filter / T_symbol)/(1-(2 * roll_off * t_filter / T_symbol)**2))
    # raised_cosine_filter[t_filter == 0] = np.pi / (4 * T_symbol)  # Set value at zero to avoid division by zero

    # Normalize the filter to maintain the signal amplitude
    normalized_raised_cosine_filter = raised_cosine_filter / np.sum(raised_cosine_filter)

    # Apply the filter using convolution
    filtered_signal = np.convolve(original_signal, normalized_raised_cosine_filter, mode='same')

    # # Plot the original and filtered signals
    # plt.figure(figsize=(10, 6))
    # plt.plot(t, original_signal, label='Original Signal')
    # plt.plot(t, filtered_signal, label='Filtered Signal', ls='-.')
    # plt.title('Original Signal and Raised Cosine Filtered Signal')
    # plt.xlabel('Time')
    # plt.ylabel('Amplitude')
    # plt.legend()
    # plt.grid(True)
    # plt.show()

    # Compute the Fourier transform of the original and filtered signals
    original_signal_fft = np.fft.fft(original_signal)
    filtered_signal_fft = np.fft.fft(filtered_signal)

    # Frequency values for the Fourier transform
    frequency = np.fft.fftfreq(len(t), t[1] - t[0])

    # Plot only the positive part of the frequency spectra
    positive_freq_mask = (frequency >= 0) & (frequency <= 2047e6)

    plt.figure(figsize=(10, 6))

    plt.subplot(2, 1, 1)
    plt.plot(frequency[positive_freq_mask], np.abs(original_signal_fft[positive_freq_mask]))
    plt.title('Positive Frequencies - Fourier Transform of Original Signal')
    plt.xlabel('Frequency')
    plt.ylabel('Amplitude')
    plt.yscale('log')
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.plot(frequency[positive_freq_mask], np.abs(filtered_signal_fft[positive_freq_mask]))
    plt.title('Positive Frequencies - Fourier Transform of Filtered Signal')
    plt.xlabel('Frequency')
    plt.ylabel('Amplitude')
    plt.yscale('log')
    plt.grid(True)

    plt.tight_layout()
    plt.show()


def fft_wakepotential():
    pass

def get_var(var, assembly):
    result = []
    for nn, res_folder in enumerate(assembly):
        result.append(pd.read_csv(fr"{res_folder}\{var}", sep='\t', header=None))

    result = pd.concat(result, ignore_index=True)
    return result


def get_field(folder, var):
    return pd.read_csv(fr"{folder}\{var}", sep='\t', header=None)


if __name__ == '__main__':

    # signal = pd.read_csv(fr"D:\CST Studio\PhD\tt\Wakefield\W_C3795_3DQW_1FPC\Export\Particle Beams_ParticleBeam1_Wake potential_Z.txt", sep='\t', header=None)
    # print(signal[1])
    # rrc_filter(signal)

    # rf0 = r'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\E_C3794\Export'
    rf0 = r'D:\CST Studio\5. tt\Eigenmode\E_C3795\Export'
    rf1 = r'D:\CST Studio\5. tt\Eigenmode\E_C3795_from_1433MHz\Export'
    rf2 = r'D:\CST Studio\5. tt\Eigenmode\E_C3795_from_1700MHz\Export'
    rf3 = r'D:\CST Studio\5. tt\Eigenmode\E_C3795_from_2100MHz\Export'
    rf4 = r'D:\CST Studio\5. tt\Eigenmode\E_C3795_from_2407MHz\Export'

    rf01hc = r'D:\CST Studio\5. tt\Eigenmode\E_C3795_4HC\Export'
    rf02hc = r'D:\CST Studio\5. tt\Eigenmode\E_C3795_4HC_from_1493MHz\Export'
    rf03hc = r'D:\CST Studio\5. tt\Eigenmode\E_C3795_4HC_from_1700MHz\Export'
    rf04hc = r'D:\CST Studio\5. tt\Eigenmode\E_C3795_4HC_from_2005MHz\Export'

    rf379401hc = r'D:\CST Studio\3. W\Assembly\E_3794_4HC_1FPC_Optimized\Export'
    rf379402hc = r'D:\CST Studio\3. W\Assembly\E_3794_4HC_1FPC_Optimized_from_1000MHz\Export'
    rf379403hc = r'D:{\CST Studio\3. W\Assembly\E_3794_4HC_1FPC_Optimized_from_1250MHz\Export'
    rf379404hc = r'D:\CST Studio\3. W\Assembly\E_3794_4HC_1FPC_Optimized_from_1425MHz\Export'
    rf379405hc = r'D:\CST Studio\3. W\Assembly\E_3794_4HC_1FPC_Optimized_from_1597MHz\Export'
    rf379406hc = r'D:\CST Studio\3. W\Assembly\E_3794_4HC_1FPC_Optimized_from_1697MHz\Export'
    rf379407hc = r'D:\CST Studio\3. W\Assembly\E_3794_4HC_1FPC_Optimized_from_1834MHz\Export'
    rf379408hc = r'D:\CST Studio\3. W\Assembly\E_3794_4HC_1FPC_Optimized_from_2000MHz\Export'

    folder = r'D:\CST Studio\PhD\tt\Module\Eigenmode'
    rf01dqw = r'D:\CST Studio\PhD\tt\Module\Eigenmode\E_C3795_2DQW_1FPC\Export'
    rf02dqw = r'D:\CST Studio\PhD\tt\Module\Eigenmode\E_C3795_2DQW_1FPC_D2\Export'
    rf03dqw = r'D:\CST Studio\PhD\tt\Module\Eigenmode\E_C3795_2DQW_1FPC_D2_2\Export'

    folder = r'D:\CST Studio\PhD\tt\Assembly'
    rf012qrwg1 = r'D:\CST Studio\PhD\tt\Assembly\E_C3794_2QRWG_1FPC\Export'
    rf012qrwg2 = r'D:\CST Studio\PhD\tt\Assembly\E_C3794_2QRWG_1FPC_1250_1450\Export'
    rf012qrwg3 = r'D:\CST Studio\PhD\tt\Assembly\E_C3794_2QRWG_1FPC_1450_1700\Export'
    rf012qrwg4 = r'D:\CST Studio\PhD\tt\Assembly\E_C3794_2QRWG_1FPC_1700_2000\Export'

     # plot all
    # names = ["C3795_4DQW", "C3795_4HC", "C3794_4HC"]
    # assemblies = [[rf0, rf1, rf2, rf3, rf4], [rf01hc, rf02hc, rf03hc, rf04hc],
    #               [rf379401hc, rf379402hc, rf379403hc, rf379404hc, rf379405hc, rf379406hc, rf379407hc, rf379408hc]]
    # n_modes_list = [[37, 100, 89, 100, 200], [35, 100, 90, 100], [100, 100, 100, 100, 100, 150, 150, 100]]
    # names = ["C3795_4DQW", "C3795_4HC", "C3794_4HC"]
    # assemblies = [[rf0, rf1], [rf01hc, rf02hc],
    #               [rf379401hc, rf379402hc, rf379403hc, rf379404hc]]
    # n_modes_list = [[37, 100], [35, 100], [100, 100, 100, 100]]

    #  # plot all
    # names = ["C3795_4DQW", "C3795_4HC"]
    # assemblies = [[rf0, rf1, rf2, rf3, rf4], [rf01hc, rf02hc, rf03hc, rf04hc]]
    # n_modes_list = [[37, 100, 89, 100, 200], [35, 100, 90, 100]]

    #  # plot all
    # names = ["C3795_4DQW"]
    # assemblies = [[rf0, rf3]]
    # n_modes_list = [[37, 100]]

    #  # plot all
    # names = ["C3795_4HC"]
    # assemblies = [[rf02hc]]
    # n_modes_list = [[100]]
    #  # plot all
    # names = ["C3795_4HC"]
    # assemblies = [[rf01hc]]
    # n_modes_list = [[35]]

    # # Plot c3794
    # names = ["C3794_4HC"]
    # assemblies = [[rf379401hc, rf379402hc, rf379403hc, rf379404hc, rf379405hc, rf379406hc, rf379407hc, rf379408hc]]
    # n_modes_list = [[100, 100, 100, 100, 100, 150, 150, 100]]

    # # Plot c3794
    # names = ["E_C3794_2QRWG_1FPC_1700_2000"]
    # assemblies = [[rf01dqw, rf02dqw, rf03dqw]]
    # n_modes_list = [[132, 100, 74]]
    #
    # Plot c3794
    names = ["E_C3794_2QRWG_1FPC_1700_2000"]
    assemblies = [[rf012qrwg1, rf012qrwg2, rf012qrwg3, rf012qrwg4]]
    n_modes_list = [[26, 8, 23, 52]]

    # rf04hc_mon_pass = r'D:\CST Studio\5. tt\Eigenmode\E_C3795_long_passband\Export'
    # names = ['C3795']
    # assemblies = [[rf04hc_mon_pass]]
    # n_modes_list = [[44]]

    for a, assembly in enumerate(assemblies):
        # assembly = [rf0]
        # result_folders = [rf0, rf1, rf2, rf3, rf4]
        modes = n_modes_list[a]
        freq_all = get_var("Frequency (Multiple Modes).txt", assembly)
        # ic(freq_all)

        mode_num = np.shape(freq_all)[0]
        freq_all = freq_all.iloc[:, 1] * 1e6  # second column is frequency: 1e6 for MHz and 1e9 for GHz

        V_z, VT_x, VT_y = [], [], []
        n = 0
        for i, freq in enumerate(freq_all):
            if i < sum(modes[0:n+1]):
                Ex = get_field(assembly[n], fr"e_X (Z)_Mode {i - sum(modes[0:n]) + 1}.txt")
                Ey = get_field(assembly[n], fr"e_Y (Z)_Mode {i - sum(modes[0:n]) + 1}.txt")
                Ez = get_field(assembly[n], fr"e_z (Z)_Mode {i - sum(modes[0:n]) + 1}.txt")
                Hx = get_field(assembly[n], fr"h_X (Z)_Mode {i - sum(modes[0:n]) + 1}.txt")
                Hy = get_field(assembly[n], fr"h_Y (Z)_Mode {i - sum(modes[0:n]) + 1}.txt")
                Hz = get_field(assembly[n], fr"h_z (Z)_Mode {i - sum(modes[0:n]) + 1}.txt")

                Z_axis = Ex.iloc[:, 0] * 1e-3  # 1e-3 to convert mm to m
                v_z = np.trapz((Ez.iloc[:, 0] + 1j * Ez.iloc[:, 1]) * np.exp(1j * 2 * np.pi * freq * Z_axis / c0), Z_axis)

                vt_x = np.trapz((Ex.iloc[:, 0] + 1j * Ex.iloc[:, 1] - c0 * mu0 * (Hy.iloc[:, 0] + 1j * Hy.iloc[:, 1]))
                                * np.exp(1j * 2 * np.pi * freq * Z_axis / c0), Z_axis)
                vt_y = np.trapz((Ey.iloc[:, 0] + 1j * Ey.iloc[:, 1] + c0 * mu0 * (Hx.iloc[:, 0] + 1j * Hx.iloc[:, 1]))
                                * np.exp(1j * 2 * np.pi * freq * Z_axis / c0), Z_axis)

                V_z.append(v_z)
                VT_x.append(vt_x)
                VT_y.append(vt_y)

                if i == sum(modes[0:n+1]) - 1:
                    n += 1

            # change

        R_Q = abs(np.array(V_z) ** 2) / (2 * np.pi * freq_all)
        R_Q_x = abs(np.array(VT_x) ** 2) / (2 * np.pi * freq_all)
        R_Q_y = abs(np.array(VT_y) ** 2) / (2 * np.pi * freq_all)
        R_Q_T = R_Q_x + R_Q_y
        # ic(R_Q_x)
        # pp = plt.scatter(freq_all, VT_x, label=f'{nn}', zorder=10, edgecolor='k')
        # mplcursors.cursor(pp)

        Qext_all = get_var("Q-Factor (lossy E) (Multiple Modes).txt", assembly)
        Qext_all = Qext_all.iloc[:, 1]
        Z_ = Qext_all * R_Q * 0.5
        Z_T = Qext_all * R_Q_T * 0.5 * 2 * np.pi * freq_all / c0
        Z_export = [freq_all, R_Q_T, Qext_all, Z_T]

        Z_x = Qext_all * R_Q_x * 0.5 * 2 * np.pi * freq_all / c0
        Z_y = Qext_all * R_Q_y * 0.5 * 2 * np.pi * freq_all / c0

        # fs = np.linspace(min(freq_all), max(freq_all), 1000)  # the frequency interval to calculate impedance from eigenmodes
        fs = np.linspace(1, max(freq_all)*1.5, 10000)  # the frequency interval to calculate impedance from eigenmodes
        # fs = freq_all

        Z_sum_z, Z_sum_x, Z_sum_y = [], [], []
        for i, freq in enumerate(freq_all):
            z_sum_x = 1 * Z_x[i] / (1 + 1j * Qext_all[i] * (fs / freq - freq / fs))
            z_sum_y = 1 * Z_y[i] / (1 + 1j * Qext_all[i] * (fs / freq - freq / fs))
            z_sum_z = 1 * Z_[i] / (1 + 1j * Qext_all[i] * (fs / freq - freq / fs))

            Z_sum_x.append(z_sum_x)
            Z_sum_y.append(z_sum_y)
            Z_sum_z.append(z_sum_z)

        Z = sum(Z_sum_z)*1e-3  # impedance in x direction, kOhm
        Z_x_all = sum(Z_sum_x)*1e-3  # impedance in x direction, kOhm/m
        Z_y_all = sum(Z_sum_y)*1e-3  # impedance y direction, kOhm/m
        ZT = 0.5*(Z_x_all + Z_y_all)  # impedance y direction

        po1 = plt.plot(fs, abs(Z), label=fr'$Z_\parallel$ {names[a]}')
        po2 = plt.plot(fs, abs(ZT), label=fr'$Z_\perp$ {names[a]}')

        # get cst calculated values
        Z_CST = get_var("Z_kOhm.txt", assembly).iloc[:, 1]
        ZT_CST = get_var("Z_T_kOhm_m.txt", assembly).iloc[:, 1]

        po3 = plt.scatter(freq_all, Z_*1e-3, label=fr'$Z_\parallel (CST)$ {names[a]}', edgecolors='k', alpha=0.5)
        po4 = plt.scatter(freq_all, Z_T*1e-3, label=fr'$Z_\perp (CST)$ {names[a]}', edgecolors='k', alpha=0.5)

        po3 = plt.scatter(freq_all, Z_CST, label=fr'$Z_\parallel (CST)$ {names[a]}', marker='x')
        po4 = plt.scatter(freq_all, ZT_CST, label=fr'$Z_\perp (CST)$ {names[a]}', marker='+')

        # save to excel
        data = np.array([fs*1e-6, abs(Z), abs(ZT)]).T
        # data = np.array([freq_all*1e-6, abs(Z_*1e-3), abs(Z_T*1e-3)]).T
        df = pd.DataFrame(data, columns=['f [MHz]', 'Z [kOhm]', 'ZT [kOhm/m]'])
        # ic(df)
        df.to_excel(fr"{folder}\Z_ZT {names[a]}.xlsx")

        mplcursors.cursor(po1)
        mplcursors.cursor(po2)
        mplcursors.cursor(po3)
        mplcursors.cursor(po4)
        # plt.scatter(freq_all, abs(Z_x))
        # plt.scatter(freq_all, abs(Z_y))
    plt.legend()
    plt.yscale('log')
    plt.show()

