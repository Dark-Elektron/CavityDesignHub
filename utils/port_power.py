from scipy.fft import fft, ifft, fftfreq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':

    # initialize constants
    c0 = 299792458
    mu0 = 4 * np.pi * 1e-7

    # import signals
    result_folder = r'D:\CST Studio\Comparison FPC_Coupler Configs\WL_C3794_HOM_Power_calc\Export'

    # iex = pd.read_csv(fr"{result_folder}\Frequency (Multiple Modes).txt", sep='\t', header=None)

    # specify port numbers and number of modes in each port
    pn = [1, 2]  # ports number
    mn = [5, 5]  # number of modes [int]

    # solve for weighting transfer function
    iex = pd.read_csv(fr"{result_folder}\Particle Beams_ParticleBeam1_Current.txt", sep='\t', header=None)
    Ifcc = pd.read_csv(fr"{result_folder}\Particle Beams_ParticleBeam1_Charge distribution spectrum.txt", sep='\t', header=None)  # beam spectrum

    Fpm = {}
    bpm = {}
    for i, p in enumerate(pn):
        for m in range(1, mn[i] + 1):
            b = pd.read_csv(fr"{result_folder}\Port signals_o{p}({m}),pb.txt", sep='\t', header=None)

            # 6.15
            # sample spacing
            T = b[0][1]
            print(T)

            N = len(list(b[1]))
            xf = fftfreq(N, T)[:N//2]
            Fpm[f"{i}{m}"] = np.divide(fft(list(b[1])), fft(list(iex[1])))  #
            plt.plot(xf, np.abs(Fpm[f"{i}{m}"][0:N//2]))
            plt.show()
            # print(list(b[1]))
            # print(Fpm[f"{i}{m}"])
            # bpm = np.multiply(Fpm[f"{i}{m}"], Ifcc[1])

    # from scipy.fft import fft, fftfreq
    # # Number of sample points
    # N = 600
    # # sample spacing
    # T = 1.0 / 9000.0
    # x = np.linspace(0.0, N*T, N, endpoint=False)
    # y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
    # plt.plot(x, y)
    # plt.show()
    # yf = fft(y)
    # xf = fftfreq(N, T)[:N//2]
    # import matplotlib.pyplot as plt
    # plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
    # plt.grid()
    # plt.show()