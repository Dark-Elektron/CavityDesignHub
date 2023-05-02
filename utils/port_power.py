from icecream import ic
from scipy.fft import fft, ifft, fftfreq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# initialize constants
from scipy.linalg import dft

c0 = 299792458
mu0 = 4 * np.pi * 1e-7


def beam_spectrum(sigma, fmax, N, q=1e-9):
    sigma = sigma*1e-3  # convert to mm

    f = np.linspace(0, fmax, N)
    w = 2*np.pi*f
    Iw = q*np.exp(-0.5*(sigma*w/c0)**2)

    # normalise current
    Iw = Iw
    print(np.sum(Iw))

    plt.plot(f*1e-6, Iw*1e3)
    plt.show()
    return Iw


if __name__ == '__main__':
    # beam_spectrum(12.1)

    # import signals
    result_folder = r'D:\CST Studio\LONGITIDUNAL_WAKE_TESLA_CAVITY_100m\Export'

    # iex = pd.read_csv(fr"{result_folder}\Frequency (Multiple Modes).txt", sep='\t', header=None)

    # specify port numbers and number of modes in each port
    pn = [1, 2]  # ports number
    mn = [20, 20]  # number of modes [int]

    # solve for weighting transfer function
    iex = pd.read_csv(fr"{result_folder}\Particle Beams_ParticleBeam1_Current.txt", sep='\t', header=None)
    qs = pd.read_csv(fr"{result_folder}\Particle Beams_ParticleBeam1_Charge distribution spectrum.txt", sep='\t', header=None)
    wake_potential_Z = pd.read_csv(fr"{result_folder}\Particle Beams_ParticleBeam1_Wake potential_Z.txt", sep='\t', header=None)
    print("total charge", sum(qs))
    power_outgoing_all_ports = pd.read_csv(fr"{result_folder}\Power_Excitation (pb)_Power Outgoing all Ports.txt", sep='\t', header=None)  # beam spectrum

    # plot impedance
    N = len(list(wake_potential_Z[0]))
    T = iex[0][1]
    xf = fftfreq(N, T)[:N//2]*1e-6
    ic(xf)
    Z = np.abs(fft(list(wake_potential_Z[1])))
    plt.plot(fft(list(iex[1]))[:N//2])
    plt.yscale('log')
    plt.show()

    Fpm = {}
    bpm = {}
    P = {}

    # number of samples
    N = len(list(iex[1]))
    # sample spacing
    T = iex[0][1]
    xf = fftfreq(N, T)[:N//2]

    Ifcc = beam_spectrum(25, max(xf), N//2)

    for i, p in enumerate(pn):
        P[p] = 0
        for m in range(1, mn[i] + 1):
            b = pd.read_csv(fr"{result_folder}\Port signals_o{p}({m}),pb.txt", sep='\t', header=None)

            Fpm[f"{i}{m}"] = np.divide(fft(list(b[1])), fft(list(iex[1])))[:N//2]

            # port signal spectrum is larger than charge distribution spectrum so the minimum of both is taken
            bpm = np.dot(Fpm[f"{i}{m}"], Ifcc)

            P[p] += 2*np.abs(bpm)**2

    ic(P)
    plt.plot(power_outgoing_all_ports[0], power_outgoing_all_ports[1])
    plt.yscale('log')

    plt.show()
