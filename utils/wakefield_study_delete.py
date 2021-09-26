
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

# constants
Tb = 300e-9
Ql = 1e4
omega_n = 400e6 # Hz
omega_n2 = 500e6 # Hz
c = 299792458 # m/s
Vq = 1e6 # V

t_list = np.linspace(0, 200 * Tb, 10000)  # 2 seconds
w_list = np.linspace(0, 2e9, 1000)  # 2 seconds

Vb = [0* a for a in Vq * np.cos(omega_n * t_list) * np.exp(-t_list / Tb)]
for n, w in enumerate(w_list):
    Vb2 = Vq *np.exp(-n*Tb*w/(2*Ql)) * (np.cos(n*w*Tb)*np.cos(w * t_list) - np.sin(n*w*Tb)*np.sin(w * t_list))

    Vb = [a+b for a, b in zip(Vb, Vb2)]

plt.plot(t_list, Vb, 'k')
plt.show()


# Perform Fourier transform using scipy
from scipy import fftpack
y_fft = fftpack.fft(Vb)

Fs = 100 # sampling rate
n = np.size(t_list)
fr = Fs/2 * np.linspace(0, 1, int(n/2))
y_m = 2/n * abs(y_fft[0:np.size(fr)])

plt.plot(fr, y_m)
plt.show()