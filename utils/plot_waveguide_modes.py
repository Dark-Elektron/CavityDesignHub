# beampipe dimensions
import numpy as np
from icecream import ic
from matplotlib import pyplot as plt
from scipy.special import *
import cmath
import matplotlib as mpl

m0 = 9.1093879e-31
q0 = 1.6021773e-19
mu0 = 4 * np.pi * 1e-7
eps0 = 8.85418782e-12
c0 = 2.99792458e8
eta0 = 376.7303134111465


def plot_settings():
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12

    mpl.rcParams['axes.labelsize'] = 14
    mpl.rcParams['axes.titlesize'] = 14
    mpl.rcParams['legend.fontsize'] = 12
    mpl.rcParams['legend.title_fontsize'] = 12

    mpl.rcParams['figure.dpi'] = 100

    # Set the desired colormap
    plt.rcParams['axes.prop_cycle'] = plt.cycler('color', plt.cm.Set2.colors)


def plot_analytical():
    pass


def calc_te_analytical(m, n, theta, radius):
    j_mn = jn_zeros(m, n)[n - 1]
    j_mn_p = jnp_zeros(m, n)[n - 1]
    A = 1
    kc = j_mn_p/R
    w = kc/np.sqrt(mu0*eps0)
    k = 0
    beta = cmath.sqrt(k**2 - kc**2)

    Er = -1j*w*mu0*m/(kc**2*radius) * A * (np.cos(m*theta) - np.sin(m*theta))*jv(m, kc*radius)
    Et = 1j*w*mu0/kc * A * (np.sin(m*theta) + np.cos(m*theta))*jvp(m, kc*radius)
    Ez = 0
    Hr = -1j*beta/kc * A * (np.sin(m*theta) + np.cos(m*theta))*jvp(m, kc*radius)
    Ht = -1j*beta*m/(kc**2*radius) * A * (np.cos(m*theta) - np.sin(m*theta))*jv(m, kc*radius)
    Hz = A * (np.sin(m*theta) + np.cos(m*theta))*jv(m, kc*radius)

    Emag = np.abs(np.sqrt(Er**2 + Et**2 + Ez**2))
    Hmag = np.abs(np.sqrt(Hr ** 2 + Ht ** 2 + Hz ** 2))

    return Emag, Hmag


def calc_tm_analytical(m, n, theta, radius, component='abs'):
    j_mn = jn_zeros(m, n)[n - 1]
    j_mn_p = jnp_zeros(m, n)[n - 1]
    A = 1
    k = 0  # no propagation in z
    kc = j_mn/R
    w = kc/np.sqrt(mu0*eps0)
    beta = cmath.sqrt(k**2 - kc**2)

    Ez = A * (np.sin(m*theta) + np.cos(m*theta))*jv(m, kc*radius)
    Er = -1j*beta/kc * A * (np.sin(m*theta) + np.cos(m*theta))*jvp(m, kc*radius)
    Et = -1j*beta*m/(kc**2*radius) * A * (np.cos(m*theta) - np.sin(m*theta))*jv(m, kc*radius)
    Hr = 1j*w*eps0*m/(kc**2*radius) * A * (np.cos(m*theta) - np.sin(m*theta))*jv(m, kc*radius)
    Ht = -1j*w*eps0/kc * A * (np.sin(m*theta) + np.cos(m*theta))*jvp(m, kc*radius)
    Hz = 0

    Emag = np.abs(np.sqrt(Er**2 + Et**2 + Ez**2))
    Hmag = np.abs(np.sqrt(Hr ** 2 + Ht ** 2 + Hz ** 2))

    if component.lower() == 'abs':
        return Emag, Hmag
    if component.lower() == 'azimuthal':
        return Et, Ht
    if component.lower() == 'radial':
        return Er, Hr
    if component.lower() == 'longitudinal':
        return Ez, Hz


plot_settings()
mpl.rcParams['figure.figsize'] = [5, 5]

R = 80e-3
m, n = [1, 1]

# create radial grid
r = np.linspace(1e-6, R, 500)
theta = np.linspace(0, 2 * np.pi, 500)
radius_matrix, theta_matrix = np.meshgrid(r, theta)
ax = plt.subplot(111, polar=True)
# ax.plot(theta_matrix, radius_matrix, color='r', ls='none', marker='.')
Emag, Hmag = calc_te_analytical(m, n, theta_matrix, radius_matrix)
ax.contour(theta_matrix, radius_matrix, Emag, levels=30, cmap='RdBu')

x_tr = 20e-3
ax.plot(np.linspace(0, 2 * np.pi, 500), x_tr*np.ones(500), c='k')

ax.plot([0, 0], [0, x_tr], [0, np.pi/(2*m)], [0, x_tr], c='r', marker='o', mfc='none', lw=2, ms=10)
ax.plot([0, np.pi/4], [0, x_tr], [0, np.pi/4+np.pi/(2*m)], [0, x_tr], c='b', ls='--', marker='+', lw=2, ms=10)

# get E field at certa
ic(calc_te_analytical(m, n, np.pi/15, x_tr))
ic(calc_te_analytical(m, n, np.pi/15 + np.pi/(2*m), x_tr))
ic(calc_te_analytical(m, n, np.pi/3, x_tr))
ic(calc_te_analytical(m, n, np.pi/4+np.pi/(2*m), x_tr))

# ax.contour(theta_matrix, radius_matrix, Hmag, linestyles='--', cmap='gray')
# ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.show()
