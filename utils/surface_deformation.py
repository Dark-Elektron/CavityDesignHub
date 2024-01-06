import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from icecream import ic


def random(n=11, sigma=0.1):
    sampl = np.random.uniform(low=0.0005, high=0.0015, size=(n,))
    return sampl


def random_deform(surface=None):
    if surface is None:
        n = 5
        bottom = np.array([np.linspace(0, 1, n), [0 for i in range(n)]])
        right = np.array([[1 for i in range(n)], np.linspace(0, 1, n)])
        top = np.array([np.flip(np.linspace(0, 1, n)), [1 for i in range(n)]])
        left = np.array([[0 for i in range(n)], np.flip(np.linspace(0, 1, n))])

        surface = np.vstack((bottom.T, right.T, top.T, left.T))

    plt.plot(surface[:, 0], surface[:, 1])
    # sigma = np.arange(1, 50, 10)
    sigma = [50]
    for s in sigma:
        deform_matrix = np.diag(random(len(surface), s)) + np.identity(len(surface))
        # deform surface
        surface_def = deform_matrix@surface
        plt.plot(surface_def[:, 0], surface_def[:, 1], label=s)

    plt.legend()
    # plt.gca().set_aspect('equal')
    plt.show()


def gauss(n=11,sigma=1, shift=0):
    r = np.linspace(-int(n/2)+0.5,int(n/2)-0.5, n)
    g = 1 / (sigma * np.sqrt(2*np.pi)) * np.exp(-(r-shift)**2/(2*sigma**2))
    return g/max(g)

def gaussian_deform(surface=None):
    # plt.plot(gauss(len(surface), 50))
    # plt.show()
    if surface is None:
        n = 5
        bottom = np.array([np.linspace(0, 1, n), [0 for i in range(n)]])
        right = np.array([[1 for i in range(n)], np.linspace(0, 1, n)])
        top = np.array([np.flip(np.linspace(0, 1, n)), [1 for i in range(n)]])
        left = np.array([[0 for i in range(n)], np.flip(np.linspace(0, 1, n))])
    
        surface = np.vstack((bottom.T, right.T, top.T, left.T))

    plt.plot(surface[:, 0], surface[:, 1], label='undeformed')
    # sigma = np.arange(1, 50, 10)
    sigma = [int(0.1*len(surface)), int(0.2*len(surface))]
    shift = [-250, 250]
    max_disp = [-0.02, 0.02]

    for s in sigma:
        ic(np.max(gauss(len(surface), s, shift=-0)), np.min(gauss(len(surface), s, shift=-0)))
        # deform_matrix = -0.02*np.diag(gauss(len(surface), s, shift=-0)) + np.identity(len(surface))
        deform_vector = np.atleast_2d(-1e-3*gauss(len(surface), s, shift=-0)).T
        ic(np.max(deform_vector), np.min(deform_vector))
        # deform surface
        # surface_def = deform_matrix@surface
        surface_def = deform_vector + surface
        plt.plot(surface_def[:, 0], surface_def[:, 1], label='gaussian')  #, marker='o', mfc='none'

    plt.legend()
    # plt.gca().set_aspect('equal')
    plt.show()


cav_surf = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP\CEPCv2_n2\contour.txt", sep=r'\s+', header=None, skiprows=3)[[1, 0]]
cav_surf_A_pert = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP\contour_A_perturbed.txt", sep=r'\s+', header=None, skiprows=3)[[1, 0]]
# plt.plot(cav_surf_A_pert.drop_duplicates().to_numpy()[:, 0], cav_surf_A_pert.drop_duplicates().to_numpy()[:, 1], label='A deformed')
# cav_surf = pd.read_csv(fr"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS\C3794_n2\contour.txt", sep=r'\s+', header=None, skiprows=3)[[1, 0]]
print(cav_surf.to_numpy())
gaussian_deform(cav_surf.drop_duplicates().to_numpy())
# random_deform(cav_surf.to_numpy()*1e3)

# gg = gauss(1000, 100, shift=100)
# plt.plot(gg/max(gg))
# plt.show()