# Define the model geometry and materials
import numpy as np
from icecream import ic
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

if __name__ == '__main__':

    L = 1.0  # Length of the waveguide
    N = 100  # Number of grid points
    dx = L / (N - 1)  # Grid spacing
    dy = L / (N - 1)  # Grid spacing
    dz = L / (N - 1)  # Grid spacing

    # Define the source particle properties
    q = 1.0  # Charge of the particle
    m = 1.0  # Mass of the particle
    c = 299792458
    v = c  # Velocity of the particle

    tt = L/v  # ns , total time
    dt = tt/N  # Time step
    ic(dt)

    # Define the electromagnetic properties of the materials
    eps0 = 8.854e-12  # Permittivity of free space
    mu0 = 4 * np.pi * 1e-7  # Permeability of free space
    epsr = 1.0  # Relative permittivity of the waveguide
    mur = 1.0  # Relative permeability of the waveguide

    # Initialize geometry
    x = np.linspace(0, 1, N)
    y = np.linspace(0, 1, N)
    X, Y = np.meshgrid(x, y)

    # Initialize field vectors as 3D array for easier calculation. Location where there are no field values are set to zero
    ex = np.zeros((N, N, N))
    ey = np.zeros((N, N, N))
    ez = np.zeros((N, N, N))
    hx = np.zeros((N, N, N))
    hy = np.zeros((N, N, N))
    hz = np.zeros((N, N, N))
    jx = np.zeros((N, N, N))
    jy = np.zeros((N, N, N))
    jz = np.zeros((N, N, N))

    E = {'0': [ex, ey, ez]}
    H = {'0': [hx, hy, hz]}
    J = {'0': [jx, jy, jz]}

    # Initial conditions
    jx[:, :, 0] = 1
    jy[:, :, 0] = 1
    # ex[<some index>, 0] = ey[<some index>, 0] = [<some vector>] initial field set up by bunch moving at speed of light
    # it has only radial components at the inlet, z=0
    # j[0, :N//2] = 1  # some analytic vector generating function, maybe using mirror charge concept
    # j[0, :N//2] = 1  # some analytic vector generating function, maybe using mirror charge concept
    # unit impulse
    # ex[N//2, N//2] = 1
    # ey[N//2, N//2] = 1
    # plt.imshow(j.T)
    # this should be thw field distribution moved through the structure

    # Initialize the time loop
    tmax = N

    fig, ax = plt.subplots()
    pcolormesh = ax.pcolormesh(np.sqrt(ex[N//2, :, :]**2))
    # fig.colorbar(pcolormesh)


    def animate(t):
        ic(t)
        # update electric and magnetic fields
        ext, eyt, ezt = np.copy(ex), np.copy(ey), np.copy(ez)
        ex[1:-1, 1:-1, 1:-1] = ex[1:-1, 1:-1, 1:-1] \
                               + dt/eps0*(((hz[:, 1:, :] - hz[:, :-1, :])/dy)[1:-1, :-1, 1:-1]
                                          - ((hy[:, :, 1:] - hy[:, :, :-1])/dz)[1:-1, 1:-1, :-1]) \
                               - dt / eps0 * jx[1:-1, 1:-1, 1:-1]

        ey[1:-1, 1:-1, 1:-1] = ey[1:-1, 1:-1, 1:-1] \
                               + dt/eps0*(((hx[:, :, 1:] - hx[:, :, :-1])/dz)[1:-1, 1:-1, :-1]
                                          - ((hz[1:, :, :] - hz[:-1, :, :])/dx)[:-1, 1:-1, 1:-1]) \
                               - dt / eps0 * jy[1:-1, 1:-1, 1:-1]

        ez[1:-1, 1:-1, 1:-1] = ez[1:-1, 1:-1, 1:-1] \
                               + dt/eps0*(((hx[:, 1:, :] - hx[:, :-1, :])/dy)[1:-1, :-1, 1:-1]
                                          - ((hy[1:, :, :] - hy[:-1, :, :])/dx)[:-1, 1:-1, 1:-1]) \
                               - dt / eps0 * jz[1:-1, 1:-1, 1:-1]

        hx[1:-1, 1:-1, 1:-1] = hx[1:-1, 1:-1, 1:-1] \
                               + dt/mu0*(((eyt[:, :, 1:] - eyt[:, :, :-1])/dz)[1:-1, 1:-1, :-1]
                                         - ((ezt[:, 1:, :] - ezt[:, :-1, :])/dy)[1:-1, :-1, 1:-1])

        hy[1:-1, 1:-1, 1:-1] = hz[1:-1, 1:-1, 1:-1] \
                               + dt/mu0*(((ezt[1:, :, :] - ezt[:-1, :, :])/dx)[:-1, 1:-1, 1:-1]
                                         - ((ext[:, :, 1:] - ext[:, :, :-1])/dz)[1:-1, 1:-1, :-1])

        hz[1:-1, 1:-1, 1:-1] = hz[1:-1, 1:-1, 1:-1] \
                               + dt/mu0*(((ext[:, 1:, :] - ext[:, :-1, :])/dy)[1:-1, :-1, 1:-1]
                                         - ((eyt[1:, :, :] - eyt[:-1, :, :])/dx)[:-1, 1:-1, 1:-1])

        # # move source along z
        jx[:, :, t+1] = 1  # some analytic vector generating function, maybe using mirror charge concept
        jy[:, :, t+1] = 1
        jx[:, :, t] = 0
        jy[:, :, t] = 0

        # plt.pcolor(X, Y, np.sqrt(ex**2 + ey**2).T)
        # plt.show()
        # pcolormesh.set_array(np.sqrt(ex[N//2, :, :]**2).flatten())
        pcolormesh.set_array(np.sqrt(ex[N//2, :, :]**2 + ey[N//2, :, :]**2 + ez[N//2, :, :]**2).flatten())
        # pcolormesh.set_array(j.T.flatten())

        return pcolormesh


    ani = FuncAnimation(fig, animate, interval=1, frames=N, repeat=False)
    plt.colorbar(pcolormesh)
    plt.show()