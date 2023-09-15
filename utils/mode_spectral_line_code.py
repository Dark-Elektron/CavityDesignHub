import numpy as np
from scipy.special import *
import matplotlib.pyplot as plt

if __name__ == '__main__':

    Ri = 0.08038
    L = 4*0.0935
    c = 299792458

    # permutations up to 4
    p = []

    M = [0, 1, 2, 3, 4, 5, 6]
    N = [1, 2, 3, 4, 4, 5, 6]
    P = [0, 1, 2, 3, 4, 5, 6]
    mode_list = []
    mode_name_list = ['TM', 'TE']
    f_list = []

    for m in M:
        for n in N:
            # get bessel
            j_mn = jn_zeros(m, n)[n-1]
            j_mn_p = jnp_zeros(m, n)[n-1]
            J = [j_mn, j_mn_p]
            for p in P:
                for mode_type, j in enumerate(J):
                    # formula
                    f = c/(2*np.pi)*((j/Ri)**2 + (p*np.pi/L)**2)**0.5

                    # apend mode name
                    mode_list.append(f'{mode_name_list[mode_type]}{m}{n}{p}')

                    # append to f list
                    f_list.append(f*1e-9)

    mode_list_sorted = [x for _, x in sorted(zip(f_list, mode_list))]
    f_list_sorted = sorted(f_list)

    print(mode_list_sorted)
    print(f_list_sorted)

    color = ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c']
    step = [0.02, 0.2, 0.4, 0.6, 0.8]
    count = 0
    for mode, f in zip(mode_list_sorted, f_list_sorted):
        plt.axvline(x=f, c=color[count%len(color)])
        plt.text(f-0.02, step[count%len(step)], f'{mode}', rotation=90, c=color[count%len(color)], clip_on=True)
        count += 1

    plt.xlim(min(f_list_sorted)-min(f_list_sorted)/5, 3)
    plt.show()
