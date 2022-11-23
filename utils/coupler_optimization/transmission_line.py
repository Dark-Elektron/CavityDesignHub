import numpy as np
import matplotlib.pyplot as plt

epsilon_0 = 8.854e-12
c = 299792458


def test(Z, C, Z1, L1):
    # this function computes the total resistance
    N = 3000
    freq = np.linspace(10, 1600, N) * 1e6
    Ztot = []

    for p in range(0, N):
        w = 2 * np.pi * freq[p]
        # Za = parallel(Z, 1 / (1j * w * C))
        Za = Z + 1 / (1j * w * C)
        Za = trans_line(Z1, Za, L1, w)
        Ztot.append(Za)

    return Ztot, freq


def hook_Z_total(Z, Z1, Z2, Z3, l1, l2, l3, Ln, Cn, Ct, C2t, M):
    # this function computes the total resistance
    N = 3000
    freq = np.linspace(100, 600, N) * 1e6
    Ztot = []

    for p in range(0, N):
        w = 2 * np.pi * freq[p]
        # parallel function gives the overall resistance for parallel resistances
        Za = parallel(Z, 1 / (1j * w * Ct))
        # trans_line function computes the resistance at the input of a
        # transmission line given (Z0, ZL, length, ang. freq)
        Zb = trans_line(Z3, Za, l3, w)
        Zc = 1 / (1j * w * C2t) + Zb
        Zd = trans_line(Z2, Zc, l2, w)
        Ze = parallel(1j * w * M, Zd)
        Zf = trans_line(Z1, Ze, l1, w)
        Zg = 1j * (w * Ln - 1 / (Cn * w)) + Zf
        Ztot.append(Zg)

    return Ztot, freq


def trans_line(Z0, ZL, l, w):
    freq = w / (2 * np.pi)
    lamda = 3e8 / freq
    beta = 2 * np.pi / lamda
    Zin = Z0 * (ZL + 1j * Z0 * np.tan(beta * l)) / (Z0 + 1j * ZL * np.tan(beta * l))

    return Zin


def parallel(Z1, Z2):
    Zpar = Z1 * Z2 / (Z1 + Z2)
    return Zpar


def parallel_plate_capacitor(r=None, d=None, C=None, epsilon_r=1):
    # r is radius
    # d is separation
    # C is capacitance
    # epsilon_r is relative permittivity
    if r is None:
        A = C*d/(epsilon_r*epsilon_0)
        r = np.sqrt(A/np.pi)
        return r

    elif d is None:
        A = np.pi * r**2
        d = epsilon_r*epsilon_0*A/C
        return d

    elif C is None:
        A = np.pi * r**2
        C = epsilon_r*epsilon_0*A/d
        return C
    else:
        print("Please enter at least 2 parameter values.")
        return 0


def coaxial_plate_capacitor(r_in=None, r_out=None, L=None, C=None, epsilon_r=1):
    if r_in is None:
        r_in = r_out/np.e**(2*np.pi*epsilon_0*epsilon_r*L/C)
        return r_in

    elif r_out is None:
        r_out = r_in*np.e**(2*np.pi*epsilon_0*epsilon_r*L/C)
        return r_out

    elif L is None:
        L = C*(np.log(r_out/r_in))/(2*np.pi*epsilon_0*epsilon_r)
        return L

    elif C is None:
        C = 2*np.pi*epsilon_0*epsilon_r*L/(np.log(r_out/r_in))
        return C
    else:
        print("Please enter at least 3 parameter values.")
        return 0


def excentric_plate_capacitor(rin, rout, separation, L=None, C=None, epsilon_r=1):
    delta = rout-(separation+rin)

    if C is None:
        C = 4*np.pi*epsilon_0*epsilon_r*L/(2*np.log((rout**2+rin**2-delta**2 + np.sqrt((rout**2+rin**2-delta**2)**2-4*rin**2*rout**2))/(2*rin*rout)))
        return C

    elif L is None:
        L = C * (2.0*np.log((rout**2+rin**2-delta**2 + np.sqrt((rout**2+rin**2-delta**2)**2-4*rin**2*rout**2))/(2*rin*rout))) / (4*np.pi*epsilon_0*epsilon_r)
        return L


if __name__ == '__main__':
    # # [Z, Z1, Z2, Z3, l1_temp, l2, l3_temp, Ln, Cn, Ct_temp, C2t, M, value1, value2, freq[peak1], freq[peak2], peak1, peak2, f1_error, f2_error]
    # # 0: [50, 96, 112, 153.3, 0.2327101638143539, 0.01, 0.07632997086791919, 5e-08, 3.15381710938313e-12, 1.998279485080013e-11, 2.012479605796781e-12, 4.083986593346089e-09, 0.0, -0.12749545498529158, 490963654.5515171, 518672890.9636546, 1909, 2186, 2853654.551517129, 1877109.0363454223]
    # # 1: [50, 96, 112, 153.3, 0.23376279539330125, 0.01, 0.07632997086791919, 5e-08, 3.15381710938313e-12, 2.0824900113958024e-11, 2.012479605796781e-12, 4.083986593346089e-09, 0.0, -0.3081677276980792, 489763254.41813934, 517672557.5191731, 1897, 2176, 1653254.4181393385, 2877442.4808269143]
    # Z_total, f = hook_Z_total(50, 96, 112, 153.3, 0.2327101638143539, 0.01, 0.07632997086791919, 5e-08, 3.15381710938313e-12, 1.998279485080013e-11, 2.012479605796781e-12, 4.083986593346089e-09)
    # Y_dB = 20 * np.log10(np.real(1 / np.array(Z_total)) / max(np.real(1 / np.array(Z_total))))
    # # plt.plot(f, Y_dB)
    # # plt.show()
    #
    # # I choose to calculate radius so my inputs are the required capacitance and the spacing between the plates: C2t, d
    # C2t, d = 2.012e-12, 3e-3
    # r = parallel_plate_capacitor(C=C2t, d=d)*1e3
    # print(f"rh4: {r*2} :: diameter")
    #
    # # Here, I choose to calculate the inner radius so I provide outer radius, length and required capacitance: r_out, L, Ct
    # Ct, r_out, L, epsilon_r = 19.983e-12, 21.2e-3, 18e-3, 10
    # r_in = coaxial_plate_capacitor(r_out=r_out, L=L, C=Ct, epsilon_r=epsilon_r)*1e3
    # print(f"r_in: {r_in*2} :: diameter")
    #
    # # # Here, I choose to calculate the inner radius so I provide outer radius, length and required capacitance: r_out, L, Ct
    # Z1 = 96 # ohm
    # C1 = 1/(c*Z1)
    # # r_out, l1, epsilon_r = 21.2e-3, 18e-3, 10
    # # r_in = coaxial_plate_capacitor(r_out=r_out, L=L, C=C1, epsilon_r=epsilon_r)*1e3
    # # print(f"r_in: {r_in*2} :: diameter, C1: {C1}")
    #
    # # Here, I choose to calculate the inner radius so I provide outer radius, length and required capacitance: r_out, L, Ct
    # Z2 = 112 # ohm
    # C2 = 1/(c*Z2)
    # r_h6_2, l2, epsilon_r = (103/2)*1e-3, 10e-3, 1
    # r_h4 = coaxial_plate_capacitor(r_out=r_h6_2, L=l2, C=C2, epsilon_r=epsilon_r)*1e3
    # print(f"r_h4: {r_h4*2} :: diameter, C2: {C2}")
    #
    # # Here, I choose to calculate the inner radius so I provide outer radius, length and required capacitance: r_out, L, Ct
    # Z3 = 153.3 # ohm
    # C3 = 1/(c*Z3)
    # r_h6_2, l3, epsilon_r = (103/2)*1e-3, 76.33e-3, 1
    # r_h1 = coaxial_plate_capacitor(r_out=r_h6_2, L=l3, C=C3, epsilon_r=epsilon_r)*1e3
    # print(f"r_h1: {r_h1*2} :: diameter, C3: {C3}")
    #
    # ############################
    # # calculate for C from l
    # C2 = coaxial_plate_capacitor(r_in=18/2, r_out=103/2, L=10)
    # Z2 = 1/(c*C2)
    # print(Z2)
    #
    # # calculate for C from l
    # C3 = coaxial_plate_capacitor(r_in=10/2, r_out=103/2, L=75)
    # Z3 = 1/(c*C3)
    # print(Z3)

    ######################################################################
    # Test

    # calculate C
    C = coaxial_plate_capacitor(r_in=16.2e-3, r_out=42.4e-3, L=50e-3)
    print(C)
    Z_total, freq = test(50, 0.09e-12, 113, 125e-3)
    Y_dB = 20 * np.log10(np.real(1 / np.array(Z_total)) / max(np.real(1 / np.array(Z_total))))

    plt.grid()
    plt.plot(freq, Y_dB)
    plt.show()
