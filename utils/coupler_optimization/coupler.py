import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps
from icecream import ic
from sympy import symbols, solve
from tqdm import tqdm
import mplcursors
from transmission_line import parallel_plate_capacitor, coaxial_plate_capacitor, excentric_plate_capacitor


class Coupler:
    def __init__(self):
        pass

    def LHCHookType(self, f_thresh):
        # This is the main script that runs all the functions needed
        # Functions needed in order to run this script:
        # hook_design.m
        # hook_Z_total.m
        # parallel.m
        # trans_line.m

        # Initial values see fig. 1 of report (url: goo.gl/z3V7Ld )###
        d_tube, dc = 103e-3, 42.4e-3  # m, Fixed value
        d1, d2, d3 = 16e-3, 16e-3, 8e-3  # m, Fixed diameter values (choice)
        Z1 = 96  # To be calculated via CST Studio because of the excentricity of the transmission line with the tube
        Z2 = self.trans_line_impedance(d2, d_tube)
        Z3 = self.trans_line_impedance(d3, d_tube)

        # f_FM, f1, f2 = 400.79e6, 488e6, 520e6  # MHz, FM freq, 1st and 2nd dipole resonance freq.
        f_FM, f1, f2 = 400.79e6, 518.89e6, 531.98e6  # MHz, FM freq, 1st and 2nd dipole resonance freq.
        print(Z2, Z3, np.sqrt(f1 * f2))

        Ln = 50e-9  # H, has to calculated from CST Studio. Related to the insertion depth
        Zt = 112  # Ohm, this line resistance is used in Z_total calculation. see Frank Gerigk's thesis

        Z = 50  # Ohm, terminating impedance
        Q = 10
        c0 = 299792458
        l2 = 10e-3  # Fixed, determined experimentally

        nk = 2
        dk_max = 0.05
        dk = 0.028  # k correction factor

        # terminating transformar
        lZ, dZ = 30e-3, 18.4e-3

        # hook_design function - uses the procedure in Gerigks thesis to give theoretical results
        # this results serve as a basis in order to begin optimization
        l1, l3, tm, Cn, Ct, C2t, M, L1n, C1n, L2, C2 = self.hook_design(f_FM, f1, f2, Ln, Z1, Z2, Z3, Zt, Z, Q, c0, l2, dk)
        print(l1, l3, Cn, Ct, C2t, M)

        # self.sweep_hook(f1, f2, Z, Z1, Z2, Z3, l1, l2, l3, Ln, Cn, Ct, C2t, M)

        # Z_total_E, freq = self.hook_Z_total_electric_couping(Z, Z1, Z2, Z3, l1, l2, l3, Ln, Cn, Ct, C2t, M)
        Z_total_M, freq = self.hook_Z_total_magnetic_couping(Z, Z1, Z2, Z3, l1, l2, l3, Ln, Cn, Ct, C2t, M)
        #

        # circuit parameters
        circuit_params = [Ln, Cn, L1n, C1n, M, L2, C2, C2t, Ct, Z]

        # transmission line parameters
        trans_params = [Ln, Cn, l1, Z1, M, l2, Z2, C2t, l3, Z3, Ct, Z]

        # convert to 3D model parameters
        self.circuit_2_3d(d1, d2, d3, d_tube, dc, l1, l2, l3, tm, C2t, Ct, Cn, Ln, lZ, dZ)

        #
        # conductance in dB (vector)
        # Y_dB = 20 * np.log10(np.real(1 / np.array(Z_total_E)) / max(np.real(1 / np.array(Z_total_E))))
        Y_dB_magnetic = 20 * np.log10(np.real(np.array(Z_total_M)) / max(np.real(np.array(Z_total_M))))

        # plt.plot(freq, Y_dB)
        plt.plot(freq, Y_dB_magnetic)
        # self.test()
        plt.show()

    def test(self):
        # circuit definition
        Z0 = 50

        # capacitor
        C = 11.21e-12
        Ct = [C, 'C', 'p']  # inductor parallel

        # tline
        zt3, lt3 = 153.3, 7.61e-2
        Zt3 = [zt3, lt3, 't']  # transmission line

        # capacitor
        C = 2e-12
        C2t = [C, 'C', 's']  # inductor parallel

        # tline
        zt2, lt2 = 112, 1e-2
        Zt2 = [zt2, lt2, 't']  # transmission line

        # inductor
        L = 4.9e-9
        M = [L, 'L', 'p']  # inductor parallel

        # tline
        zt1, lt1 = 96, 21.85e-2
        Zt1 = [zt1, lt1, 't']  # transmission line

        # inductor
        L = 50e-9
        Ln = [L, 'L', 's']  # inductor parallel

        # capacitor
        C = 3.15e-12
        Cn = [C, 'C', 's']  # inductor parallel

        hook_coupler_circuit = [Ct, Zt3, C2t, Zt2, M, Zt1, Ln, Cn]

        # # ####################### test circuit
        # Z0 = 1e-12
        # do = 103
        # di = 16
        # zt1, lt1 = 60*np.log(do/di), 200e-3
        # Zt1 = [zt1, lt1, 't']
        #
        # L = 50e-9
        # L1 = [L, 'L', 'p']  # inductor parallel
        # Zt2 = [zt1, 145*1e-3, 't']
        #
        # # loop capacitor and inductor
        # d_loop, d_wire = 40e-3, 16e-3
        # L = 4*np.pi*1e-7*(d_loop/2)*(np.log(8*d_loop/d_wire)-2)
        # Ln = [L, 'L', 's']
        # print("Ln", Ln)
        #
        # f_rejected = 400.79e6 # Hz
        # C = 1/L*(1/(2*np.pi*f_rejected))**2
        # print("Cn", C)
        # Cn = [C, 'C', 's']
        #
        # # calculate corresponding length
        # rout, rin, separation = 103e-3/2, 16e-3/2, 3e-3
        # L = excentric_plate_capacitor(rin, rout, separation, C=C)
        # print("L", L)
        #
        # r = parallel_plate_capacitor(C=C, d=5e-2)
        # print(r)
        #
        # test_circuit = [Zt1, [Cn, Ln]]
        #
        # ######################################################
        # hook_coupler_circuit = [Ct, Zt3, C2t, Zt2, M, Zt1, [Ln, Cn]]

        # transmission curve circuit
        Z_tot, freq = self.test_Z_total(Z0, hook_coupler_circuit)
        Y_dB = 20 * np.log10(np.real(1 / np.array(Z_tot)) / max(np.real(1 / np.array(Z_tot))))

        mplcursors.cursor(plt.plot(freq, Y_dB))
        plt.show()

    def test_Z_total(self, Z0, circuit):
        # this function computes the total impedance
        N = 1001
        freq = np.linspace(1, 800, N) * 1e6
        Ztot = []

        for f in freq:
            w = 2 * np.pi * f

            # solve for total Z in a loop
            zz = Z0
            for element in circuit:
                # check if part of sub-chain
                if isinstance(element[0], list):
                    # process sub-chain first
                    z_sc = 0
                    for el in element:
                        z_sc = self.chain(z_sc, el, w)

                    zz = self.parallel(zz, z_sc)
                else:
                    zz = self.chain(zz, element, w)

            Ztot.append(zz)

        return np.array(Ztot), freq

    def chain(self, zz, element, w):
        if element[-1] == 't':
            zz = self.trans_line(element[0], zz, element[1], w)
        elif element[-1] == 'p':
            if element[-2] == 'L':
                zz = self.parallel(1j * w * element[0], zz)
            elif element[-2] == 'C':
                zz = self.parallel(1 / (1j * w * element[0]), zz)
            elif element[-2] == 'R':
                zz = self.parallel(element[0], zz)
            else:
                print("Unrecognized circuit element")
        else:
            if element[-2] == 'L':
                zz = 1j * w * element[0] + zz
            elif element[-2] == 'C':
                zz = 1 / (1j * w * element[0]) + zz
            elif element[-2] == 'R':
                zz = element[0] + zz
            else:
                print("Unrecognized circuit element")

        return zz

    @staticmethod
    def trans_line(Z0, ZL, l, w):
        freq = w / (2 * np.pi)
        lamda = 3e8 / freq
        beta = 2 * np.pi / lamda

        if ZL == np.infty:  #Open circuit
            Zin = -1j*Z0/np.tan(beta*l)
        elif ZL == 0:  # Short circuit
            Zin = 1j*Z0*np.tan(beta*l)
        else:
            Zin = Z0 * (ZL + 1j * Z0 * np.tan(beta * l)) / (Z0 + 1j * ZL * np.tan(beta * l))

        return Zin

    def sweep_hook(self, f1, f2, Z, Z1, Z2, Z3, l1, l2, l3, Ln, Cn, Ct, C2t, M):
        dB_thresh = 1  # a threshold in dB error (concerning the 7 dB difference)
        f_thresh = 3e6  # a threshold in frequency position of the resonance
        dB_difference = 7  # the difference in dB between the two resonances (the Ytotal difference in dB)

        # creation of Ct values matrix for the optimization procedure
        nCt = 20  # number of different values to be checked
        dCt = 16e-12  # range of the values that will be checked around the calculated value
        delta_Ct = np.linspace(-dCt / 2, dCt / 2, nCt)

        # creation of l1 values matrix for the optimization procedure
        nl1 = 20  # number of different values to be checked
        dl1 = 2e-2  # range of the values that will be checked around the calculated value
        delta_l1 = np.linspace(-dl1 / 2, dl1 / 2, nl1)

        #  creation of l3 values matrix for the optimization procedure
        nl3 = 20  # number of different values to be checked
        dl3 = 2e-2  # range of the values that will be checked around the calculated value
        delta_l3 = np.linspace(-dl3 / 2, dl3 / 2, nl3)

        # initiation of matrixes in order to save the generated data

        max_solutions_given = 30
        # Results = np.zeros(max_solutions_given, 12)

        counter = 0
        solution = 0

        ## denp.ping on the width of the values for Ct l1 and l3 and the thresholds there might be more than one results
        # that will satisfy the criteria that were given. This value is used in
        # order to put a maximum on the amount of this results - if there are too
        # many results this means that the threshold and the values used are not
        # working very well
        max_solutions_given = 30

        remaining_minutes = 0

        ## here, the optimization begins. It is a brut force type so there are three for loops
        ## each one for every variable that is being examined (Ct, l1, l3)
        peaks1 = []
        peaks2 = []
        freq1 = []
        freq2 = []
        values1 = []
        values2 = []
        cr_freq = []
        cr_Y = []
        result = []

        count = 0

        for i in tqdm(range(0, nCt), desc='outer'):
            for j in range(0, nl1):
                for p in range(0, nl3):
                    # this is just for the Loading window
                    counter = counter + 1

                    # this temp (temporary) variables are used in order to put a
                    # small change every time in the actual value (actual_value + Delta_actual_value)
                    Ct_temp = Ct + delta_Ct[i]
                    l1_temp = l1 + delta_l1[j]
                    l3_temp = l3 + delta_l3[p]

                    # hook_Z_total function computes the total complex resistance
                    # for different freqs (vector)
                    # of the equivalent circuit (it is just resistances in series and in parallel)
                    Z_total, freq = self.hook_Z_total_electric_couping(Z, Z1, Z2, Z3, l1_temp, l2, l3_temp, Ln, Cn,
                                                                       Ct_temp, C2t, M)

                    # conductance in dB (vector)
                    Y_dB = 20 * np.log10(np.real(1 / np.array(Z_total)) / max(np.real(1 / np.array(Z_total))))
                    # plt.plot(freq, Y_dB)
                    # plt.show()

                    # plot(Y_dB)
                    # because of the need to have two resonances, findpeaks function is used
                    # in order to make the check if the maxima are on the level and
                    # the frequency that is needed
                    peaks, _ = sps.find_peaks(Y_dB)

                    # if there are no or only one peak, the algorithm skips this
                    # step
                    if len(peaks) <= 1:
                        continue

                    value1 = Y_dB[peaks][0]
                    value2 = Y_dB[peaks][1]
                    peak1 = peaks[0]
                    peak2 = peaks[1]

                    # save peaks for plot later
                    peaks1.append(peak1)
                    peaks2.append(peak2)
                    freq1.append(freq[peak1])
                    freq2.append(freq[peak2])
                    values1.append(value1)
                    values2.append(value2)
                    count = count + 1

                    # check if value of first and second peaks is zero
                    if 0 <= abs(value1) <= 1 and 0 <= abs(value2) <= 1:

                        f1_error = abs(f1 - freq[peak1])
                        f2_error = abs(f2 - freq[peak2])

                        if f1_error < f_thresh and f2_error < f_thresh:
                            print("Found one solution")
                            cr_freq.append(freq)
                            cr_Y.append(Y_dB)

                            result.append([Z, Z1, Z2, Z3, l1_temp, l2, l3_temp, Ln, Cn, Ct_temp, C2t, M, value1, value2,
                                           freq[peak1], freq[peak2], peak1, peak2, f1_error, f2_error])

                    # if there are many solutions that satisfies our criteria (odd situation)
                    # the algorithm breaks the loop (enough is enough :))
                    if max_solutions_given <= solution:
                        break

                    # the errors between the given values and the original ones are
                    # calculated here
                    f1_error = abs(f1 - freq[peak1])
                    f2_error = abs(f2 - freq[peak2])

                    dB_dist = value2 - value1
                    dB_error = abs(dB_dist - dB_difference)

                    # if the sum of the errors exceed the sum of the error
                    # threshold given by the user the algorithm skips this step
                    if f1_error / f_thresh > 1:
                        continue
                    if f2_error / f_thresh > 1:
                        continue
                    if dB_error / dB_thresh > 1:
                        continue

                    # counter given in order to save multiple solutions
                    solution = solution + 1

                    # solutions saving (consult this in order to distinct what is what in the results)
                    # Results(solution,:) = [Ct_temp, l1_temp, l3_temp, i, j, p, value1, value2, freq(peak1), freq(peak2),peak1,peak2]

                    # plt.plot(freq, Y_dB) # if there is a result, it is being plotted

                    # just some calculations for the remaining minutes of the
                    # optimization. Info is being displayed at the command window.
        # print(cr_freq)
        # print(cr_Y)
        for i, r in enumerate(result):
            print(f"{i}: {r}")

        for freq, Y in zip(cr_freq, cr_Y):
            mplcursors.cursor(plt.plot(freq, Y))
        plt.show()

    def LHCProbeType(self):
        # This is the main script that runs all the functions needed
        # same approach as the hook coupler also here (more detailed commenting at hook_model.m)
        # Functions needed in order to run this script:
        # probe_design.m
        # probe_Z_total.m
        # parallel.m
        # trans_line.m

        ### Initial values see fig. 7 and table 3 of report (url: goo.gl/z3V7Ld )###
        f1 = 1615e6
        fmid = 2030e6
        f2 = 2330e6
        l3 = 1.e-2
        l4 = 3.5e-2 / 2
        c0 = 3e8
        f0 = np.sqrt(f1 * f2)
        w0 = 2 * np.pi * f0
        k = (f2 - f1) / f0
        wFM = 801.364e6 * 2 * np.pi
        Z = 50
        Zt = 112

        ## probe_design function - uses the procedure in Gerigks thesis to give theoretical results
        # this results serve as a basis in order to begin optimization
        M23, C3, l1, l2, Cn, Ln = self.probe_design(Z, Zt, w0, l3, l4, wFM, k)

        # threshold in error values - difference from the needed value for freq.
        # resonances and db_height
        f_thresh = 2e6
        dB_thresh = 3

        # number of different values
        N = 4

        # these factors here in the case of the probe will be multiplied by the
        # given value in order to change and different values to be investigated
        # There is no particular reason for the choise of multiplication - just an other approach.
        # I chose it subjectively because it made my investigation easier.
        # the boundaries here are a bit synthesized according to the final results
        # in order to approach the best results. If a new investigation (with
        # different criteria) is needed then someone has to increase N and to
        # change these boundaries
        dl = np.linspace(0.9, 1.2, N)
        dM = np.linspace(1, 1.1, N)
        df1 = np.linspace(0.91, 0.97, N)
        df2 = np.linspace(1.01, 1.05, N)

        results = np.zeros(20, 4)

        counter = 0
        solution = 0

        remaining_minutes = 0
        remaining_secs = 0

        for a in range(1, N):
            for b in range(1, N):
                for c in range(1, N):
                    for d in range(1, N):
                        counter = counter + 1

                        f1 = 1615e6 * df1[a]  # variable under investigation
                        fmid = 2030e6
                        f2 = 2330e6 * df2[b]  # variable under investigation
                        l3 = 1.2e-2
                        l4 = 3.5e-2 / 2
                        c0 = 3e8
                        f0 = np.sqrt(f1 * f2)
                        w0 = 2 * np.pi * f0
                        k = (f2 - f1) / f0
                        wFM = 801.364e6 * 2 * np.pi
                        Z = 50
                        Zt = 112
                        M23, C3, l1, l2, Cn, Ln = self.probe_design(Z, Zt, w0, l3, l4, wFM, k)

                        l2_temp = l2 * dl[c]  # variable under investigation
                        M23_temp = M23 * dM[d]  # variable under investigation

                        [Z_total, freq] = self.probe_Z_total(Z, Zt, l1, l2_temp, l3, l4, Ln, Cn, C3, M23_temp)
                        Z_dB = 20 * np.log10(np.real(Z_total) / max(np.real(Z_total)))

                        # In probe coupler it is not very well defined (as it is in
                        # the hook coupler) how many peaks there are and where - and
                        # also many solutions could be satisfactory. So the
                        # threshold values take under consideration the value of
                        # Z_total in dB at the frequencies we are interested at.
                        if Z_dB(912) < -1:
                            continue
                        if Z_dB(1348) < -2:
                            continue
                        if Z_dB(1664) < -1:
                            continue

                        solution = solution + 1  # counting when there are any solutions

                        results[solution, :] = [M23_temp, l2_temp, c, d]

                        # plt.plot(freq,Z_dB,'Linewidth',1.5,'Color',[np.random('unif',1) random('unif',1) random('unif',1)]) # if there is a result it is plotted. For multiple results there are different colors

    @staticmethod
    def hook_design(f_FM, f1, f2, Ln, Z1, Z2, Z3, Zt, Z, Q, c0, l2, dk):
        # this algorithm follows step-by-step the procedure given in the master
        # thesis of frank gerick - the same notation applies so you can follow it
        # easily

        # base circuit
        f0 = np.sqrt(f1 * f2)
        w0 = 2 * np.pi * f0
        wFM = 2 * np.pi * f_FM
        Df = f2 - f1
        k = Df / f0  # + dk
        k = 0.1  # + dk
        ic(k)
        L1 = Zt * np.pi / (2 * w0 * (1 + k))
        ic(L1)
        M = k * L1
        C1 = 1 / ((L1 + M) * w0 ** 2)
        ic(C1)
        L2 = L1
        C2 = 1 / ((L2 + M) * w0 ** 2)
        C1 = C2
        Ze = (L1 + M) * w0 / Q  # Ze in thesis What actually is Q? Why was 10 chosen?
        ic(Ze)
        # notch filter
        Cn = 1 / (Ln * wFM ** 2)
        L1n = L1 - Ln
        ic(L1n)
        C1n = 1 / (1 / C1 - 1 / Cn)
        ic(C1n)

        # transformer
        t = Z / Ze
        ic(t)
        Ctp = -t / (Z * w0 * np.sqrt(t - 1))
        Ct = np.sqrt(t - 1) / (w0 * Z)
        C2t = (1 / C2 + 1 / Ctp) ** (-1)

        # transmission line circuit from lumped element circuit
        Lr = (1 / w0) * (w0 * Ln - 1 / (w0 * Cn))
        ic(Lr)
        Cr = 1 / (w0 ** 2 * Lr)
        ic(Cr)
        lr = c0 / w0 * np.arctan(-1 / (Zt * w0 * Cr)) + 0.5 * c0 / f0
        lm = c0 / w0 * np.arctan(w0 * M / Zt)
        ic(lm)
        l1 = lr - lm
        Xt = w0 * (M + L2) - 1 / (w0 * C2t)
        lm2 = c0 / w0 * np.arctan(w0 * M / Z2)
        Xh = Zt * np.tan(w0 * (lm2 + l2) / c0) - 1 / (C2t * w0)
        ls = c0 / w0 * np.arctan(Xh / Z3)
        l3 = c0 / w0 * np.arctan(Xt / Z3) - ls

        ic(Ln, M, Ct, Cn, C2t, Z, l1, l2, l3, Z1, Z2, Z3)

        return l1, l3, lm, Cn, Ct, C2t, M, L1n, C1n, L2, C2

    @staticmethod
    def probe_design(R, Zt, w0, l3, l4, wFM, k):
        p = 2
        c0 = 3e8
        f0 = w0 / (2 * np.pi)
        t = (R / Zt) ** 2
        Lt = R / w0 * np.sqrt((1 - t) / t)
        c = Zt / (4 * f0)

        x, y = symbols('x y')
        a, b = solve([x + k * np.sqrt(x * y) - c, y + k * np.sqrt(Lt * y) + k * np.sqrt(y * x) - c], [x, y])

        L3 = np.double(a[p])
        L2 = np.double(b[p])

        M23 = (c - L3) * 1
        M12 = c - M23 - L2
        C3 = 1 / (w0 ** 2 * (L3 + M23))

        lm12 = c0 / w0 * np.arctan(w0 * M12 / Zt)
        lm23 = c0 / w0 * np.arctan(w0 * M23 / Zt)
        l1 = 0.25 * c0 / f0 - lm12
        l2 = 0.5 * c0 / f0 - lm12 - lm23
        Ln = w0 * M12 / (w0 - wFM ** 2 / w0)
        Cn = 1 / (wFM ** 2 * Ln)
        t1 = l1 / c0
        t2 = l2 / c0
        t3 = l3 / c0
        t4 = l4 / c0

        return M23, C3, l1, l2, Cn, Ln

    def hook_Z_total_electric_couping(self, Z, Z1, Z2, Z3, l1, l2, l3, Ln, Cn, Ct, C2t, M):
        # this function computes the total resistance
        N = 3000
        freq = np.linspace(1, 800, N) * 1e6
        Ztot = []

        for p in range(0, N):
            w = 2 * np.pi * freq[p]
            # parallel function gives the overall resistance for parallel resistances
            Za = self.parallel(Z, 1 / (1j * w * Ct))
            # trans_line function computes the resistance at the input of a
            # transmission line given (Z0, ZL, length, ang. freq)
            Zb = self.trans_line(Z3, Za, l3, w)
            Zc = 1 / (1j * w * C2t) + Zb
            Zd = self.trans_line(Z2, Zc, l2, w)
            Ze = self.parallel(1j * w * M, Zd)
            Zf = self.trans_line(Z1, Ze, l1, w)
            Zg = 1j * (w * Ln - 1 / (Cn * w)) + Zf
            Ztot.append(Zg)

        return Ztot, freq

    def hook_Z_total_magnetic_couping(self, Z, Z1, Z2, Z3, l1, l2, l3, Ln, Cn, Ct, C2t, M):
        # this function computes the total resistance
        N = 3000
        freq = np.linspace(300, 600, N) * 1e6
        Ztot = []

        for p in range(0, N):
            w = 2 * np.pi * freq[p]
            # parallel function gives the overall resistance for parallel resistances
            Za = self.parallel(Z, 1 / (1j * w * Ct))
            # trans_line function computes the resistance at the input of a
            # transmission line given (Z0, ZL, length, ang. freq)
            Zb = self.trans_line(Z3, Za, l3, w)
            Zc = 1 / (1j * w * C2t) + Zb
            Zd = self.trans_line(Z2, Zc, l2, w)
            Ze = self.parallel(1j * w * M, Zd)
            Zf = self.trans_line(Z1, Ze, l1, w)
            Zg = self.parallel(Zf, 1j * (w * Ln - 1 / (Cn * w)))
            Ztot.append(Zg)

        return Ztot, freq

    def probe_Z_total_electric_couping(self, Z, Zt, l1, l2, l3, l4, Ln, Cn, C3, M):
        # this function computes the total resistance
        N = 3000
        freq = np.linspace(300, 600, N) * 1e6
        Ztot = []

        for p in range(N):
            w = 2 * np.pi * freq[p]
            Za = self.trans_line(Zt, Z, l4, w)
            Zb = 1 / (1j * w * C3) + Za
            Zc = self.trans_line(Zt, Zb, l3, w)
            Zd = self.parallel(1j * w * M, Zc)
            Ze = self.trans_line(Zt, Zd, l2, w)
            Zf = self.parallel(1j * (w * Ln - 1 / (Cn * w)), Ze)
            Zg = self.trans_line(Zt, Zf, l1, w)
            Ztot.append(Zg)

        return Ztot, freq

    def probe_Z_total_magnetic_couping(self, Z, Zt, l1, l2, l3, l4, Ln, Cn, C3, M):
        # this function computes the total resistance
        N = 3000
        freq = np.linspace(300, 600, N) * 1e6
        Ztot = []

        for p in range(N):
            w = 2 * np.pi * freq[p]
            Za = self.trans_line(Zt, Z, l4, w)
            Zb = 1 / (1j * w * C3) + Za
            Zc = self.trans_line(Zt, Zb, l3, w)
            Zd = self.parallel(1j * w * M, Zc)
            Ze = self.trans_line(Zt, Zd, l2, w)
            Zf = self.parallel(1j * (w * Ln - 1 / (Cn * w)), Ze)
            Zg = self.trans_line(Zt, Zf, l1, w)
            Ztot.append(Zg)

        return Ztot, freq

    def tesla_Z_total_electric_couping(self, Z, Zt, l1, l2, l3, l4, Ln, Cn, C3, M):
        # this function computes the total resistance
        N = 3000
        freq = np.linspace(300, 600, N) * 1e6
        Ztot = []

        for p in range(N):
            w = 2 * np.pi * freq[p]
            Za = self.trans_line(Zt, Z, l4, w)
            Zb = 1 / (1j * w * C3) + Za
            Zc = self.trans_line(Zt, Zb, l3, w)
            Zd = self.parallel(1j * w * M, Zc)
            Ze = self.trans_line(Zt, Zd, l2, w)
            Zf = self.parallel(1j * (w * Ln - 1 / (Cn * w)), Ze)
            Zg = self.trans_line(Zt, Zf, l1, w)
            Ztot.append(Zg)

        return Ztot, freq

    def tesla_Z_total_magnetic_couping(self, Z, Zt, l1, l2, l3, l4, Ln, Cn, C3, M):
        # this function computes the total resistance
        N = 3000
        freq = np.linspace(300, 600, N) * 1e6
        Ztot = []

        for p in range(N):
            w = 2 * np.pi * freq[p]
            Za = self.trans_line(Zt, Z, l4, w)
            Zb = 1 / (1j * w * C3) + Za
            Zc = self.trans_line(Zt, Zb, l3, w)
            Zd = self.parallel(1j * w * M, Zc)
            Ze = self.trans_line(Zt, Zd, l2, w)
            Zf = self.parallel(1j * (w * Ln - 1 / (Cn * w)), Ze)
            Zg = self.trans_line(Zt, Zf, l1, w)
            Ztot.append(Zg)

        return Ztot, freq

    @staticmethod
    def parallel(Z1, Z2):
        Zpar = Z1 * Z2 / (Z1 + Z2)
        return Zpar

    @staticmethod
    def trans_line_impedance(din, dout):
        Z = 60 * np.log(dout / din)
        return Z

    def circuit_2_3d(self, d1, d2, d3, d_tube, dc, l1, l2, l3, tm, C2t, Ct, Cn, Ln, lZ, dZ):
        # Calculate radius from the required capacitance and the spacing between the plates: C2t, d
        g2t = 3e-3  # separation, fixed by choice
        r2t = parallel_plate_capacitor(C=C2t, d=g2t) * 1e3
        dcap = 2 * r2t

        # mutual inductance
        alpha_m = 120  # fixed

        # Calculate the inner radius from outer radius, length and required capacitance: r_out, L, Ct
        # Ct, r_out, L, epsilon_r = 19.983e-12, 21.2e-3, 18e-3, 10
        # rt = coaxial_plate_capacitor(r_out=r_out, L=L, C=Ct, epsilon_r=epsilon_r)*1e3

        ln = 69.4
        sn = 2

        # Inductance couplint
        gn = 5  # fixed by choice
        dn = 62  # fixed by choice

        # calculate diameter of expansion from Ct, length of expansion varied, other variables fixed
        # Ct = 2*np.pi/ln(dtube_top/dexp)*(n*er*dd + (l4 - n*dd))
        # n is number of rings, limited to 20
        # dd is thickness of each ring, dd = 1.1mm, d_exp is fixed by choice
        n, dd, er, dt = 9, 1.1, 1, 30
        e0 = 8.854e-12
        lt = Ct / (2 * np.pi * e0 / np.log(dc / dt)) + n * dd * (1 - er)
        tc = n*dd

        ic(d_tube, d1, d2, d3, l1, l2, l3, ln, dn, gn, tm, alpha_m, dcap, g2t, lt, dt, tc, dc, lZ, dZ)

        # self.translate_to_shahnam(d1, d2, d3, d_tube, l1, l2, l3, lm, s2t, d2t, ln, sn)

    @staticmethod
    def translate_to_shahnam(d1, d2, d3, d_tube, l1, l2, l3, lm, s2t, d2t, ln, sn):  # s - capacitor plates separation

        print(f"L1:: rh2: {d1}, rh6: {d_tube}, lh1+lh2: {l1}")
        print(f"L2:: rh2: {d2}, lh3-dh3: {l2}")
        print(f"L3:: rh1: {d3}, lh4: {l3}")
        print(f"M:: dh3: {lm}, alpha_h: {60}")  # update alpha later
        print(f"Ln:: dh1: {60}")  # update later
        print(f"Cn:: dh2: {ln}, ch1: {sn}")
        print(f"C2t:: rh4: {d2t}, ch2: {s2t}, cap_thickness: {6}")
        # print(f"Ct:: rh3: {st}, lh6: {6}")

    def parallel_plate_capacitor(self):
        Z0 = 50

        d_tube, d1, l_tube = 103e-3, 16e-3, 300e-3
        zt = self.trans_line_impedance(d1, d_tube)
        lt = l_tube / 2 - (5 + 3 / 2) * 1e-3
        print(lt)

        # for l in np.linspace(0.05, 1, 20):
        for l in [lt]:
            Zt = [zt, l, 't']  # transmission line
            ic(Zt)

            # capacitor
            dcap, gc = 30e-3, 3e-3
            C = parallel_plate_capacitor(r=dcap/2, d=gc)
            Ct = [C, 'C', 's']  # capacitor series
            ic(Ct)

            circuit = [Zt]
            # transmission curve circuit
            N = 1001
            freq = np.linspace(1, 2000, N) * 1e6
            # this function computes the total resistance
            Z_tot = []
            print(l2, l3, lt)
            for p in range(0, N):
                w = 2 * np.pi * freq[p]
                # parallel function gives the overall resistance for parallel resistances
                Z = self.trans_line(zt, Z0, lt, w)
                Z = 1 / (1j * w * C) + Z
                Z = self.trans_line(zt, Z, lt, w)

                Z_tot.append(Z)

            # Z_tot, freq = self.test_Z_total(Z0, [Zt, Ct, Zt])
            # Z_tot = self.trans_line(zt, np.infty, 2*lt, 2 * np.pi*freq)
            # print(Z_tot)
            Z_tot = np.array(Z_tot)
            # Z_tot = np.sqrt(Z_tot.real**2 + Z_tot.imag**2)

            Y_dB = 20 * np.log10(np.real(1 / np.array(Z_tot)) / max(np.real(1 / np.array(Z_tot))))

            mplcursors.cursor(plt.plot(freq, Y_dB, label=l))
        # mplcursors.cursor(plt.plot(freq, np.sqrt(Z_tot.real**2 + Z_tot.imag**2)))
        plt.legend()
        plt.show()


if __name__ == '__main__':
    coupler = Coupler()
    # coupler.LHCHookType(3e6)
    # coupler.test()

    d_tube, dc = 103e-3, 42.4  # m, Fixed value
    d1, d2, d3 = 16e-3, 16e-3, 8e-3  # m, Fixed diameter values (choice)
    l1, l2, l3 = 0.21437864378362825, 0.01, 0.07602281297625253

    # coupler.circuit_2_3d(d1, d2, d3, d_tube, dc, l1, l2, l3, tm, C2t, Ct, Cn, Ln, lZ, dZ)

    coupler.parallel_plate_capacitor()