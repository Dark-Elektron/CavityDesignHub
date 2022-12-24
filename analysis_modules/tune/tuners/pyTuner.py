import json

from scipy.optimize import fsolve
from termcolor import colored
from analysis_modules.eigenmode.SLANS.slans_geom_par import SLANSGeometry
from utils.file_reader import FileReader
import numpy as np
from itertools import groupby

slans_geom = SLANSGeometry()
fr = FileReader()
file_color = 'cyan'


DEBUG = True


def print_(*arg):
    if DEBUG:
        print(colored(f'\t\t\t{arg}', file_color))


class PyTune:
    def __init__(self):
        self.plot = None

    def tuneL(self, par_mid, par_end, target_freq, beampipes, bc,
              parentDir, projectDir, iter_set, proc=0, conv_list=None):
        # tv => tune variable
        indx = 5

        if proc == '':
            fid = '_process_0'
        else:
            fid = f'_process_{proc}'

        # get parameters
        freq_list = []
        tv_list = []
        error = 1

        slans_geom.cavity(1, 1, par_mid, par_end, par_mid, f_shift=0, bc=bc,
                          beampipes=beampipes, proc=proc, fid=fid,
                          parentDir=parentDir, projectDir=projectDir)
        dirc = fr'{projectDir}\SimulationData\SLANS\Cavity{fid}\cavity_{bc}.svl'
        d = fr.svl_reader(dirc)

        tv = par_end[indx]
        freq = d['FREQUENCY'][0]
        freq_list.append(freq)
        tv_list.append(tv)

        # first shot
        tv = tv + 2

        par_end[indx] = tv

        tol = iter_set[1]
        # max_iter = iter_set[2]
        n = 0

        f_min, f_max = 0, 0
        # tv_min, tv_max = 0, 0
        # control_switch = True

        while abs(error) > tol:
            # print(par_mid, par_end)
            # run slans cavity code
            slans_geom.cavity(1, 1, par_mid, par_end, par_mid, f_shift=0, bc=bc,
                              beampipes=beampipes, proc=proc, fid=fid,
                              parentDir=parentDir, projectDir=projectDir)

            # get results and compare with set value
            d = fr.svl_reader(dirc)
            freq = d['FREQUENCY'][0]
            freq_list.append(freq)
            tv_list.append(tv)
            # update bounds
            if f_min < freq < target_freq and freq > f_min:
                f_min = freq
                # tv_min = tv

            if target_freq < freq < f_max and freq < f_max:
                f_max = freq
                # tv_max = tv

            # update bounds
            # if self.f_min < freq <= target_freq and freq > self.f_min:
            #     self.f_min = freq
            #     self.tv_min = tv
            #
            # if freq >= target_freq:
            #     if self.f_max == 0:
            #         self.f_max = freq
            #         self.tv_max = tv
            #     elif freq <= self.f_max:
            #         self.f_max = freq
            #         self.tv_max = tv

            # calculate slope of line from base point to new point

            if freq_list[n] - freq_list[n + 1] != 0:
                m = (tv_list[n + 1] - tv_list[n]) / (freq_list[n + 1] - freq_list[n])

                if iter_set[0] == "Linear Interpolation":
                    # calculate new Req with straight line formula
                    step = m * (target_freq - freq_list[n+1])
                    if step > 10:
                        step = 10
                    if step < -10:
                        step = -10

                    tv = tv_list[n+1] + step

            # if freq_list[n] - freq_list[n+1] != 0 and control_switch and self.f_max - self.f_min != 0:
            #     if iter_set[0] == "Linear Interpolation":
            #         m = (tv_list[n+1] - tv_list[n])/(freq_list[n+1] - freq_list[n])
            #         # calculate new tv with straight line formula
            #         tv = tv_list[n+1] + m*(target_freq - freq_list[n+1])

                    # # check if new tv is between the current min and max tv
                    # if self.tv_min < abs(tv) < self.tv_max or (self.tv_min == 0 or self.tv_max == 0):
                    #     if self.tv_min == 0:
                    #         tv = self.tv_max - 5
                    #     pass
                    # else:
                    #     if self.tv_min*self.tv_max != 0:
                    #         tv = self.tv_min +
                #         (self.tv_max-self.tv_min)*(target_freq-self.f_min)/(self.f_max - self.f_min)
                    #
                    #     elif self.tv_min != 0:
                    #         # check the side of the curve the new value falls to
                    #         if abs(tv) > self.tv_min:
                    #             tv = max(self.tv_min + 1, self.tv_min*(target_freq/self.f_min))
                    #         else:
                    #             tv = max(self.tv_min - 1, self.tv_min/(target_freq/self.f_min))
                    #     elif self.tv_max != 0:
                    #         # check the side of the curve the new value falls to
                    #         if abs(tv) < self.tv_max:
                    #             tv = self.tv_max - 1
                    #         else:
                    #             tv = self.tv_max + 1

                else:
                    # newton interpolation
                    a_s = self.divided_diff(freq_list, tv_list)[0, :]

                    # use curve to compute for x
                    tv = self.newton_poly(a_s, freq_list, target_freq)

            # elif self.f_max - self.f_min != 0:
            #     # switch control to just this part of the code for the end parts of the simulation for convergence
            #     control_switch = False
            #     if self.tv_min*self.tv_max != 0:
            #         tv = (self.tv_min + self.tv_max)/2

            # print(f'{proc} L', tv_list)
            # print(f'{proc} freq:', freq_list)

            # change tv
            par_end[indx] = tv

            # if equal, break else continue with new shape
            error = target_freq - freq_list[-1]

            if n > 20:
                print('here at breaking')
                break

            # condition for repeated last four values
            if self.all_equal(freq_list[-2:]):
                print("Break because last two elements are equal")
                break

            if tv_list[-1] < 0:
                print("Negative value encountered. It is possible that there no solution for the parameter input set.")
                break

            # check if Req is less than lower limit
            if tv < par_end[0]:
                print("Break because L is less than A")
                break

            # check if alpha is less than or greater than 90.5
            if par_mid == par_end:
                alpha, error_msg = self.calculate_alpha(par_mid[0], par_mid[1], par_mid[2],
                                                        par_mid[3], par_mid[4], tv, par_mid[6], 0)
                if alpha < 90.0 or error_msg != 1:
                    break
            else:
                alpha, error_msg = self.calculate_alpha(par_end[0], par_end[1], par_end[2],
                                                        par_end[3], par_end[4], tv, par_mid[6], 0)
                if alpha < 90.0 or error_msg != 1:
                    break

                alpha, error_msg = self.calculate_alpha(par_mid[0], par_mid[1], par_mid[2],
                                                        par_mid[3], par_mid[4], tv, par_mid[6], 0)
                if alpha < 90.0 or error_msg != 1:
                    break

            # update convergence list
            if conv_list:
                conv_list[proc] = [tv_list, freq_list]
            # self.write_output(tv_list, freq_list, fid, projectDir)

            n += 1

        # return best answer from iteration
        min_error = [abs(x-target_freq) for x in freq_list]
        key = min_error.index(min(min_error))
        print(freq_list, tv_list)
        import matplotlib.pyplot as plt
        plt.scatter(tv_list, freq_list)
        plt.show()

        return tv_list[key], freq_list[key]

    def tuneR(self, par_mid, par_end, target_freq, beampipes, bc, parentDir, projectDir,
              iter_set, proc=0, conv_list=None):

        if proc == '':
            fid = '_process_0'
        else:
            fid = f'_process_{proc}'

        # get parameters
        freq_list = []
        Req_list = []
        error = 1
        slans_geom.cavity(1, 1, par_mid, par_end, par_mid, f_shift=0, bc=bc, beampipes=beampipes, proc=proc, fid=fid,
                          parentDir=parentDir, projectDir=projectDir)
        dirc = fr'{projectDir}\SimulationData\SLANS\Cavity{fid}\cavity_{bc}.svl'
        d = fr.svl_reader(dirc)

        Req = par_end[6]
        freq = d['FREQUENCY'][0]
        freq_list.append(freq)
        Req_list.append(Req)

        # first shot
        Req = Req + 5
        par_end[6], par_mid[6] = Req, Req

        # iteration settings
        tol = iter_set[1]
        max_iter = iter_set[2]
        n = 0
        # Req_min, Req_max = 0, 0
        f_min, f_max = 0, 0

        while abs(error) > tol:
            # run slans cavity code
            slans_geom.cavity(1, 1, par_mid, par_end, par_mid, f_shift=0, bc=bc, beampipes=beampipes, proc=proc,
                              fid=fid, parentDir=parentDir, projectDir=projectDir)

            # get results and compare with set value
            d = fr.svl_reader(dirc)
            freq = d['FREQUENCY'][0]
            freq_list.append(freq)
            Req_list.append(Req)

            # update bounds
            if f_min < freq < target_freq and freq > f_min:
                f_min = freq
                # Req_min = Req

            if target_freq < freq < f_max and freq < f_max:
                f_max = freq
                # Req_max = Req

            # calculate slope of line from base point to new point
            if freq_list[n] - freq_list[n + 1] != 0:
                m = (Req_list[n + 1] - Req_list[n]) / (freq_list[n + 1] - freq_list[n])

                if iter_set[0] == "Linear Interpolation":
                    # calculate new Req with straight line formula
                    Req = Req_list[n+1] + m * (target_freq - freq_list[n+1])

                    # # check if new Req is between the current min and max Req
                    # if self.Req_min < abs(Req) < self.Req_max:
                    #     pass
                    # else:
                    #     if self.Req_min*self.Req_max != 0:
                    #         Req = (self.Req_min + self.Req_max)/2
                    #
                    #     elif self.Req_min != 0:
                    #         # check the side of the curve the new value falls to
                    #         if abs(Req) > self.Req_min:
                    #             Req = max(self.Req_min + 1, self.Req_min*(target_freq/self.f_min))
                    #         else:
                    #             Req = max(self.Req_min - 1, self.Req_min/(target_freq/self.f_min))
                    #     elif self.Req_max != 0:
                    #         # check the side of the curve the new value falls to
                    #         if abs(Req) < self.Req_max:
                    #             Req = self.Req_max - 1
                    #         else:
                    #             Req = self.Req_max + 1

                else:
                    # newton interpolation
                    a_s = self.divided_diff(freq_list, Req_list)[0, :]

                    # use curve to compute for x
                    Req = self.newton_poly(a_s, freq_list, target_freq)
            # else:
            #     Req = Req_list[n+1] + 0.5

            # change R
            par_end[6], par_mid[6] = Req, Req

            # if equal, break else continue with new shape
            error = target_freq - freq

            # avoid infinite loop-1
            if n > max_iter:
                break

            # # check if Req is less than lower limit
            # if Req < par_mid[1] + par_mid[3] + par_mid[4] or Req < par_end[1] + par_end[3] + par_end[4]:
            #     break

            # check if alpha is less or greater than 90.5
            if par_mid == par_end:
                alpha, error_msg = self.calculate_alpha(par_mid[0], par_mid[1], par_mid[2],
                                                        par_mid[3], par_mid[4], par_mid[5], Req, 0)
                if alpha < 90.0 or error_msg != 1:
                    break
            else:
                alpha, error_msg = self.calculate_alpha(par_end[0], par_end[1], par_end[2],
                                                        par_end[3], par_end[4], par_end[5], Req, 0)
                if alpha < 90.0 or error_msg != 1:
                    break

                alpha, error_msg = self.calculate_alpha(par_mid[0], par_mid[1], par_mid[2],
                                                        par_mid[3], par_mid[4], par_mid[5], Req, 0)
                if alpha < 90.0 or error_msg != 1:
                    break

            # update convergence list
            if conv_list:
                conv_list[proc] = [Req_list, freq_list]

            # print('Req', Req_list)
            # print('freq:', freq_list)

            n += 1

        # return best answer from iteration
        min_error = [abs(x-target_freq) for x in freq_list]
        key = min_error.index(min(min_error))

        # import matplotlib.pyplot as plt
        # plt.scatter(Req_list, freq_list)
        # plt.show()

        return Req_list[key], freq_list[key]

    @staticmethod
    def divided_diff(x, y):
        """
        function to calculate the divided
        differences table
        """
        n = len(y)
        coef = np.zeros([n, n])
        # the first column is y
        coef[:, 0] = y

        for j in range(1, n):
            for i in range(n - j):
                coef[i][j] = (coef[i + 1][j - 1] - coef[i][j - 1]) / (x[i + j] - x[i]) if (x[i + j] - x[i]) != 0 else 0

        return coef

    @staticmethod
    def newton_poly(coef, x_data, x):
        """
        evaluate the newton polynomial
        at x
        """
        n = len(x_data) - 1
        p = coef[n]
        for k in range(1, n + 1):
            p = coef[n - k] + (x - x_data[n - k]) * p
        return p

        # return Shape([self.A, self.B, self.a, self.b, self.r, self.l, self.R])

    @staticmethod
    def all_equal(iterable):
        g = groupby(iterable)
        return next(g, True) and not next(g, False)

    @staticmethod
    def write_output(tv_list, freq_list, fid, projectDir):
        dd = {"tv": tv_list, "freq": freq_list}

        with open(fr"{projectDir}\SimulationData\SLANS\Cavity{fid}\convergence_output.json", "w") as outfile:
            json.dump(dd, outfile, indent=4, separators=(',', ': '))

    def calculate_alpha(self, A, B, a, b, Ri, L, Req, L_bp):

        data = ([0 + L_bp, Ri + b, L + L_bp, Req - B],
                [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        df = fsolve(self.ellipse_tangent,
                    np.array([0.5*a + L_bp, Ri + 0.5 * b, L - A + L_bp, Req - 0.5 * B]), args=data, full_output=True)
        x1, y1, x2, y2 = df[0]
        error_msg = df[-2]

        m = (y2 - y1) / (x2 - x1)
        alpha = 180 - np.arctan(m) * 180 / np.pi
        return alpha, error_msg

    @staticmethod
    def ellipse_tangent(z, *data):
        coord, dim = data
        h, k, p, q = coord
        a, b, A, B = dim
        x1, y1, x2, y2 = z

        f1 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
        f2 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
        f3 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
        f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

        return f1, f2, f3, f4
