import json

from termcolor import colored
from simulation_codes.SLANS.slans_geom_par import SLANSGeometry
from utils.file_reader import FileReader
import numpy as np
from itertools import groupby

slans_geom = SLANSGeometry()
fr = FileReader()
file_color = 'cyan'


DEBUG = True
def print_(*arg):
    if DEBUG: print(colored(f'\t\t\t{arg}', file_color))


class PyTune:
    def __init__(self):
        self.plot = None

    def tuneL(self, par_mid, par_end, target_freq, beampipes, bc, parentDir, projectDir, iter_set, proc=0):
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

        slans_geom.cavity(1, 1, par_mid, par_end, par_mid, f_shift=0, bc=bc, beampipes=beampipes, proc=proc, fid=fid, parentDir=parentDir, projectDir=projectDir)
        dirc = f'{projectDir}\SimulationData\SLANS\Cavity{fid}\cavity_{bc}.svl'
        d = fr.svl_reader(dirc)

        tv = par_end[indx]
        freq = d['FREQUENCY'][0]
        freq_list.append(freq)
        tv_list.append(tv)

        # first shot
        tv = tv + 10

        par_end[indx], par_mid[indx] = tv, tv

        tol = iter_set[1]
        max_iter = iter_set[2]
        n = 0

        self.f_min, self.f_max = 0, 0
        self.tv_min, self.tv_max = 0, 0
        control_switch = True

        while abs(error) > tol:
            # run slans cavity code
            slans_geom.cavity(1, 1, par_mid, par_end, par_mid, f_shift=0, bc=bc, beampipes=beampipes, proc=proc, fid=fid, parentDir=parentDir, projectDir=projectDir)

            # get results and compare with set value
            d = fr.svl_reader(dirc)
            freq = d['FREQUENCY'][0]
            freq_list.append(freq)
            tv_list.append(tv)
            print(tv_list)
            print(freq_list, self.f_min, self.f_max)

            # update bounds
            if self.f_min < freq <= target_freq and freq > self.f_min:
                self.f_min = freq
                self.tv_min = tv

            if freq >= target_freq:
                if self.f_max == 0:
                    print("First max set")
                    self.f_max = freq
                    self.tv_max = tv
                elif freq <= self.f_max:
                    print("Max update")
                    self.f_max = freq
                    self.tv_max = tv

            # calculate slope of line from base point to new point
            if freq_list[n] - freq_list[n+1] != 0 and control_switch and self.f_max - self.f_min != 0:
                if iter_set[0] == "Linear Interpolation":
                    print("\t\t\t\t 1")
                    m = (tv_list[n+1] - tv_list[n])/(freq_list[n+1] - freq_list[n])
                    # calculate new tv with straight line formula
                    tv = tv_list[n+1] + m*(target_freq - freq_list[n+1])

                    # check if new tv is between the current min and max tv
                    if self.tv_min < abs(tv) < self.tv_max or (self.tv_min == 0 or self.tv_max == 0):
                        print("\t\t\t\t 2")
                        if self.tv_min == 0:
                            print("\t\t\t\t 2a")
                            tv = self.tv_max - 5
                            print("tv, max tv", tv, self.tv_max)
                        pass
                    else:
                        print("\t\t\t\t 3")
                        if self.tv_min*self.tv_max != 0:
                            print("\t\t\t\t  4")
                            tv = self.tv_min + (self.tv_max-self.tv_min)*(target_freq-self.f_min)/(self.f_max - self.f_min)

                        elif self.tv_min != 0:
                            print("\t\t\t\t 5")
                            # check the side of the curve the new value falls to
                            if abs(tv) > self.tv_min:
                                print("\t\t\t\t 6")
                                tv = max(self.tv_min + 1, self.tv_min*(target_freq/self.f_min))
                            else:
                                print("\t\t\t\t 7")
                                tv = max(self.tv_min - 1, self.tv_min/(target_freq/self.f_min))
                        elif self.tv_max != 0:
                            print("\t\t\t\t 8")
                            # check the side of the curve the new value falls to
                            if abs(tv) < self.tv_max:
                                print("\t\t\t\t 9")
                                tv = self.tv_max - 1
                            else:
                                print("\t\t\t\t 10")
                                tv = self.tv_max + 1

                else:
                    # newton interpolation
                    a_s = self.divided_diff(freq_list, tv_list)[0, :]

                    # use curve to compute for x
                    tv = self.newton_poly(a_s, freq_list, target_freq)

            elif self.f_max - self.f_min != 0:
                # switch control to just this part of the code for the end parts of the simulation for convergence
                print("\t\t\t\t 11")
                control_switch = False
                if self.tv_min*self.tv_max != 0:
                    print("\t\t\t\t 12")
                    tv = (self.tv_min + self.tv_max)/2

            # change tv
            par_end[indx], par_mid[indx] = tv, tv

            # if equal, break else continue with new shape
            error = target_freq - freq_list[-1]

            if n > max_iter:
                print('here at breaking')
                break

            # condition for repeated last four values
            if self.all_equal(freq_list[-2:]):
                print("Break because last two elements are equal")
                break

            if tv_list[-1] < 0:
                print("Negative value encountered. It is possible that there no solution for the parameter input set.")
                break

            # write output
            self.write_output(tv_list, freq_list, fid, projectDir)
            n += 1

        # return best answer from iteration
        min_error = [abs(x-target_freq) for x in freq_list]
        key = min_error.index(min(min_error))

        # import matplotlib.pyplot as plt
        # plt.scatter(tv_list, freq_list)
        # plt.show()

        return tv_list[key], freq_list[key]

    def tuneR(self, par_mid, par_end, target_freq, beampipes, bc, parentDir, projectDir, iter_set, proc=0):

        if proc == '':
            fid = '_process_0'
        else:
            fid = f'_process_{proc}'

        # get parameters
        freq_list = []
        Req_list = []
        error = 1
        print("in tuner ", projectDir)
        slans_geom.cavity(1, 1, par_mid, par_end, par_mid, f_shift=0, bc=bc, beampipes=beampipes, proc=proc, fid=fid, parentDir=parentDir, projectDir=projectDir)
        dirc = f'{projectDir}\SimulationData\SLANS\Cavity{fid}\cavity_{bc}.svl'
        d = fr.svl_reader(dirc)

        Req = par_end[6]
        freq = d['FREQUENCY'][0]
        freq_list.append(freq)
        Req_list.append(Req)

        # first shot
        Req = Req + 20
        par_end[6], par_mid[6] = Req, Req

        # iteration settings
        tol = iter_set[1]
        max_iter = iter_set[2]
        n = 0
        self.Req_min, self.Req_max = 0, 0
        self.f_min, self.f_max = 0, 0

        while abs(error) > tol:
            # run slans cavity code
            slans_geom.cavity(1, 1, par_mid, par_end, par_mid, f_shift=0, bc=bc, beampipes=beampipes, proc=proc, fid=fid, parentDir=parentDir, projectDir=projectDir)

            # get results and compare with set value
            d = fr.svl_reader(dirc)
            freq = d['FREQUENCY'][0]
            freq_list.append(freq)
            Req_list.append(Req)

            # update bounds
            if self.f_min < freq < target_freq and freq > self.f_min:
                self.f_min = freq
                self.Req_min = Req

            if target_freq < freq < self.f_max and freq < self.f_max:
                self.f_max = freq
                self.Req_max = Req

            # calculate slope of line from base point to new point
            if freq_list[n] - freq_list[n + 1] != 0:
                m = (Req_list[n + 1] - Req_list[n]) / (freq_list[n + 1] - freq_list[n])

                if iter_set[0] == "Linear Interpolation":
                    # calculate new Req with straight line formula
                    Req = Req_list[n+1] + m * (target_freq - freq_list[n+1])

                    # check if new Req is between the current min and max Req
                    if self.Req_min < abs(Req) < self.Req_max:
                        pass
                    else:
                        if self.Req_min*self.Req_max != 0:
                            Req = (self.Req_min + self.Req_max)/2

                        elif self.Req_min != 0:
                            # check the side of the curve the new value falls to
                            if abs(Req) > self.Req_min:
                                Req = max(self.Req_min + 1, self.Req_min*(target_freq/self.f_min))
                            else:
                                Req = max(self.Req_min - 1, self.Req_min/(target_freq/self.f_min))
                        elif self.Req_max != 0:
                            # check the side of the curve the new value falls to
                            if abs(Req) < self.Req_max:
                                Req = self.Req_max - 1
                            else:
                                Req = self.Req_max + 1

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
        '''
        function to calculate the divided
        differences table
        '''
        n = len(y)
        coef = np.zeros([n, n])
        # the first column is y
        coef[:, 0] = y

        for j in range(1, n):
            for i in range(n - j):
                coef[i][j] = (coef[i + 1][j - 1] - coef[i][j - 1]) / (x[i + j] - x[i]) if (x[i + j] - x[i]) != 0  else 0

        return coef

    @staticmethod
    def newton_poly(coef, x_data, x):
        '''
        evaluate the newton polynomial
        at x
        '''
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
        print("Writing output")
        dd = {"tv": tv_list, "freq": freq_list}

        with open(f"{projectDir}\SimulationData\SLANS\Cavity{fid}\convergence_output.json", "w") as outfile:
            json.dump(dd, outfile, indent=4, separators=(',', ': '))