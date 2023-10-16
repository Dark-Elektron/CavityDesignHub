from slans_geom_par import SLANSGeometry
from itertools import groupby
from utils.shared_functions import *

slans_geom = SLANSGeometry()


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
                          parentDir=parentDir, projectDir=projectDir, opt=True)
        dirc = fr'{projectDir}\{fid}\cavity_{bc}.svl'

        d = svl_reader(dirc)

        tv = par_end[indx]
        freq = d['FREQUENCY'][0]
        freq_list.append(freq)
        tv_list.append(tv)

        # first shot
        tv = tv + 2

        par_end[indx] = tv

        tol = iter_set[1]
        n = 0

        f_min, f_max = 0, 0

        while abs(error) > tol:
            # print('this point', par_mid, par_end)
            # run slans cavity code
            slans_geom.cavity(1, 1, par_mid, par_end, par_mid, f_shift=0, bc=bc,
                              beampipes=beampipes, proc=proc, fid=fid,
                              parentDir=parentDir, projectDir=projectDir, opt=True)

            # get results and compare with set value
            d = svl_reader(dirc)
            freq = d['FREQUENCY'][0]
            freq_list.append(freq)
            tv_list.append(tv)
            # update bounds
            if f_min < freq < target_freq and freq > f_min:
                f_min = freq
                # tv_min = tv

            if target_freq < freq < f_max and freq < f_max:
                f_max = freq

            # calculate slope of line from base point to new point

            if freq_list[n] - freq_list[n + 1] != 0:
                m = (tv_list[n + 1] - tv_list[n]) / (freq_list[n + 1] - freq_list[n])

                if iter_set[0] == "Linear Interpolation":
                    # calculate new Req with straight line formula
                    step = m * (target_freq - freq_list[n + 1])
                    if step > 10:
                        step = 10
                    if step < -10:
                        step = -10

                    tv = tv_list[n + 1] + step
                else:
                    # newton interpolation
                    a_s = self.divided_diff(freq_list, tv_list)[0, :]

                    # use curve to compute for x
                    tv = self.newton_poly(a_s, freq_list, target_freq)
            # change tv
            par_end[indx] = tv

            # if equal, break else continue with new shape
            error = target_freq - freq_list[-1]
            # print(error, tol)

            if n > 20:
                print('Maximum number of iterations exceeded. No solution found.')
                break

            # condition for repeated last four values
            if self.all_equal(freq_list[-2:]):
                print("Converged. Solution found.")
                break

            if tv_list[-1] < 0:
                print("Negative value encountered. It is possible that there no solution for the parameter input set.")
                break

            # check if Req is less than lower limit
            if tv < par_end[0]:
                print("Error: L < A. No solution.")
                break

            # check if alpha is less than or greater than 90.5
            # print('here')
            if par_mid == par_end:
                alpha, error_msg = calculate_alpha(par_mid[0], par_mid[1], par_mid[2],
                                                   par_mid[3], par_mid[4], tv, par_mid[6], 0)
                if alpha < 90.0 or error_msg != 1:
                    print("Mid cell alpha is less than 90", alpha, error_msg)
                    break
            else:
                alpha, error_msg = calculate_alpha(par_end[0], par_end[1], par_end[2],
                                                   par_end[3], par_end[4], tv, par_mid[6], 0)
                if alpha < 90.0 or error_msg != 1:
                    print("End cell alpha is less than 90", alpha, error_msg)
                    break

                alpha, error_msg = calculate_alpha(par_mid[0], par_mid[1], par_mid[2],
                                                   par_mid[3], par_mid[4], tv, par_mid[6], 0)
                if alpha < 90.0 or error_msg != 1:
                    print("Mid cell alpha is less than 90", alpha, error_msg)
                    break
            # print(alpha, error_msg)
            # update convergence list
            if conv_list:
                conv_list[proc] = [tv_list, freq_list]
            # self.write_output(tv_list, freq_list, fid, projectDir)

            n += 1

        # return best answer from iteration
        min_error = [abs(x - target_freq) for x in freq_list]
        key = min_error.index(min(min_error))

        # print(tv_list, freq_list)
        # import matplotlib.pyplot as plt
        # plt.scatter(tv_list, freq_list)
        # plt.show()

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
                          parentDir=parentDir, projectDir=projectDir, opt=True)
        dirc = fr'{projectDir}\SimulationData\SLANS_opt\{fid}\cavity_{bc}.svl'
        d = svl_reader(dirc)

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
                              fid=fid, parentDir=parentDir, projectDir=projectDir, opt=True)

            # get results and compare with set value
            d = svl_reader(dirc)
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
                    Req = Req_list[n + 1] + m * (target_freq - freq_list[n + 1])

                else:
                    # newton interpolation
                    a_s = self.divided_diff(freq_list, Req_list)[0, :]

                    # use curve to compute for x
                    Req = self.newton_poly(a_s, freq_list, target_freq)
            # change R
            par_end[6], par_mid[6] = Req, Req

            # if equal, break else continue with new shape
            error = target_freq - freq

            # avoid infinite loop-1
            if n > max_iter:
                break

            # check if alpha is less or greater than 90.5
            if par_mid == par_end:
                alpha, error_msg = calculate_alpha(par_mid[0], par_mid[1], par_mid[2],
                                                   par_mid[3], par_mid[4], par_mid[5], Req, 0)
                if alpha < 90.0 or error_msg != 1:
                    break
            else:
                alpha, error_msg = calculate_alpha(par_end[0], par_end[1], par_end[2],
                                                   par_end[3], par_end[4], par_end[5], Req, 0)
                if alpha < 90.0 or error_msg != 1:
                    break

                alpha, error_msg = calculate_alpha(par_mid[0], par_mid[1], par_mid[2],
                                                   par_mid[3], par_mid[4], par_mid[5], Req, 0)
                if alpha < 90.0 or error_msg != 1:
                    break

            # update convergence list
            if conv_list:
                conv_list[proc] = [Req_list, freq_list]

            n += 1

        # return best answer from iteration
        min_error = [abs(x - target_freq) for x in freq_list]
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

        with open(fr"{projectDir}\SimulationData\SLANS_opt\{fid}\convergence_output.json", "w") as outfile:
            json.dump(dd, outfile, indent=4, separators=(',', ': '))


def svl_reader(filename):
    dict = {
        'TITLE': [],
        'CAVITY RADIUS': [],
        'LENGTH': [],
        'FREQUENCY': [],
        'LENGTH OF WAVE': [],
        'WAVE VALUE': [],
        'QUALITY FACTOR': [],
        'STORED ENERGY': [],
        'TRANSIT TIME': [],
        'EFFECTIVE IMPEDANCE': [],
        'SHUNT IMPEDANCE': [],
        'MAXIMUM MAG. FIELD': [],
        'MAXIMUM ELEC. FIELD': [],
        'ACCELERATION': [],
        'ACCELERATION RATE': [],
        'AVERAGE E.FIELD ON AXIS': [],
        'KM (Emax/Accel.rate)': [],
        'KH (Hmax*Z0/Accel.rate)': [],
    }
    with open(filename, 'r') as f:
        data = f.readlines()
        for i, line in enumerate(data):
            if '*SLANS*' in line:
                dict['TITLE'].append(line)

            if 'CAVITY RADIUS' in line:
                dict['CAVITY RADIUS'].append(process_line(line, 'CAVITY RADIUS'))
                dict['LENGTH'].append(process_line(line, 'LENGTH'))

            if 'FREQUENCY' in line:
                dict['FREQUENCY'].append(process_line(line, 'FREQUENCY'))

            if 'LENGTH OF WAVE' in line:
                dict['LENGTH OF WAVE'].append(process_line(line, 'LENGTH OF WAVE'))

            if 'WAVE VALUE' in line:
                dict['WAVE VALUE'].append(process_line(line, 'WAVE VALUE'))

            if 'QUALITY FACTOR' in line:
                dict['QUALITY FACTOR'].append(process_line(line, 'QUALITY FACTOR'))

            if 'STORED ENERGY' in line:
                dict['STORED ENERGY'].append(process_line(line, 'STORED ENERGY'))

            if 'TRANSIT TIME' in line:
                dict['TRANSIT TIME'].append(process_line(line, 'TRANSIT TIME'))

            if 'EFFECTIVE IMPEDANCE' in line:
                dict['EFFECTIVE IMPEDANCE'].append(process_line(line, 'EFFECTIVE IMPEDANCE'))

            if 'SHUNT IMPEDANCE' in line:
                dict['SHUNT IMPEDANCE'].append(process_line(line, 'SHUNT IMPEDANCE'))

            if 'MAXIMUM MAG. FIELD' in line:
                dict['MAXIMUM MAG. FIELD'].append(process_line(line, 'MAXIMUM MAG. FIELD'))

            if 'MAXIMUM ELEC.FIELD' in line:
                dict['MAXIMUM ELEC. FIELD'].append(process_line(line, 'MAXIMUM ELEC.FIELD'))

            if 'ACCELERATION' in line and not 'RATE' in line:
                dict['ACCELERATION'].append(process_line(line, 'ACCELERATION'))

            if 'ACCELERATION RATE' in line:
                dict['ACCELERATION RATE'].append(process_line(line, 'ACCELERATION RATE'))

            if 'AVERAGE E.FIELD ON AXIS' in line:
                dict['AVERAGE E.FIELD ON AXIS'].append(process_line(line, 'AVERAGE E.FIELD ON AXIS'))

            if 'KM (Emax/Accel.rate)' in line:
                dict['KM (Emax/Accel.rate)'].append(process_line(line, 'KM (Emax/Accel.rate)'))

            if 'KH (Hmax*Z0/Accel.rate)' in line:
                dict['KH (Hmax*Z0/Accel.rate)'].append(process_line(line, 'KH (Hmax*Z0/Accel.rate)'))

    return dict


def process_line(line, request):
    # select substring from index
    line = line[line.index(request):]
    line = line.strip().split(" ")
    # print(line)
    res = 0
    for val in line:
        try:
            res = float(val)
            break
        except:
            continue
    return res
