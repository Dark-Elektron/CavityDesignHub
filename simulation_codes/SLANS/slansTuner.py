import subprocess
import numpy as np
from scipy.optimize import fsolve


class SLANSTune:
    def __init__(self, parentDir, projectDir):
        self.parentDir = parentDir
        self.projectDir = projectDir
        self.filename = 'tuned'

    def mid_cell_tune(self, A, B, a, b, Ri, L, Req_0, freq0, proc=0, f_shift=0, n_modes=1, beta=1):
        self.tuner = fr'{self.projectDir}\SimulationData\SLANS\Cavity_process_{proc}\SLANS_exe\MidCellTune\Superlans_Files\TunedCell.exe'
        self.slans_files = fr'{self.projectDir}\SimulationData\SLANS\Cavity_process_{proc}\SLANS_exe\MidCellTune\Superlans_Files'

        self.freq0 = freq0
        # print("\t\tSLANS_TUNE:: Started tuning, initial values -> ", A, B, a, b, Req_0, Ri, L, freq0)
        self.write_geometry_parameters_mid_tune(A, B, a, b, Req_0, Ri, L)

        # print("\t\tDone writing geometry parameters")
        self.write_beta_mid_Tune(beta, f_shift, n_modes)
        # print("\t\tDone writing tuned parameters")

        filepath = fr'{self.slans_files}\{self.filename}'
        # print("\t\tCalling subprocess, TunedCell\n")

        # the next two lines suppress pop up windows from the slans codes
        # the slans codes, however, still disrupts windows operation, sadly. This is the case even for the slans tuner
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        cwd = fr'{self.projectDir}\SimulationData\SLANS\Cavity_process_{proc}\SLANS_exe\MidCellTune'
        subprocess.call([self.tuner, filepath, '-b'], cwd=cwd, startupinfo=startupinfo)
        # print("\t\t\tDone tuning")

        R_eq, freq, alpha, h, e = self.read_tuned_data()
        # print("\t\t", R_eq, freq, alpha)

        return R_eq, freq, alpha, h, e

    def end_cell_tune(self, par_in, par_out, freq0, proc, f_shift=0, n_modes=1, beta=1):
        print("here in end cell tuner", proc)
        self.tuner_end = fr'{self.projectDir}\SimulationData\SLANS\Cavity_process_{proc}\SLANS_exe\EndCellTune\Superlans_Files\TunedCellEnd_120421.exe'
        self.slans_files_end = fr'{self.projectDir}\SimulationData\SLANS\Cavity_process_{proc}\SLANS_exe\EndCellTune\Superlans_Files'

        self.freq0 = freq0
        print("\t\tSLANS_TUNE:: Started tuning, initial values -> ", par_in, par_out, freq0)
        self.write_geometry_parameters_end_tune(par_in, par_out)

        print("\t\tDone writing geometry parameters")
        self.write_beta_end_Tune(beta, f_shift, n_modes)
        print("\t\tDone writing tuned parameters")

        filepath = fr'{self.slans_files_end}/{self.filename}'

        print("\t\tCalling subprocess, TunedCell\n")
        cwd = fr'{self.projectDir}\SimulationData\SLANS\Cavity_process_{proc}\SLANS_exe\EndCellTune\Superlans_Files'
        print(filepath)
        subprocess.call([self.tuner_end, filepath], cwd=cwd)
        print("\t\t\tDone tuning")

        L, freq, alpha, h, e = self.read_tuned_data_end()
        print("\t\t", L, freq, alpha)

        return L, freq, alpha, h, e

    def write_geometry_parameters_mid_tune(self, A, B, a, b, Req0, Rbp, L):
        file_text = []
        filename = fr'{self.slans_files}\{self.filename}.tgp'
        with open(filename, 'r') as f:
            tline = f.readline().rstrip('\n')
            while isinstance(tline, str):
                k = tline.find('------------------------')
                t = tline.find('delta(mm)')
                file_text.append(tline)
                if t != -1:
                    tline = f.readline()
                    file_text.append(f"{A: 0.8f}")
                    tline = f.readline()
                    file_text.append(f"{B: 0.8f}")
                    tline = f.readline()
                    file_text.append(f"{a: 0.8f}")
                    tline = f.readline()
                    file_text.append(f"{b: 0.8f}")

                tline = f.readline().rstrip('\n')
                if k != -1:
                    break

            file_text.append(tline)
            tline = f.readline()
            file_text.append(tline)

        file_text[0] = '{:g}    {:g}   {:g} : Req(mm), Rbp(mm), L(mm)'.format(Req0, Rbp, L)
        file_text[1] = '{:0.2f}     1 : Freq(MHz), mode'.format(self.freq0)

        filename = fr'{self.slans_files}\tuned.tgp'
        with open(filename, 'w') as f:
            f.write('{}\n'.format(file_text[0]))
        with open(filename, 'a') as f:
            for i in range(1, len(file_text)):
                f.write('{}\n'.format(file_text[i]))

    def write_geometry_parameters_end_tune(self, par_in, par_out):
        A, B, a, b, Ri, L, Req_i, _ = par_in
        li = self.calculate_li(A, B, a, b, Ri, L, Req_i, 0)

        Ae, Be, ae, be, Rie, Le, Req_e, _ = par_out
        at, bt, c = 0, 0, 0
        with open(fr'{self.slans_files_end}\{self.filename}.tgpe', 'w') as f:
            f.write('{:g} {:g}  {:g} : Req(mm), Rae(mm), Le(mm) - parameters of the end half-cell \n'.format(Req_i, Rie, Le))
            f.write('{:g}  {:g}  {:g}  {:g}  {:g}  {:g}  {:g} : Al(mm), Bl(mm), al(mm), bl(mm), ll(mm), Ral(mm), Ll(mm) - parameters of the inner half-cell\n'.format(A, B, a, b, li, Ri, L))
            f.write('{:g}    1            : Freq(MHz), mode\n'.format(self.freq0))
            f.write('1.           16                  : delta(mm), m\n')
            f.write('{:g}     90      0      : ZA1(mm), A2(mm), nA\n'.format(Ae))
            f.write('{:g}      80     0      : B1(mm), B2(mm), nB\n'.format(Be))
            f.write('{:g}        20          0 : a1(mm), a2(mm), na\n'.format(ae))
            f.write('{:g}    13.4      0: b1(mm), b2(mm), nb\n'.format(be))
            f.write('2.00         2.20     -1 : l1(mm), l2(mm), nl\n')
            f.write('1.e-7  50  12: delta freq, f_iter, n_l\n')
            f.write('-1.0        : r - reentrant 1.0 or not -1.0 or 0.0 for automatic detection\n')
            f.write('54.9914  50.7779   1 : y_begin(mm), y_end(mm), ny; automatic if ny = 1\n')
            f.write('{:g}   {:g}     {:g}     4      : at(mm), bt(mm), c(mm), Lbp / Rbp\n'.format(at, bt, c))

            f.write('----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
            f.write('      A           B          a          b           l              L           Freq          h            rh        zh           e               re        ze        Q         Rsh         ee             hh             it  r   alpha\n')
            f.write('-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')

    def write_beta_mid_Tune(self, beta, f_shift, n_modes):
        filename = fr'{self.slans_files}\{self.filename}.dtr'
        with open(filename, 'w') as f:
            f.write(':          Date:02/04/16 \n')
            f.write('{:g} :number of iterative modes 1-10\n'.format(n_modes))
            f.write('{:g} :number of search modes\n'.format(n_modes - 1))
            f.write('9.99999997E-007 :convergence accuracy\n')
            f.write('50 :maximum number of iterations\n')
            f.write('0 :continue iterations or not 1,0\n')
            f.write(' {:g}. :initial frequency shift MHz\n'.format(f_shift))
            f.write('1 :wave type 1-E, 2-H\n')
            f.write(' 1 :struct. 1-cav,2-per.str,3-w.guid.,4-l.-hom.\n')
            f.write('1 :symmetry yes or not 1,0\n')
            f.write(' 1 :number of met.surfaces, then:sign and sigma\n')
            f.write('5  1.\n')
            f.write('0 : number of mark volumes,then:sign,EPS,MU,TGE,TGM\n')
            f.write('{:g} : beta (v/c)\n'.format(beta))

    def write_beta_end_Tune(self, beta, f_shift, n_modes):
        filename = fr'{self.slans_files_end}\{self.filename}.dtr'
        with open(filename, 'w') as f:
            f.write(':          Date:02/04/16 \n')
            f.write('{:g} :number of iterative modes 1-10\n'.format(n_modes+1))
            f.write('{:g} :number of search modes\n'.format(n_modes))
            f.write('9.99999997E-007 :convergence accuracy\n')
            f.write('50 :maximum number of iterations\n')
            f.write('0 :continue iterations or not 1,0\n')
            f.write(' {:g}. :initial frequency shift MHz\n'.format(f_shift))
            f.write('1 :wave type 1-E, 2-H\n')
            f.write(' 1 :struct. 1-cav,2-per.str,3-w.guid.,4-l.-hom.\n')
            f.write('0 :symmetry yes or not 1,0\n')
            f.write(' 1 :number of met.surfaces, then:sign and sigma\n')
            f.write('5  1.\n')
            f.write('0 : number of mark volumes,then:sign,EPS,MU,TGE,TGM\n')
            f.write('{:g} : beta (v/c)\n'.format(beta))

    def read_tuned_data(self):
        filename = fr'{self.slans_files}\{self.filename}.tgp'
        with open(filename, 'r') as f:

            tline = f.readline()
            while isinstance(tline, str):
                k = tline.find('------------------------')
                tline = f.readline()

                if k != -1:
                    break
            var_names = tline.strip().split(' ')
            data = {}
            for name in var_names:
                if name != '':
                    data[name] = 0

            tline = f.readline()
            tline = f.readline()
            tline = f.readline()

            var_values = tline.strip().split(' ')
            v = 0
            for key, val in data.items():
                try:
                    data[key] = float(var_values[v])
                    v += 1
                except:
                    continue
        # This could be modified to return more parameter values
        return data['Req'], data['Freq'], data['alpha'], data['h'], data['e']

    def read_tuned_data_end(self):
        filename = fr'{self.slans_files_end}\{self.filename}.tgpe'
        print("Reading data from:: ", filename)
        with open(filename, 'r') as f:

            tline = f.readline()
            while isinstance(tline, str):
                k = tline.find('------------------------')
                tline = f.readline()

                if k != -1:
                    break
            var_names = tline.strip().split(' ')
            data = {}
            for name in var_names:
                if name != '':
                    data[name] = 0

            tline = f.readline()
            tline = f.readline()

            var_values = tline.strip().split(' ')
            v = 0
            for key, val in data.items():
                try:
                    data[key] = float(var_values[v])
                    v += 1
                except:
                    continue

        print(data)
        print(data["L"], data['h'], data['e'], data['alpha'])

        # This could be modified to return more parameter values
        return data['L'], data['Freq'], data['alpha'], data['h'], data['e']

    def calculate_li(self, A, B, a, b, Ri, L, Req, L_bp):
        data = ([0 + L_bp, Ri + b, L + L_bp, Req - B],
                [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        x1, y1, x2, y2 = fsolve(self.ellipse_tangent,
                                np.array([a + L_bp, Ri + 0.85 * b, L - A + L_bp, Req - 0.85 * B]),
                                args=data)

        li = np.sqrt((x1-x2)**2 + (y1-y2)**2)
        return li

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


if __name__ == '__main__':
    parentDir = "D:\Dropbox\CavityDesignHub" # '<parentDir>\SimulationData\SLANS\Cavity_process_{proc}\SLANS_exe\MidCellTune\Superlans_Files\TunedCell.exe'
    projectDir = "D:\Dropbox\CavityDesignHub\SampleProject" # '<projectDir>\SimulationData\SLANS\Cavity_process_{proc}\SLANS_exe\MidCellTune\Superlans_Files'
    "D:\Dropbox\CavityDesignHub\SampleProject\SimulationData\SLANS\Cavity_process_0\SLANS_exe\EndCellTune\Superlans_Files\TunedCellEnd_120421.exe"
    # Create tuner object
    tune = SLANSTune(parentDir, projectDir)

    # # Call mid_cell_tune function
    # L, freq, alpha, h, e = tune.mid_cell_tune(67.72, 67.72, 21.75, 21.75, 60, 93.5, 160, 801.58)
    # # L, freq, alpha, h, e = tune.end_cell_tune([52, 52, 28.5, 28.5, 55, 93.5, 169.476], [52.0, 52.0, 28.5, 28.5, 82, 93.5, 169.476], 801.58) # par_out = [Ae, Be, ae, be, Rie, Le, at, bt, c]
    #
    # print(L, freq, alpha, h, e)


    mid_cell_par = [71.25, 49.87, 20, 28.07, 65, 93.5, 165.355, 0] # par_in = [A, B, a, bi, R, L, Req]
    end_cell_par = [72.500, 52.031, 19.966, 18.202, 75, 93.5, 165.355, 0] # par_out = [Ae, Be, ae, be, Rie, Le, Req]
    target_freq = 801.58 # MHz
    L, freq, alpha, h, e = tune.end_cell_tune(mid_cell_par, end_cell_par, target_freq, 0)
    print(L, freq, alpha, h, e)

# D:\Dropbox\CavityDesignHub\SampleProject\SimulationData\SLANS\Cavity_process_0\SLANS_exe\EndCellTune\Superlans_Files\tuned
# D:\Dropbox\CavityDesignHub\SampleProject\SimulationData\SLANS\Cavity_process_0\SLANS_exe\EndCellTune\Superlans_Files\tuned