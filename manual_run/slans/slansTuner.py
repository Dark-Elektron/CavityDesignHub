import os
import subprocess


class SLANSTune:
    def __init__(self):
        self.filename = 'tuned'

    def mid_cell_tune(self, A, B, a, b, Ri, L, Req_0, freq0, f_shift=0, n_modes=1, beta=1):
        cwd = os.getcwd()
        self.tuner = fr'{cwd}\SLANS_exe\MidCellTune\Superlans_Files\TunedCell.exe'
        self.slans_files = fr'{cwd}\SLANS_exe\MidCellTune\Superlans_Files'

        self.freq0 = freq0
        print("\t\tSLANS_TUNE:: Started tuning, initial values -> ", A, B, a, b, Req_0, Ri, L, freq0)
        self.write_geometry_parameters_mid_tune(A, B, a, b, Req_0, Ri, L)

        print("\t\tDone writing geometry parameters")
        self.write_beta_mid_Tune(beta, f_shift, n_modes)
        print("\t\tDone writing tuned parameters")

        filepath = fr'{self.slans_files}\tuned'
        print("\t\tCalling subprocess, TunedCell\n")
        subprocess.run([self.tuner, filepath, '-b'], cwd=self.slans_files)
        print("\t\t\tDone tuning")

        R_eq, freq, alpha, h, e = self.read_tuned_data()
        print("\t\t", R_eq, freq, alpha)

        return R_eq, freq, alpha, h, e

    def end_cell_tune(self, par_in, par_out, freq0, f_shift=0, n_modes=1, beta=1):
        cwd = os.getcwd()
        self.tuner_end = fr'{cwd}\SLANS_exe\EndCellTune\Superlans_Files\TunedCellEnd_120421.exe'
        self.slans_files_end = fr'{cwd}\SLANS_exe\EndCellTune\Superlans_Files'

        self.freq0 = freq0
        print("\t\tSLANS_TUNE:: Started tuning, initial values -> ", par_in, par_out, freq0)
        self.write_geometry_parameters_end_tune(par_in, par_out)

        print("\t\tDone writing geometry parameters")
        self.write_beta_end_Tune(beta, f_shift, n_modes)
        print("\t\tDone writing tuned parameters")

        filepath = fr'{self.slans_files_end}\tesla_end2'

        print("\t\tCalling subprocess, TunedCell\n")
        print("\t\t", filepath)
        subprocess.run([self.tuner_end, filepath, '-b'], cwd=self.slans_files_end)
        print("\t\t\tDone tuning")

        L, freq, alpha, h, e = self.read_tuned_data_end()
        print("\t\t", L, freq, alpha)

        return L, freq, alpha, h, e

    def write_geometry_parameters_mid_tune(self, A, B, a, b, Req0, Rbp, L):
        file_text = []
        filename = fr'{self.slans_files}\tuned_reference.tgp'
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
        A, B, a, b, Ri, L, Req_i = par_in
        Ae, Be, ae, be, Rie, Le, Req_e = par_out
        at, bt, c = 0, 0, 0
        with open(fr'{self.slans_files_end}\tesla_end2.tgpe', 'w') as f:
            f.write('{:g} {:g}  {:g} : Req(mm), Rae(mm), Le(mm) - parameters of the end half-cell \n'.format(Req_i, Rie, Le))
            f.write('{:g}  {:g}  {:g}  {:g}  {:g}  {:g}  {:g} : Al(mm), Bl(mm), al(mm), bl(mm), ll(mm), Ral(mm), Ll(mm) - parameters of the inner half-cell\n'.format(A, B, a, b, 32.408, Ri, L))
            f.write('{:g}    5            : Freq(MHz), mode\n'.format(self.freq0))
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
        filename = fr'{self.slans_files_end}\tesla_end2.dtr'
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
                    data[key]= float(var_values[v])
                    v += 1
                except:
                    continue
        # This could be modified to return more parameter values
        return data['Req'], data['Freq'], data['alpha'], data['h'], data['e']

    def read_tuned_data_end(self):
        filename = fr'{self.slans_files_end}\tesla_end2.tgpe'
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


if __name__ == '__main__':
    # Create tuner object
    tune = SLANSTune()

    # Uncomment to run mid cell tune
    # Call mid_cell_tune function
    Req, freq, alpha, h, e = tune.mid_cell_tune(67.72, 67.72, 21.75, 21.75, 60, 93.5, 160, 801.58)
    print(Req, freq, alpha, h, e)

    # Note: The tune for end cell still gives zeros for many inpus. I think that this might be because of the variable
    # "ll" on line 100 (...Al(mm), Bl(mm), al(mm), bl(mm), ll(mm), Ral(mm), Ll(mm) .... I don't know what small "ll" stands for.
    # It is set to a constant value.
    # # Uncomment to run end cell tune
    # # Call end_cell_tune function
    # L, freq, alpha, h, e = tune.end_cell_tune([67.72, 67.72, 21.75, 21.75, 60, 93.5, 171.381], [52.0, 52.0, 28.5, 28.5, 82, 93.5, 169.476], 801.58) # par_out = [Ae, Be, ae, be, Rie, Le]
    # print(L, freq, alpha, h, e)
