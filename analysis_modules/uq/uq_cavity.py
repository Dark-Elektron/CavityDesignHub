import pandas as pd
from analysis_modules.eigenmode.SLANS.slans_geometry import SLANSGeometry
from utils.file_reader import FileReader
from utils.shared_classes import *
from utils.shared_functions import *
import numpy as np
from icecream import ic

slans_geom = SLANSGeometry()
fr = FileReader()

file_color = 'green'
DEBUG = True


def print_(*arg):
    if DEBUG:
        print(colored(f'\t{arg}', file_color))


class UQ:
    def __init__(self):
        pass

    def uq(self, key, shape, qois, f_shift=801.58e6, bc='mm', parentDir='', projectDir='', rand_vars=None, constraint=None):
        """

        Parameters
        ----------
        key: str

        shape: list, array like

        qois: list, array like
            Quantities of interest

        f_shift: float
            Frequency shift
        bc: str
            Bounday condition, 1:inner contour, 2:Electric wall Et = 0, 3:Magnetic Wall En = 0, 4:Axis, 5:metal
        parentDir: str
            Parent directory
        projectDir: str
            Project directory
        rand_vars: list, array like
            Random variables
        constraint: list, array like
            Constraints

        Returns
        -------

        """
        if rand_vars is None:
            rand_vars = ['A', 'B', 'a', 'b']
            rand_vars = ['A', 'Ri', "Req"]
            rand_vars = ['A']

        err = False
        result_dict_slans = {}
        slans_obj_list = qois
        for o in qois:
            result_dict_slans[o] = {'expe': [], 'stdDev': []}

        # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
        p_true = shape['GEOM']
        df = pd.DataFrame(p_true, columns=["A", "B", "a", 'b', 'Ri', 'L', 'Req', 'alpha'])
        # print(df)
        # get the random variables from the list
        df_rv = df[rand_vars]

        # print(df_rv)
        # print(df_rv.to_numpy())
        # print(df_rv.to_numpy().flatten())
        df_rv_shape = df_rv.shape
        p_true = df_rv.to_numpy().flatten()
        init_len = len(p_true)

        # apply constraint(s)
        pop_indx, non_pop_indx = [], []
        all_indxs = []
        if constraint:
            for k, val in constraint.items():
                # print(k)
                # get index from rand_vars
                ind = rand_vars.index(k)
                # print(ind)
                # print(np.array(val).flatten()*len(rand_vars) + ind)
                #
                all_indx = np.array(val).flatten()*len(rand_vars) + ind
                pop_indx.extend([p for p in all_indx if p % 2 != 0])
                non_pop_indx.extend([p for p in all_indx if p % 2 == 0])
                all_indxs.extend(all_indx)
                # print(pop_indx)
                # p_true.pop()
            print(pop_indx)
            print(non_pop_indx)

        pop_indx.sort(), non_pop_indx.sort()
        print(f'\t\t', pop_indx)
        # remove constraint index from p_true
        ic(p_true, len(p_true))
        p_true = np.delete(p_true, pop_indx)
        ic(p_true, len(p_true))
        rdim = len(p_true)  # How many variables will be considered as random in our case 5
        degree = 1

        #  for 1D opti you can use stroud5 (please test your code for stroud3 less quadrature nodes 2rdim)
        flag_stroud = 1
        if flag_stroud == 1:
            nodes, weights, bpoly = quad_stroud3(rdim, degree)
            nodes = 2. * nodes - 1.
        elif flag_stroud == 2:
            nodes, weights, bpoly = quad_stroud3(rdim, degree)  # change to stroud 5 later
            nodes = 2. * nodes - 1.
        else:
            ic('flag_stroud==1 or flag_stroud==2')

        #  mean value of geometrical parameters
        # p_init = np.zeros(np.shape(p_true))

        no_parm, no_sims = np.shape(nodes)
        # ic(nodes)
        delta = 0.03  # or 0.1

        Ttab_val_f = []
        # print_('3')
        sub_dir = fr'{key}'  # the simulation runs at the quadrature points are saved to the key of mean value run
        for i in range(no_sims):
            skip = False
            p_init = p_true * (1 + delta * nodes[:, i])
            # ic(p_init, len(p_init))

            # reintroduce constraints
            pp = []
            inc = 0
            # ic(pop_indx, non_pop_indx)
            for nn in range(init_len):
                if nn in pop_indx:
                    pp.insert(pop_indx[inc], pp[non_pop_indx[inc]])
                    # print("\t\t\t", pop_indx[inc], pp[non_pop_indx[inc]], inc)
                    # print("\t\t\t", pp)
                    inc += 1
                else:
                    pp.append(p_init[nn - inc])

                # print(f'{len(pp)}', nn - inc)

            # ic(pp, len(pp))
            p_init = np.array(pp)

            # reshape p_init
            p_init = p_init.reshape(df_rv_shape)
            # print(p_init)

            # update df with new values
            df[rand_vars] = p_init
            ic(df)

            cells_par = df.to_numpy()

            # perform checks on geometry
            for inds, row in df.iterrows():
                ok = self.perform_geometry_checks(row)

            # print_("OK", ok)
            if not ok:
                err = True
                break
            fid = fr'{key}_Q{i}'

            # check if folder exists and skip if it does
            # print_(fr'{projectDir}\SimulationData\SLANS\{key}\{fid}')
            if os.path.exists(fr'{projectDir}\SimulationData\SLANS\{key}\{fid}'):
                skip = True
                # ic("Skipped: ", fid, fr'{projectDir}\SimulationData\ABCI\{key}\{fid}')

            # skip analysis if folder already exists.
            if not skip:
                #  run model using SLANS or CST
                # # create folders for all keys
                slans_geom.createFolder(fid, parentDir, subdir=sub_dir)

                slans_geom.cavity_multicell_full(no_of_modules=1, cells_par=cells_par, fid=fid, bc=33,
                                                 f_shift='default', beta=1, n_modes=6, beampipes="both",
                                                 parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)

            filename = fr'{parentDir}\SimulationData\SLANS\{key}\{fid}\qois.json'
            # ic(filename)
            if os.path.exists(filename):
                # n_cells = 5
                # params = fr.svl_reader(filename)
                # norm_length = 2 * n_cells * shape['IC'][5]
                # ic(n_cells, norm_length)
                # qois_result = self.get_qoi_value(params, slans_obj_list, n_cells, norm_length)

                with open(fr'{parentDir}\SimulationData\SLANS\{key}\{fid}\qois.json') as json_file:
                    results = json.load(json_file)

                # print(results)
                qois_result = [results[x] for x in qois]
                # print_(qois_result)
                # sometimes some degenerate shapes are still generated and the solver returns zero
                # for the objective functions, such shapes are considered invalid
                for objr in qois_result:
                    if objr == 0:
                        # skip key
                        err = True
                        break

                tab_val_f = qois_result

                Ttab_val_f.append(tab_val_f)
                # ic(Ttab_val_f)
            else:
                err = True

        # # add original point
        # filename = fr'{projectDir}\SimulationData\SLANS\{key}\cavity_33.svl'
        # params = fr.svl_reader(filename)
        # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
        # tab_val_f = obj_result
        # Ttab_val_f.append(tab_val_f)

        print_("Error: ", err)
        # import matplotlib.pyplot as plt
        ic(Ttab_val_f)
        if not err:
            v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights)
            # ic(v_expe_fobj, v_stdDev_fobj)
            # append results to dict
            ic(v_expe_fobj, v_stdDev_fobj)
            for i, o in enumerate(slans_obj_list):
                result_dict_slans[o]['expe'].append(v_expe_fobj[i])
                result_dict_slans[o]['stdDev'].append(v_stdDev_fobj[i])

                # pdf = normal_dist(np.sort(np.array(Ttab_val_f).T[i]), v_expe_fobj[i], v_stdDev_fobj[i])
                # plt.plot(np.sort(np.array(Ttab_val_f).T[i]), pdf)

            # plt.show()

            with open(fr"{parentDir}\SimulationData\SLANS\{key}\uq.json", 'w') as file:
                file.write(json.dumps(result_dict_slans, indent=4, separators=(',', ': ')))
        else:
            print_(fr"There was a problem running UQ analysis for {key}")

    @staticmethod
    def get_qoi_value(d, obj, n_cells, norm_length):
        """Get SLANS quantities of itnerest

        Parameters
        ----------
        d: dict

        obj:

        n_cells: int

        norm_length: float
            Normalisation length :math: `E_\\mathrm{acc} = \frac{V_{acc}}{L_{norm}`.

        Returns
        -------

        """
        Req = d['CAVITY RADIUS'][n_cells - 1] * 10  # convert to mm
        Freq = d['FREQUENCY'][n_cells - 1]
        E_stored = d['STORED ENERGY'][n_cells - 1]
        Rsh = d['SHUNT IMPEDANCE'][n_cells - 1]  # MOhm
        Q = d['QUALITY FACTOR'][n_cells - 1]
        Epk = d['MAXIMUM ELEC. FIELD'][n_cells - 1]  # MV/m
        Hpk = d['MAXIMUM MAG. FIELD'][n_cells - 1]  # A/m
        # Vacc = dict['ACCELERATION'][0]
        Eavg = d['AVERAGE E.FIELD ON AXIS'][n_cells - 1]  # MV/m
        Rsh_Q = d['EFFECTIVE IMPEDANCE'][n_cells - 1]  # Ohm

        Vacc = np.sqrt(
            2 * Rsh_Q * E_stored * 2 * np.pi * Freq * 1e6) * 1e-6
        # factor of 2, remember circuit and accelerator definition
        # Eacc = Vacc / (374 * 1e-3)  # factor of 2, remember circuit and accelerator definition
        Eacc = Vacc / (norm_length * 1e-3)  # for 1 cell factor of 2, remember circuit and accelerator definition
        Epk_Eacc = Epk / Eacc
        Bpk_Eacc = (Hpk * 4 * np.pi * 1e-7) * 1e3 / Eacc

        d = {
            "Req": Req,
            "freq": Freq,
            "Q": Q,
            "E": E_stored,
            "R/Q": 2 * Rsh_Q,
            "Epk/Eacc": Epk_Eacc,
            "Bpk/Eacc": Bpk_Eacc
        }

        objective = []

        # append objective functions
        for o in obj:
            if o in d.keys():
                objective.append(d[o])

        return objective

    @staticmethod
    def perform_geometry_checks(par_half_cell):
        """

        Parameters
        ----------
        par_half_cell: list, array like
            Geometric parameters of half call of cavity: [A, B, a, b, Ri, L, Req, alpha]

        Returns
        -------

        """
        # # check if Req is less than lower limit
        # if par_mid[6] < par_mid[1] + par_mid[3] + par_mid[4] or par_end[6] < par_end[1] + par_end[3] + par_end[4]:
        #     return False

        # check if alpha is less than 90.5
        alpha, error_msg = calculate_alpha(par_half_cell[0], par_half_cell[1], par_half_cell[2], par_half_cell[3],
                                           par_half_cell[4], par_half_cell[5],
                                           par_half_cell[6], 0)
        if alpha < 90.0 or error_msg != 1:
            print("1:", alpha, error_msg)
            return False

        # check if L is less than lower limit
        if par_half_cell[5] < par_half_cell[0]:
            print("4:", alpha, error_msg)
            return False

        return True


if __name__ == '__main__':
    uq = UQ()

    key = "Multicell_UQ_test"
    shape = {
        "GEOM": [
            [64, 58, 17, 12, 81, 96.45, 171.36, 0],
            [62, 66, 30, 23, 72, 93.5, 171.36, 0],
            [62, 66, 30, 23, 72, 93.5, 171.36, 0],
            [62, 66, 30, 23, 72, 93.5, 171.36, 0],
            [62, 66, 30, 23, 72, 93.5, 171.36, 0],
            [62, 66, 30, 23, 72, 93.5, 171.36, 0],
            [62, 66, 30, 23, 72, 93.5, 171.36, 0],
            [62, 66, 30, 23, 72, 93.5, 171.36, 0],
            [62, 66, 30, 23, 72, 93.5, 171.36, 0],
            [64, 58, 17, 12, 81, 96.45, 171.36, 0]
        ],
        "BP": "both",
        "FREQ": 801.5796
    }

    with open(fr'D:\Dropbox\CavityDesignHub\Cavity800\Cavities\multicell.json') as json_file:
        shapes = json.load(json_file)

    parentDir = r"D:\Dropbox\Files\Test_multicell"
    projectDir = r"SimulationData\SLANS"
    qois = ['freq [MHz]', 'R/Q [Ohm]', "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]"]

    no_of_cells = len(shape['GEOM'])/2
    Req_pairs = [[2 * i, 2 * i + 1] for i in range(int(no_of_cells))]
    Ri_pairs = [[2 * i - 1, 2 * i] for i in range(1, int(no_of_cells))]
    print(Ri_pairs)
    print(Req_pairs)
    constraint = {
        "Ri": Ri_pairs,
        "Req": Req_pairs
    }
    constraint = None
    for key, shape in shapes.items():
        uq.uq(key, shape, qois, parentDir=r"D:\Dropbox\Files\Test_multicell", projectDir=r"SimulationData\SLANS",
              constraint=constraint)
