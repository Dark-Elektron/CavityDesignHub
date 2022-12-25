import time
from analysis_modules.tune.tuners.pyTuner import PyTune
from analysis_modules.eigenmode.SLANS.slansTuner import SLANSTune
from utils.shared_classes import *
from utils.shared_functions import *
from utils.file_reader import FileReader
from analysis_modules.eigenmode.SLANS.slans_geom_par import SLANSGeometry

fr = FileReader()
slans_geom = SLANSGeometry()


class Tuner:
    def __init__(self):
        pass

    def tune(self, pseudo_shape_space, bc, parentDir, projectDir, filename, resume="No",
             proc=0, tuner_option='SLANS', tune_variable='Req', iter_set=None, cell_type='Mid Cell',
             progress_list=None, convergence_list=None, save_last=True, n_cell_last_run=1):

        # tuner
        pytune = PyTune()

        if tuner_option == 'SLANS':
            slans_tune = SLANSTune(parentDir, projectDir)
        else:
            slans_tune = None

        start = time.time()
        population = {}
        total_no_of_shapes = len(list(pseudo_shape_space.keys()))

        # check for already processed shapes
        existing_keys = []

        if resume == "Yes":
            # check if value set is already written. This is to enable continuation in case of break in program
            if os.path.exists(fr'{projectDir}/Cavities/{filename}'):
                population = json.load(open(fr'{projectDir}/Cavities/{filename}', 'r'))

                existing_keys = list(population.keys())
                # print(f'Existing keys: {existing_keys}')

        progress = 0
        error_msg1 = 1
        error_msg2 = 1

        for key, pseudo_shape in pseudo_shape_space.items():
            # print(pseudo_shape['IC'])
            A_i, B_i, a_i, b_i, Ri_i, L_i, Req, _, _ = pseudo_shape['IC']  # Req here is none but required
            # since the shape space consists of all variables
            A_o, B_o, a_o, b_o, Ri_o, L_o, Req_o, _, _ = pseudo_shape['OC']  # Req here is none but required
            # since the shape space consists of all variables
            beampipes = pseudo_shape['BP']
            target_freq = pseudo_shape['FREQ']

            # Check if simulation is already run
            freq = 0
            alpha_i = 0
            alpha_o = 0
            if resume == "Yes" and os.path.exists(fr'{projectDir}/SimulationData/SLANS/Cavity{key}'):
                # if folder exist, read value
                filename = fr'{projectDir}/SimulationData/SLANS/Cavity{key}/cavity_{bc}.svl'
                try:
                    data_dict = fr.svl_reader(filename)
                    # print(data_dict)
                    if tune_variable == 'Req':
                        Req = data_dict['CAVITY RADIUS'][0]*10
                        freq = data_dict['FREQUENCY'][0]
                    else:
                        L = data_dict['LENGTH'][0]*10
                        freq = data_dict['FREQUENCY'][0]

                        if cell_type == 'Mid Cell':
                            L_i, L_o = L, L
                        else:
                            L_o = L
                except FileNotFoundError:
                    if tune_variable == 'Req':
                        Req = 0
                        freq = 0
                    else:
                        L = 0
                        freq = 0

                        if cell_type == 'Mid Cell':
                            L_i, L_o = L, L
                        else:
                            L_o = L

                alpha_i, error_msg1 = calculate_alpha(A_i, B_i, a_i, b_i, Ri_i, L_i, Req, 0)
                alpha_o, error_msg2 = calculate_alpha(A_o, B_o, a_o, b_o, Ri_o, L_o, Req, 0)

                inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req, alpha_i]
                outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req, alpha_o]
            else:

                # new mid cell and end cell with initial Req guess
                if cell_type == 'End Cell':  # change later for two types of end cells,
                    # one with half mid cell and the other with same dimensions
                    Req_o = Req

                if cell_type == 'End-End Cell':
                    if Req < Ri_i + B_i + b_i or Req < Ri_o + B_o + b_o:
                        Req = Ri_i + B_i + b_i
                    Req_o = Req

                inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req, 0]
                outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req_o, 0]

                # edit to check for key later
                if key not in existing_keys:
                    if tuner_option == 'SLANS' and slans_tune:
                        if tune_variable == 'Req':
                            # Tune cell to get Req
                            Req, freq, alpha, h, e = slans_tune.mid_cell_tune(A_i, B_i, a_i, b_i, Ri_i, L_i, Req,
                                                                              target_freq, proc=proc)
                        else:
                            L, freq, alpha, h, e = slans_tune.end_cell_tune(inner_cell, outer_cell,
                                                                            target_freq, proc=proc)

                            if cell_type == 'Mid Cell':
                                L_i, L_o = L, L
                            else:
                                L_o = L

                        inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req, alpha]
                        outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req, alpha]
                    else:
                        if tune_variable == 'Req':
                            try:
                                Req, freq = pytune.tuneR(inner_cell, outer_cell, target_freq, beampipes, bc,
                                                         parentDir, projectDir, iter_set=iter_set, proc=proc,
                                                         conv_list=convergence_list)
                            except FileNotFoundError:
                                Req, freq = 0, 0
                            # round
                            Req, freq = round(Req, 4), round(freq, 2)
                        else:
                            try:
                                L, freq = pytune.tuneL(inner_cell, outer_cell, target_freq, beampipes, bc,
                                                       parentDir, projectDir, iter_set=iter_set, proc=proc,
                                                       conv_list=convergence_list)
                            except FileNotFoundError:
                                L, freq = 0, 0

                            if cell_type == 'Mid Cell':
                                L_i, L_o = L, L
                            else:
                                L_o = L

                        if Req == 0 or L_i == 0 or L_o == 0:
                            alpha_i = 0
                            alpha_o = 0
                        else:
                            alpha_i, error_msg1 = calculate_alpha(A_i, B_i, a_i, b_i, Ri_i, L_i, Req, 0)
                            alpha_o, error_msg2 = calculate_alpha(A_o, B_o, a_o, b_o, Ri_o, L_o, Req, 0)

                        inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req, alpha_i]
                        outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req, alpha_o]

            result = "Failed"
            # # round
            # L, freq = round(L, 2), round(freq, 2)
            if (1 - 0.001)*target_freq < round(freq, 2) < (1 + 0.001)*target_freq \
                    and (90.5 <= alpha_i <= 180) \
                    and (90.5 <= alpha_o <= 180) and error_msg1 == 1 and error_msg2 == 1:
                result = f"Success: {target_freq, freq}"

                if cell_type == 'Mid Cell':
                    population[key] = {"IC": inner_cell, "OC": outer_cell, "BP": 'none', 'FREQ': freq}
                else:
                    population[key] = {"IC": inner_cell, "OC": outer_cell, "BP": 'both', 'FREQ': freq}

                # save last slans run if activated. Save only when condition is fulfilled
                if save_last:
                    # self.save_last(projectDir, proc, key)
                    # print("n_cell tuner", n_cell_last_run)
                    if 'End' in cell_type:
                        beampipes = 'both'

                    slans_geom.cavity(n_cell_last_run, 1, inner_cell, outer_cell, outer_cell, f_shift=0, bc=bc,
                                      beampipes=beampipes, proc=proc,
                                      n_modes=n_cell_last_run+1, fid=key, parentDir=parentDir, projectDir=projectDir)

                    # write cst_studio parameters
                    write_cst_paramters(key, inner_cell, outer_cell, projectDir, cell_type)

                    # write tune results
                    if cell_type == 'Mid Cell':
                        d_tune_res = {'Req': Req, 'L': L_i, 'alpha_i': alpha_i, 'alpha_o': alpha_o, 'freq': freq}
                    else:
                        d_tune_res = {'Req': Req, 'L': L_o, 'alpha_i': alpha_i, 'alpha_o': alpha_o, 'freq': freq}
                    self.save_tune_result(d_tune_res, projectDir, key)

            print(f'Done Tuning Cavity {key}: {result}')

            # clear folder after every run. This is to avoid copying of wrong values to save folder
            # processor folder
            proc_fold = fr'{projectDir}/SimulationData/SLANS/Cavity_process_{proc}'
            keep = ['SLANS_exe']
            for item in os.listdir(proc_fold):
                if item not in keep:  # If it isn't in the list for retaining
                    try:
                        os.remove(item)
                    except FileNotFoundError:
                        continue

            # update progress
            progress_list.append((progress + 1) / total_no_of_shapes)

            # Update progressbar
            progress += 1

            # print("Saving Dictionary", f"shape_space{proc}.json")◙
            # print("Done saving")

        end = time.time()

        runtime = end - start
        print(f'Processor {proc} runtime: {runtime}')

    # def calculate_alpha(self, A, B, a, b, Ri, L, Req, L_bp):
    #
    #     data = ([0 + L_bp, Ri + b, L + L_bp, Req - B],
    #             [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    #     df = fsolve(self.ellipse_tangent,
    #                             np.array([0.5*a + L_bp, Ri + 0.5 * b, L - A + L_bp, Req - 0.5 * B]),
    #                             args=data, full_output=True)
    #     x1, y1, x2, y2 = df[0]
    #     error_msg = df[-2]
    #
    #     m = (y2 - y1) / (x2 - x1)
    #     alpha = 180 - np.arctan(m) * 180 / np.pi
    #     return alpha, error_msg

    # @staticmethod
    # def save_last(projectDir, key):
    #     # copy slans files from folder with key
    #     # Source path
    #     # src = fr'{projectDir}/SimulationData/SLANS/Cavity_process_{process_num}'
    #     # Destination path
    #     dest = fr'{projectDir}/SimulationData/SLANS/Cavity{key}'
    #     if os.path.exists(dest):
    #         shutil.rmtree(dest)
    #     # Copy the content of source to destination
    #     # destination = shutil.copytree(src, dest)
    #
    #     # delete SLANS_EXE folder
    #     if os.path.exists(fr'{projectDir}/SimulationData/SLANS/Cavity{key}/SLANS_exe'):
    #         shutil.rmtree(fr'{projectDir}/SimulationData/SLANS/Cavity{key}/SLANS_exe')

    @staticmethod
    def save_tune_result(d, projectDir, key):
        with open(fr"{projectDir}\SimulationData\SLANS\Cavity{key}\tune_res.json", 'w') as file:
            file.write(json.dumps(d, indent=4, separators=(',', ': ')))

    # if __name__ == '__main__':
#     #
#     tune = Tuner()
#
#     tune_var =
#     par_mid =
#     par_end =
#     target_freq = 400  # MHz
#     beampipes =
#     bc =  # boundary conditions
#     parentDir = ""  # location of slans code. See folder structure in the function above
#     projectDir = ""  # location to write results to
#     iter_set =
#     proc = 0
#     tune.tune(self, tune_var, par_mid, par_end, target_freq, beampipes, bc, parentDir, projectDir, iter_set, proc=0):
