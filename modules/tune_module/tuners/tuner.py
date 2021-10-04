import json
import os
import time
import numpy as np
from scipy.optimize import fsolve
from modules.tune_module.tuners.pyTuner import PyTune
from simulation_codes.SLANS.slansTuner import SLANSTune


class Tuner:
    def __init__(self):
        pass

    def tune(self, pseudo_shape_space, bc, parentDir, projectDir, filename, resume="No",
              proc=0, tuner_option='SLANS', tune_variable='Req', iter_set=None, cell_type='Mid Cell', progress_list=None):

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
                print(f'Existing keys: {existing_keys}')

        progress = 0

        for key, pseudo_shape in pseudo_shape_space.items():

            A_i, B_i, a_i, b_i, Ri_i, L_i, Req, _ = pseudo_shape['IC'] # Req here is none but required since the shape space consists of all variables
            A_o, B_o, a_o, b_o, Ri_o, L_o, Req, _ = pseudo_shape['OC'] # Req here is none but required since the shape space consists of all variables
            beampipes = pseudo_shape['BP']
            freq = pseudo_shape['FREQ']

            # new mid cell and end cell with initial Req guess
            inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req, 0]
            outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req, 0]

            # edit to check for key later
            if key not in existing_keys:
                if tuner_option == 'SLANS' and slans_tune:
                    if tune_variable == 'Req':
                        # Tune cell to get Req
                        Req, freq, alpha, h, e = slans_tune.mid_cell_tune(A_i, B_i, a_i, b_i, Ri_i, L_i, Req, freq, proc=proc)
                    else:
                        L, freq, alpha, h, e = slans_tune.end_cell_tune(inner_cell, outer_cell, freq, proc=proc)

                        if cell_type == 'Mid Cell':
                            L_i, L_o = L, L
                        else:
                            L_o = L

                    inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req, alpha]
                    outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req, alpha]
                else:
                    if tune_variable == 'Req':
                        Req, freq = pytune.tuneR(inner_cell, outer_cell, freq, beampipes, bc, parentDir, projectDir, iter_set=iter_set, proc=proc)
                        # round
                        Req, freq = round(Req, 4), round(freq, 2)
                    else:
                        L, freq = pytune.tuneL(inner_cell, outer_cell, freq, beampipes, bc, parentDir, projectDir, iter_set=iter_set, proc=proc)

                        # round
                        L, freq = round(L, 2), round(freq, 2)

                        if cell_type == 'Mid Cell':
                            L_i, L_o = L, L
                        else:
                            L_o = L

                    alpha_i = self.calculate_alpha(A_i, B_i, a_i, b_i, Ri_i, L_i, Req, 0)
                    alpha_o = self.calculate_alpha(A_o, B_o, a_o, b_o, Ri_o, L_o, Req, 0)

                    inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req, alpha_i]
                    outer_cell = [A_o, B_o, a_o, b_o, Ri_o, L_o, Req, alpha_o]

                print(f'Done Tuning Cavity {key}')

                # update progress
                progress_list.append((progress + 1) / total_no_of_shapes)

                if cell_type == 'Mid Cell':
                    population[key] = {"IC": inner_cell, "OC": outer_cell, "BP": 'none', 'FREQ': freq}
                else:
                    population[key] = {"IC": inner_cell, "OC": outer_cell, "BP": 'both', 'FREQ': freq}

            # Update progressbar
            progress += 1

            # print("Saving Dictionary", f"shape_space{proc}.json")
            with open(fr"{projectDir}\Cavities\shape_space{proc}.json", 'w') as file:
                file.write(json.dumps(population, indent=4, separators=(',', ': ')))
            # print("Done saving")

        end = time.time()

        runtime = end - start
        print(f'Processor {proc} runtime: {runtime}')

    def calculate_alpha(self, A, B, a, b, Ri, L, Req, L_bp):

        data = ([0 + L_bp, Ri + b, L + L_bp, Req - B],
                [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
        x1, y1, x2, y2 = fsolve(self.ellipse_tangent,
                                np.array([a + L_bp, Ri + 0.85 * b, L - A + L_bp, Req - 0.85 * B]),
                                args=data)
        m = (y2 - y1) / (x2 - x1)
        alpha = 180 - np.arctan(m) * 180 / np.pi
        return alpha

    def ellipse_tangent(self, z, *data):
        coord, dim = data
        h, k, p, q = coord
        a, b, A, B = dim
        x1, y1, x2, y2 = z

        f1 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
        f2 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
        f3 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
        f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

        return f1, f2, f3, f4


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
#
#     tune.tune(self, tune_var, par_mid, par_end, target_freq, beampipes, bc, parentDir, projectDir, iter_set, proc=0):