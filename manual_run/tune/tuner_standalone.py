import sys
import time
import pandas as pd
from icecream import ic
from scipy.optimize import fsolve
from pyTuner_standalone import PyTune
import json
import os
import numpy as np


class Tuner:
    def __init__(self):
        pass

    def tune(self, pseudo_shape_space, bc, parentDir, projectDir, filename, resume="No",
             proc=0, tune_variable='Req', iter_set=None, cell_type='Mid Cell'):
        """

        Parameters
        ----------
        pseudo_shape_space:
        bc: int {33, 22}
            Boundary condition: 33->magnetic, 22->electric
        parentDir: str, path
            Directory containing SLANS exe files
        projectDir: str, path
            Directory to write results to
        filename: str
            Filename
        resume: str {'Yes', 'No'}
            If 'Yes', previously analysed geometries are not rerun. Previously analysed geometries are determined using
            filenames
        proc: 0
        tune_variable: str {'Req', 'L'}
            Use 'Req' to tune the equator radius for mid cells and 'L' to tune the half length for end cells
        iter_set: list [<Interpolation type>, <tolerance>, <maximum number of iterations>]
            Number of iterations to run
        cell_type: str {'Mid Cell', 'End Cell'}
            Type of cell to be tuned

        Returns
        -------

        """

        # tuner
        pytune = PyTune()

        start = time.time()
        population = {}
        total_no_of_shapes = len(list(pseudo_shape_space.keys()))

        # check for already processed shapes
        existing_keys = []

        if resume == "Yes":
            # check if value set is already written. This is to enable continuation in case of break in program
            if os.path.exists(fr'{projectDir}{filename}'):
                population = json.load(open(fr'{projectDir}/{filename}', 'r'))

                existing_keys = list(population.keys())

        for key, pseudo_shape in pseudo_shape_space.items():
            A_i, B_i, a_i, b_i, Ri_i, L_i, Req = pseudo_shape['IC'][:7]
            A_o, B_o, a_o, b_o, Ri_o, L_o, Req_o = pseudo_shape['OC'][:7]

            beampipes = pseudo_shape['BP']
            target_freq = pseudo_shape['FREQ']

            # Check if simulation is already run
            freq = 0
            alpha_i = 0
            alpha_o = 0
            if resume == "Yes" and os.path.exists(fr'{projectDir}/{key}'):
                # if folder exist, read value
                filename = fr'{projectDir}/{key}/cavity_{bc}.svl'
                try:
                    data_dict = svl_reader(filename)
                    # print(data_dict)
                    if tune_variable == 'Req':
                        Req = data_dict['CAVITY RADIUS'][0] * 10
                        freq = data_dict['FREQUENCY'][0]
                    else:
                        L = data_dict['LENGTH'][0] * 10
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
                    if tune_variable == 'Req':
                        try:
                            Req, freq = pytune.tuneR(inner_cell, outer_cell, target_freq, beampipes, bc,
                                                     parentDir, projectDir, iter_set=iter_set, fid=key)
                        except FileNotFoundError as e:
                            Req, freq = 0, 0
                            print('File not found error:: ', e)
                        # round
                        Req, freq = Req, freq
                    else:
                        try:
                            L, freq = pytune.tuneL(inner_cell, outer_cell, target_freq, beampipes, bc,
                                                   parentDir, projectDir, iter_set=iter_set, fid=key)
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
            if (1 - 0.001) * target_freq < round(freq, 2) < (1 + 0.001) * target_freq \
                    and (90.0 <= alpha_i <= 180) \
                    and (90.0 <= alpha_o <= 180) and error_msg1 == 1 and error_msg2 == 1:
                result = f"Success: {target_freq, freq}"

                if cell_type == 'Mid Cell':
                    population[key] = {"IC": inner_cell, "OC": outer_cell, "BP": 'none', 'FREQ': freq}
                else:
                    population[key] = {"IC": inner_cell, "OC": outer_cell, "BP": 'both', 'FREQ': freq}

            # write tune results
            if cell_type == 'Mid Cell':
                d_tune_res = {'Req': Req, 'L': L_i, 'alpha_i': alpha_i, 'alpha_o': alpha_o, 'freq': freq}
            else:
                d_tune_res = {'Req': Req, 'L': L_o, 'alpha_i': alpha_i, 'alpha_o': alpha_o, 'freq': freq}

            self.save_tune_result(d_tune_res, projectDir, key)
            print(f'Done Tuning Cavity {key}: {result}')

            # clear folder after every run. This is to avoid copying of wrong values to save folder
            # processor folder
            proc_fold = fr'{projectDir}/_process_{proc}'
            keep = ['SLANS_exe']
            for item in os.listdir(proc_fold):
                if item not in keep:  # If it isn't in the list for retaining
                    try:
                        os.remove(item)
                    except FileNotFoundError:
                        continue

        end = time.time()

        runtime = end - start
        print(f'Processor {proc} runtime: {runtime}')

    @staticmethod
    def save_tune_result(d, projectDir, key):
        with open(fr"{projectDir}\{key}\tune_res.json", 'w') as file:
            file.write(json.dumps(d, indent=4, separators=(',', ': ')))


def update_alpha(cell):
    """
    Update geometry json file variables to include the value of alpha

    Parameters
    ----------
    cell:
        Cavity geometry parameters

    Returns
    -------
    List of cavity geometry parameters

    """
    if len(cell) == 8:
        A, B, a, b, Ri, L, Req = cell[:7]
    else:
        A, B, a, b, Ri, L, Req = cell[:7]
    alpha = calculate_alpha(A, B, a, b, Ri, L, Req, 0)
    cell = [A, B, a, b, Ri, L, Req, alpha[0]]

    return cell


def calculate_alpha(A, B, a, b, Ri, L, Req, L_bp):
    """
    Calculates the largest angle the tangent line of two ellipses makes with the horizontal axis

    Parameters
    ----------
    A: float
    B: float
    a: float
    b: float
    Ri: float
    L: float
    Req: float
    L_bp: float

    Returns
    -------
    alpha: float
        Largest angle the tangent line of two ellipses makes with the horizontal axis
    error_msg: int
        State of the iteration, failed or successful. Refer to

    """

    df = tangent_coords(A, B, a, b, Ri, L, Req, L_bp)
    x1, y1, x2, y2 = df[0]
    error_msg = df[-2]

    alpha = 180 - np.arctan2(y2 - y1, (x2 - x1)) * 180 / np.pi
    return alpha, error_msg


def tangent_coords(A, B, a, b, Ri, L, Req, L_bp, tangent_check=False):
    """
    Calls to :py:func:`utils.shared_function.ellipse_tangent`
    Parameters
    ----------
    A: float
        Equator ellipse dimension
    B: float
        Equator ellipse dimension
    a: float
        Iris ellipse dimension
    b: float
        Iris ellipse dimension
    Ri: float
        Iris radius
    L: float
        Cavity half cell length
    Req: float
        Cavity equator radius
    L_bp: float
        Cavity beampipe length
    tangent_check: bool
        If set to True, the calculated tangent line as well as the ellipses are plotted and shown

    Returns
    -------
    df: pandas.Dataframe
        Pandas dataframe containing information on the results from fsolve
    """
    # data = ([0 + L_bp, Ri + b, L + L_bp, Req - B],
    #         [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    #
    # df = fsolve(ellipse_tangent,
    #             np.array([a + L_bp, Ri + f[0] * b, L - A + L_bp, Req - f[1] * B]),
    #             args=data, fprime=jac, xtol=1.49012e-12, full_output=True)
    #
    #     # ic(df)

    data = ([0 + L_bp, Ri + b, L + L_bp, Req - B], [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    checks = {"non-reentrant": [0.45, -0.45],
              "reentrant": [0.85, -0.85],
              "expansion": [0.15, -0.01]}

    for key, ch in checks.items():
        # check for non-reentrant cavity
        df = fsolve(ellipse_tangent,
                    np.array([a + L_bp, Ri + ch[0] * b, L - A + L_bp, Req + ch[1] * B]),
                    args=data, fprime=jac, xtol=1.49012e-12, full_output=True)
        # ic(df)
        x1, y1, x2, y2 = df[0]
        alpha = 180 - np.arctan2(y2 - y1, (x2 - x1)) * 180 / np.pi

        if key == 'non-reentrant':
            if 90 < alpha < 180:
                return df
        elif key == 'reentrant':
            if 0 < alpha < 90:
                return df
        elif key == 'expansion':
            if 90 < alpha < 180:
                return df

    # error_msg = df[-2]
    # if error_msg == 4:
    #     df = fsolve(ellipse_tangent, np.array([a + L_bp, Ri + 1.15 * b, L - A + L_bp, Req - 1.15 * B]),
    #                 args=data, fprime=jac, xtol=1.49012e-12, full_output=True)

    return df


def ellipse_tangent(z, *data):
    """
    Calculates the coordinates of the tangent line that connects two ellipses

    .. _ellipse tangent:

    .. figure:: ../images/ellipse_tangent.png
       :alt: ellipse tangent
       :align: center
       :width: 200px

    Parameters
    ----------
    z: list, array like
        Contains list of tangent points coordinate's variables ``[x1, y1, x2, y2]``.
        See :numref:`ellipse tangent`
    data: list, array like
        Contains midpoint coordinates of the two ellipses and the dimensions of the ellipses
        data = ``[coords, dim]``; ``coords`` = ``[h, k, p, q]``, ``dim`` = ``[a, b, A, B]``


    Returns
    -------
    list of four non-linear functions

    Note
    -----
    The four returned non-linear functions are

    .. math::

       f_1 = \\frac{A^2b^2(x_1 - h)(y_2-q)}{a^2B^2(x_2-p)(y_1-k)} - 1

       f_2 = \\frac{(x_1 - h)^2}{a^2} + \\frac{(y_1-k)^2}{b^2} - 1

       f_3 = \\frac{(x_2 - p)^2}{A^2} + \\frac{(y_2-q)^2}{B^2} - 1

       f_4 = \\frac{-b^2(x_1-x_2)(x_1-h)}{a^2(y_1-y_2)(y_1-k)} - 1
    """

    coord, dim = data
    h, k, p, q = coord
    a, b, A, B = dim
    x1, y1, x2, y2 = z

    f1 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
    f2 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
    f3 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
    f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

    return f1, f2, f3, f4


def jac(z, *data):
    """
    Computes the Jacobian of the non-linear system of ellipse tangent equations

    Parameters
    ----------
    z: list, array like
        Contains list of tangent points coordinate's variables ``[x1, y1, x2, y2]``.
        See :numref:`ellipse tangent`
    data: list, array like
        Contains midpoint coordinates of the two ellipses and the dimensions of the ellipses
        data = ``[coords, dim]``; ``coords`` = ``[h, k, p, q]``, ``dim`` = ``[a, b, A, B]``

    Returns
    -------
    J: array like
        Array of the Jacobian

    """
    coord, dim = data
    h, k, p, q = coord
    a, b, A, B = dim
    x1, y1, x2, y2 = z

    # f1 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
    # f2 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
    # f3 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
    # f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

    df1_dx1 = A ** 2 * b ** 2 * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k))
    df1_dy1 = - A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k) ** 2)
    df1_dx2 = - A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) ** 2 * (y1 - k))
    df1_dy2 = A ** 2 * b ** 2 * (x1 - h) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k))

    df2_dx1 = 2 * (x1 - h) / a ** 2
    df2_dy1 = 2 * (y1 - k) / b ** 2
    df2_dx2 = 0
    df2_dy2 = 0

    df3_dx1 = 0
    df3_dy1 = 0
    df3_dx2 = 2 * (x2 - p) / A ** 2
    df3_dy2 = 2 * (y2 - q) / B ** 2

    df4_dx1 = -b ** 2 * ((x1 - x2) + (x1 - h)) / (a ** 2 * (y1 - y2) * (y1 - k))
    df4_dy1 = -b ** 2 * (x1 - x2) * (x1 - h) * ((y1 - y2) + (y1 - k)) / (a ** 2 * ((y1 - y2) * (y1 - k)) ** 2)
    df4_dx2 = b ** 2 * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k))
    df4_dy2 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) ** 2 * (y1 - k))

    J = [[df1_dx1, df1_dy1, df1_dx2, df1_dy2],
         [df2_dx1, df2_dy1, df2_dx2, df2_dy2],
         [df3_dx1, df3_dy1, df3_dx2, df3_dy2],
         [df4_dx1, df4_dy1, df4_dx2, df4_dy2]]

    return J


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


if __name__ == '__main__':
    tune = Tuner()

    beampipes = 'none'  # {'none', 'both', 'left', 'right'}
    bc = 33  # boundary conditions
    parentDir = r"D:\Dropbox\CavityDesignHub\manual_run\slans\SLANS_exe"  # location of slans exe files
    projectDir = r"D:\Dropbox\CavityDesignHub\manual_run\tune\Tests"  # location to write results to
    iter_set = ["Linear Interpolation", 1e-3, 10]  # accuracy cannot be better than 1e-3

    # load shape space
    data_A_b = pd.read_csv(r'D:\Dropbox\piotr\scee\initial_points_Sos.dat', names=['A', 'B', 'a', 'b '],
                           header=None, sep='\s+')  #I
    data_Ri_Req = pd.DataFrame(index=range(len(data_A_b)), columns=['Ri', 'L', 'Req'])
    data_Ri_Req.loc[:, ['Ri', 'L', 'Req']] = [35, 57.7, 122.1067]
    df = pd.concat([data_A_b, data_Ri_Req], axis=1)

    pseudo_shape_space = {}
    # pseudo_shape_space = {"f_opt": {
    #     "IC": [10.16698881, 58.71482539, 6.26252505, 27.0240481, 35, 57.7, 122.1067],  # mid cell [A, B, a, b, Ri, L, Req]
    #     "OC": [10.16698881, 58.71482539, 6.26252505, 27.0240481, 35, 57.7, 122.1067],  # left end cell
    #     "OC_R": [10.16698881, 58.71482539, 6.26252505, 27.0240481, 35, 57.7, 122.1067],  # right end cell
    #     "BP": "none",  # beampipe
    #     "FREQ": 1300  # tune frequency
    # }}  # format
    for key, row in df.iterrows():
        pseudo_shape_space[key] = {}
        pseudo_shape_space[key]['IC'] = list(row)
        pseudo_shape_space[key]['OC'] = list(row)
        pseudo_shape_space[key]['OC_R'] = list(row)
        pseudo_shape_space[key]["BP"] = "none"
        pseudo_shape_space[key]["FREQ"] = 1300

    ic(pseudo_shape_space)
    filename = 'Test'
    tune.tune(pseudo_shape_space, bc, parentDir, projectDir, filename, resume="No",
              proc=0, tune_variable='Req', iter_set=iter_set, cell_type='Mid Cell')
