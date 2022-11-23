import ast

import quadpy
from icecream import ic
from scipy.optimize import fsolve
import numpy as np


def calculate_alpha(A, B, a, b, Ri, L, Req, L_bp):
    # data = ([0 + L_bp, Ri + b, L + L_bp, Req - B],
    #         [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    # x1, y1, x2, y2 = fsolve(ellipse_tangent,
    #                         np.array([a + L_bp, Ri + 0.85 * b, L - A + L_bp, Req - 0.85 * B]),
    #                         args=data)
    # m = (y2 - y1) / (x2 - x1)
    # alpha = 180 - np.arctan(m) * 180 / np.pi
    # return alpha

    data = ([0 + L_bp, Ri + b, L + L_bp, Req - B],
            [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    df = fsolve(ellipse_tangent,
                np.array([a + L_bp, Ri + 0.85 * b, L - A + L_bp, Req - 0.85 * B]),
                args=data, full_output=True)
    x1, y1, x2, y2 = df[0]
    error_msg = df[-2]

    m = (y2 - y1) / (x2 - x1)
    alpha = 180 - np.arctan(m) * 180 / np.pi
    return alpha, error_msg


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


def perform_geometry_checks(par_mid, par_end):
    # # check if Req is less than lower limit
    # if par_mid[6] < par_mid[1] + par_mid[3] + par_mid[4] or par_end[6] < par_end[1] + par_end[3] + par_end[4]:
    #     return False

    # check if alpha is less than 90.5
    if par_mid == par_end:
        alpha, error_msg = calculate_alpha(par_mid[0], par_mid[1], par_mid[2], par_mid[3], par_mid[4], par_mid[5],
                                           par_mid[6], 0)
        if alpha < 90.0 or error_msg != 1:
            print("1:", alpha, error_msg)
            return False
    else:
        alpha, error_msg = calculate_alpha(par_end[0], par_end[1], par_end[2], par_end[3], par_end[4], par_end[5],
                                           par_mid[6], 0)
        if alpha < 90.0 or error_msg != 1:
            print("2:", alpha, error_msg)
            return False

        alpha, error_msg = calculate_alpha(par_mid[0], par_mid[1], par_mid[2], par_mid[3], par_mid[4], par_mid[5],
                                           par_mid[6], 0)
        if alpha < 90.0 or error_msg != 1:
            print("3:", alpha, error_msg)
            return False

    # check if L is less than lower limit
    if par_mid[5] < par_mid[0] or par_end[5] < par_end[0]:
        print("4:", alpha, error_msg)
        return False

    return True


def write_cst_paramters(key, ic_, oc, projectDir, cell_type):
    if cell_type is None:
        # print("Writing parameters to file")
        path = fr'{projectDir}/SimulationData/SLANS/Cavity{key}/{key}.txt'

        # print(path)
        with open(path, 'w') as f:
            name_list = ['Aeq', 'Beq', 'ai', 'bi', 'Ri', 'L', 'Req', 'alpha', 'Aeq_e', 'Beq_e', 'ai_e', 'bi_e', 'Ri_e',
                         'L_e', 'Req', 'alpha_e', 'key']

            value_list = [ic_[0], ic_[1], ic_[2], ic_[3], ic_[4], ic_[5], ic_[6], ic_[7],
                          oc[0], oc[1], oc[2], oc[3], oc[4], oc[5], oc[6], oc[7], key]

            for i in range(len(name_list)):
                if name_list[i] == 'key':
                    f.write(f'{name_list[i]} = "{0}" "{value_list[i]}"\n')
                else:
                    f.write(f'{name_list[i]} = "{value_list[i]}" ""\n')

    else:
        # print("Writing parameters to file")
        path = fr'{projectDir}/SimulationData/SLANS/Cavity{key}/{key}.txt'
        path_mc = fr'{projectDir}/SimulationData/SLANS/Cavity{key}/{key}_Multicell.txt'

        # print(path)
        with open(path, 'w') as f:
            name_list = ['Aeq', 'Beq', 'ai', 'bi', 'Ri', 'L', 'Req', 'alpha', 'Aeq_e', 'Beq_e', 'ai_e', 'bi_e', 'Ri_e',
                         'L_e', 'Req_e', 'alpha_e', 'key']

            if cell_type == 'Mid Cell':
                value_list = [ic_[0], ic_[1], ic_[2], ic_[3], ic_[4], ic_[5], ic_[6], ic_[7],
                              'Aeq', 'Beq', 'ai', 'bi', 'Ri', 'L', 'Req', 'alpha', key]
            else:
                value_list = [ic_[0], ic_[1], ic_[2], ic_[3], ic_[4], ic_[5], ic_[6], ic_[7],
                              oc[0], oc[1], oc[2], oc[3], oc[4], oc[5], oc[6], oc[7], key]

            for i in range(len(name_list)):
                if name_list[i] == 'key':
                    f.write(f'{name_list[i]} = "{0}" "{value_list[i]}"\n')
                else:
                    f.write(f'{name_list[i]} = "{value_list[i]}" ""\n')

        with open(path_mc, 'w') as f:
            name_list = ['Aeq', 'Beq', 'ai', 'bi', 'Ri', 'L', 'Req', 'alpha',
                         'Aeq_er', 'Beq_er', 'ai_er', 'bi_er', 'Ri_er', 'L_er', 'Req_er', 'alpha_er',
                         'Aeq_el', 'Beq_el', 'ai_el', 'bi_el', 'Ri_el', 'L_el', 'Req_el', 'alpha_el', 'key']

            if cell_type == 'Mid Cell':
                value_list = [ic_[0], ic_[1], ic_[2], ic_[3], ic_[4], ic_[5], ic_[6], ic_[7],
                              'Aeq', 'Beq', 'ai', 'bi', 'Ri', 'L', 'Req', 'alpha',
                              'Aeq', 'Beq', 'ai', 'bi', 'Ri', 'L', 'Req', 'alpha',
                              key]
            else:
                value_list = [ic_[0], ic_[1], ic_[2], ic_[3], ic_[4], ic_[5], ic_[6], ic_[7],
                              oc[0], oc[1], oc[2], oc[3], oc[4], oc[5], oc[6], oc[7],
                              'Aeq_er', 'Beq_er', 'ai_er', 'bi_er', 'Ri_er', 'L_er', 'Req_er', 'alpha_er',
                              key]

            for i in range(len(name_list)):
                if name_list[i] == 'key':
                    f.write(f'{name_list[i]} = "{0}" "{value_list[i]}"\n')
                else:
                    f.write(f'{name_list[i]} = "{value_list[i]}" ""\n')

        # print("Writing to file complete.")


def stroud(p):
    # Stroud-3 method
    #
    # Input parameters:
    #  p   number of dimensions
    # Output parameters:
    #  nodes   nodes of quadrature rule in [0,1]**p (column-wise)
    #

    nodes = np.zeros((p, 2 * p))
    coeff = np.pi / p
    fac = np.sqrt(2 / 3)

    for i in range(2 * p):
        for r in range(int(np.floor(0.5 * p))):
            k = 2 * r
            nodes[k, i] = fac * np.cos((k + 1) * (i + 1) * coeff)
            nodes[k + 1, i] = fac * np.sin((k + 1) * (i + 1) * coeff)

        if 0.5 * p != np.floor(0.5 * p):
            nodes[-1, i] = ((-1) ** (i + 1)) / np.sqrt(3)

    # transform nodes from [-1,+1]**p to [0,1]**p
    nodes = 0.5 * nodes + 0.5

    return nodes


def quad_stroud3(rdim, degree):
    # data for Stroud-3 quadrature in [0,1]**k
    # nodes and weights
    nodes = stroud(rdim)
    nodestr = 2. * nodes - 1.
    weights = (1 / (2 * rdim)) * np.ones((2 * rdim, 1))

    # evaluation of Legendre polynomials
    bpoly = np.zeros((degree + 1, rdim, 2 * rdim))
    for ll in range(rdim):
        for j in range(2 * rdim):
            bpoly[0, ll, j] = 1
            bpoly[1, ll, j] = nodestr[ll, j]
            for i in range(1, degree):
                bpoly[i + 1, ll, j] = ((2 * (i + 1) - 1) * nodestr[ll, j] * bpoly[i, ll, j] - i * bpoly[
                    i - 1, ll, j]) / (i + 1)

    # standardisation of Legendre polynomials
    for i in range(1, degree + 1):
        bpoly[i, :, :] = bpoly[i, :, :] * np.sqrt(2 * (i + 1) - 1)

    return nodes, weights, bpoly


def quad_stroud5(rdim, degree):
    # data for Stroud-5 quadrature in [0,1]**rdim
    # nodes and weights
    o, nodestr, weights = quadpy.cn.stroud_cn_5_2(rdim)
    nodes = 0.5*nodestr + 0.5
    weights = weights/(2**rdim)
    dummy, nnodes = np.size(nodes)
    
    # evaluation of Legendre polynomials
    bpoly = np.zeros((degree+1, rdim, nnodes))
    for l in range(rdim):
        for j in range(nnodes):
            bpoly[1, l, j] = 1.
            bpoly[2, l, j] = nodestr(l,j)
            for i in range(1, degree):
                bpoly[i+1, l, j] = ((2*i-1)*nodestr(l,j)*bpoly[i,l,j]-(i-1)*bpoly[i-1,l,j])/i
    
    # standardisation of Legendre polynomials
    for i in range(1, degree+1):
        bpoly[i, :, :] = bpoly[i, :, :]*np.sqrt(2*i-1)

    return nodes, weights, bpoly


def weighted_mean_obj(tab_var, weights):
    rows_sims_no, cols = np.shape(tab_var)
    no_weights, dummy = np.shape(weights)  # z funckji quadr_stroud wekt columnowy

    if rows_sims_no == no_weights:
        expe = np.zeros((cols, 1))
        outvar = np.zeros((cols, 1))
        for i in range(cols):
            expe[i, 0] = np.dot(tab_var[:, i], weights)
            outvar[i, 0] = np.dot(tab_var[:, i] ** 2, weights)

        stdDev = np.sqrt(abs(outvar - expe ** 2))
    else:
        expe = 0
        stdDev = 0
        ic('Cols_sims_no != No_weights')

    return list(expe.T[0]), list(stdDev.T[0])


def normal_dist(x, mean, sd):
    prob_density = (np.pi*sd) * np.exp(-0.5*((x-mean)/sd)**2)
    return prob_density


def text_to_list(ll):
    if ll == '':
        return None
    else:
        ll = ast.literal_eval(ll)
        if isinstance(ll, int) or isinstance(ll, float):
            return [ll, 2e10]
        else:
            return list(ll)


def text_to_list_(ll):
    if ll == '':
        return None
    else:
        ll = ast.literal_eval(ll)
        if isinstance(ll, int) or isinstance(ll, float):
            return [ll]
        else:
            return list(ll)
