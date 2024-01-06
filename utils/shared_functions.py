import ast
import copy
import json
import math
import os
# from datetime import time
import pandas as pd
from PyQt5.QtGui import QDoubleValidator, QValidator
from PyQt5.QtWidgets import *
from icecream import ic
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from qtpy import QtWidgets
from scipy.linalg import eigh
from scipy.optimize import fsolve, minimize_scalar
import numpy as np
from PyQt5.QtCore import *
from utils.file_reader import FileReader
# from analysis_modules.eigenmode.SLANS.slans_geometry import SLANSGeometry

fr = FileReader()
# slans_geom = SLANSGeometry()

AN_DURATION = 250
global animate


def update_alpha(cell, cell_type='normal'):
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
    A, B, a, b, Ri, L, Req = cell[:7]
    alpha = calculate_alpha(A, B, a, b, Ri, L, Req, 0)
    if cell_type == 'normal':
        cell = [A, B, a, b, Ri, L, Req, alpha[0]]
    else:
        cell = [A, B, a, b, Ri, L, Req, cell[7], alpha[0]]

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

        # ic(data)
        if tangent_check:
            h, k, p, q = data[0]
            a, b, A, B = data[1]
            el_ab = Ellipse((h, k), 2 * a, 2 * b, alpha=0.5)
            el_AB = Ellipse((p, q), 2 * A, 2 * B, alpha=0.5)

            fig, ax = plt.subplots()
            ax.add_artist(el_ab)
            ax.add_artist(el_AB)

            x1, y1, x2, y2 = df[0]
            # ax.set_xlim(-1.1 * a, 1.1 * L)
            # ax.set_ylim(Ri, 1.1 * Req)
            ax.plot([x1, x2], [y1, y2])

            plt.show()

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


def perform_geometry_checks(par_mid: list, par_end: list):
    """
    Checks geometry to filter out low loss cavity geometries

    Parameters
    ----------
    par_mid: list, ndarray
        Mid cell parameters of cavity
    par_end: list, ndarray
        End cell parameters of cavity

    Returns
    -------
    bool

    """
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


def write_cst_paramters(key, ic_, oc_l, oc_r, projectDir, cell_type, opt=False, solver='slans'):
    """
    Writes cavity geometric data that can be imported into CST Studio

    Parameters
    ----------
    key: str, int
        Cavity marker
    ic_: list, array like
        Inner cavity cell geometric variables
    oc_l: list, array like
        Outer cavity cell geometric variables
    projectDir: str
        Project directory
    cell_type: str
        Single cell or multicell

    Returns
    -------

    """
    ic_ = update_alpha(ic_)
    oc_l = update_alpha(oc_l)
    if solver.lower() == 'slans':
        if opt:
            folder = 'SLANS_opt'
        else:
            folder = 'SLANS'
    else:
        folder = 'NGSolveMEVP'

    if cell_type is None:
        # print("Writing parameters to file")
        path = fr'{projectDir}/SimulationData/{slans_folder}/{key}/{key}.txt'

        # print(path)
        with open(path, 'w') as f:
            name_list = ['Aeq', 'Beq', 'ai', 'bi', 'Ri', 'L', 'Req', 'alpha', 'Aeq_e', 'Beq_e', 'ai_e', 'bi_e', 'Ri_e',
                         'L_e', 'Req', 'alpha_e', 'key']

            value_list = [ic_[0], ic_[1], ic_[2], ic_[3], ic_[4], ic_[5], ic_[6], ic_[7],
                          oc_l[0], oc_l[1], oc_l[2], oc_l[3], oc_l[4], oc_l[5], oc_l[6], oc_l[7], key]

            for i in range(len(name_list)):
                if name_list[i] == 'key':
                    f.write(f'{name_list[i]} = "{0}" "{value_list[i]}"\n')
                else:
                    f.write(f'{name_list[i]} = "{value_list[i]}" ""\n')

    else:
        # print("Writing parameters to file")
        path = fr'{projectDir}/SimulationData/{folder}/{key}/{key}.txt'
        path_mc = fr'{projectDir}/SimulationData/{folder}/{key}/{key}_Multicell.txt'

        # print(path)
        with open(path, 'w') as f:
            name_list = ['Aeq', 'Beq', 'ai', 'bi', 'Ri', 'L', 'Req', 'alpha', 'Aeq_e', 'Beq_e', 'ai_e', 'bi_e', 'Ri_e',
                         'L_e', 'Req_e', 'alpha_e', 'key']

            if cell_type == 'Mid Cell':
                value_list = [ic_[0], ic_[1], ic_[2], ic_[3], ic_[4], ic_[5], ic_[6], ic_[7],
                              'Aeq', 'Beq', 'ai', 'bi', 'Ri', 'L', 'Req', 'alpha', key]
            else:
                value_list = [ic_[0], ic_[1], ic_[2], ic_[3], ic_[4], ic_[5], ic_[6], ic_[7],
                              oc_r[0], oc_r[1], oc_r[2], oc_r[3], oc_r[4], oc_r[5], ic_[6], oc_l[7],
                              oc_l[0], oc_l[1], oc_l[2], oc_l[3], oc_l[4], oc_l[5], ic_[6], oc_l[7],
                              key]

            for i in range(len(name_list)):
                if name_list[i] == 'key':
                    f.write(f'{name_list[i]} = "{0}" "{value_list[i]}"\n')
                else:
                    f.write(f'{name_list[i]} = "{value_list[i]}" ""\n')

        with open(path_mc, 'w') as f:
            name_list = ['Aeq', 'Beq', 'ai', 'bi', 'Ri', 'L', 'Req', 'alpha',
                         'Aeq_er', 'Beq_er', 'ai_er', 'bi_er', 'Ri_er', 'L_er', 'Req', 'alpha_er',
                         'Aeq_el', 'Beq_el', 'ai_el', 'bi_el', 'Ri_el', 'L_el', 'Req', 'alpha_el', 'key']

            if cell_type == 'Mid Cell':
                value_list = [ic_[0], ic_[1], ic_[2], ic_[3], ic_[4], ic_[5], ic_[6], ic_[7],
                              'Aeq', 'Beq', 'ai', 'bi', 'Ri', 'L', 'Req', 'alpha',
                              'Aeq', 'Beq', 'ai', 'bi', 'Ri', 'L', 'Req', 'alpha',
                              key]
            else:
                value_list = [ic_[0], ic_[1], ic_[2], ic_[3], ic_[4], ic_[5], ic_[6], ic_[7],
                              oc_r[0], oc_r[1], oc_r[2], oc_r[3], oc_r[4], oc_r[5], ic_[6], oc_r[7],
                              oc_l[0], oc_l[1], oc_l[2], oc_l[3], oc_l[4], oc_l[5], ic_[6], oc_l[7],
                              key]

            for i in range(len(name_list)):
                if name_list[i] == 'key':
                    f.write(f'{name_list[i]} = "{0}" "{value_list[i]}"\n')
                else:
                    f.write(f'{name_list[i]} = "{value_list[i]}" ""\n')

        # print("Writing to file complete.")


def stroud(p):
    """
    Stroud-3 method

    Parameters
    ----------
    p: int
        Dimension

    Returns
    -------
    Nodes of quadrature rule in [0,1]**p (column-wise)
    """
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
    """
    Stroud-3 quadrature in :math:`[0,1]^k`

    Parameters
    ----------
    rdim: int
        Dimension of variables
    degree: int
        Degree

    Returns
    -------
    Nodes and corresponding weights
    """
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


def c1_leg_monomial_integral(expon):
    if expon < 0:
        print("\n")
        print("C1_LEG_MONOMIAL_INTEGRAL - Fatal error!")
        print("EXPON < 0.")
        raise ValueError("C1_LEG_MONOMIAL_INTEGRAL - Fatal error!")

    if expon % 2 == 1:
        return 0.0

    value = 2.0 / (expon + 1)
    return value


def cn_leg_03_xiu(n):
    o = 2 * n

    x = np.zeros((n, o))
    w = np.zeros(o)

    expon = 0
    volume = c1_leg_monomial_integral(expon)
    volume = volume ** n

    for j in range(1, o + 1):

        i = 0
        for r in range(1, math.floor(n / 2) + 1):
            arg = (2 * r - 1) * j * np.pi / n
            i += 1
            x[i - 1, j - 1] = np.sqrt(2.0) * np.cos(arg) / np.sqrt(3.0)
            i += 1
            x[i - 1, j - 1] = np.sqrt(2.0) * np.sin(arg) / np.sqrt(3.0)

        if i < n:
            i += 1
            x[i - 1, j - 1] = np.sqrt(2.0) * (-1) ** j / np.sqrt(3.0)
            if n == 1:
                x[i - 1, j - 1] = x[i - 1, j - 1] / np.sqrt(2.0)

    w[0:o] = volume / o

    return x, np.atleast_2d(w).T/np.sum(w)


def cn_leg_03_1(n):
    o = 2 * n

    w = np.zeros(o)
    x = np.zeros((n, o))

    expon = 0
    volume = c1_leg_monomial_integral(expon)
    volume = volume ** n

    for j in range(1, o + 1):

        i = 0

        for r in range(1, math.floor(n / 2) + 1):
            arg = (2 * r - 1) * j * np.pi / n
            i += 1
            x[i - 1, j - 1] = np.sqrt(2.0) * np.cos(arg) / np.sqrt(3.0)
            i += 1
            x[i - 1, j - 1] = np.sqrt(2.0) * np.sin(arg) / np.sqrt(3.0)

        if i < n:
            i += 1
            if n == 1:
                x[i - 1, j - 1] = r8_mop(j) / np.sqrt(3.0)
            else:
                x[i - 1, j - 1] = np.sqrt(2.0) * r8_mop(j) / np.sqrt(3.0)

    w[0:o] = volume / o

    return x, np.atleast_2d(w).T/np.sum(w)


def r8_mop(i):
    if i % 2 == 0:
        value = 1.0
    else:
        value = -1.0

    return value


def cn_leg_05_1(n, option):
    # Check if the value of n is 4, 5, or 6
    if n not in [4, 5, 6]:
        print("\n")
        print("CN_LEG_05_1 - Fatal error!")
        print("The value of N must be 4, 5, or 6.")
        raise ValueError("CN_LEG_05_1 - Fatal error!")

    # Check for valid option when n = 4 or 5
    if n in [4, 5] and option not in [1, 2]:
        print("\n")
        print("CN_LEG_05_1 - Fatal error!")
        print("When N = 4 or 5, OPTION must be 1 or 2.")
        raise ValueError("CN_LEG_05_1 - Fatal error!")

    o = n**2 + n + 2
    w = np.zeros(o)
    x = np.zeros((n, o))

    expon = 0
    volume = c1_leg_monomial_integral(expon)
    volume = volume ** n

    if (n == 4 and option == 1):
        eta = 0.778984505799815
        lmbda = 1.284565137874656
        xsi = -0.713647298819253
        mu = -0.715669761974162
        gamma = 0.217089151000943
        a = 0.206186096875899e-01 * volume
        b = 0.975705820221664e-02 * volume
        c = 0.733921929172573e-01 * volume
    elif ( n == 4 and option == 2 ):
        eta    =   0.546190755827425E+00
        lmbda =   0.745069130115661E+00
        xsi =    - 0.413927294508700E+00
        mu =     - 0.343989637454535E+00
        gamma =    1.134017894600344E+00
        a =        0.853094758323323E-01 * volume
        b =        0.862099000096395E-01 * volume
        c =        0.116418206881849E-01 * volume
    elif ( n == 5 and option == 1 ):
        eta    =   0.522478547481276E+00
        lmbda =   0.936135175985774E+00
        xsi =    - 0.246351362101519E+00
        mu =     - 0.496308106093758E+00
        gamma =    0.827180176822930E+00
        a =        0.631976901960153E-01 * volume
        b =        0.511464127430166E-01 * volume
        c =        0.181070246088902E-01 * volume
    elif ( n == 5 and option == 2 ):
        eta    =   0.798317301388741E+00
        lmbda =   0.637344273885728E+00
        xsi =    - 0.455245909918377E+00
        mu =     - 1.063446229997311E+00
        gamma =    0.354482076665770E+00
        a =        0.116952384292206E-01 * volume
        b =        0.701731258612708E-01 * volume
        c =        0.137439132264426E-01 * volume
    else:
        eta    =   0.660225291773525E+00
        lmbda =   1.064581294844754E+00
        xsi =      0.000000000000000E+00
        mu =     - 0.660225291773525E+00
        gamma =    0.660225291773525E+00
        a =        0.182742214532872E-01 * volume
        b =        0.346020761245675E-01 * volume
        c =        0.182742214532872E-01 * volume

    # Set x and w based on parameters
    k = 0
    # k += 1
    for i in range(n):
        x[i, k] = eta
    w[k] = a

    # k += 1
    for i in range(n):
        x[i, k] = -eta
    w[k] = a

    for i1 in range(n):
        for i in range(1, n):
            x[i,k] = xsi
        x[i1,k] = lmbda
        w[k] = b
        k = k + 1

    for i1 in range(n):
        for i in range(n):
            x[i, k] = - xsi
        x[i1, k] = - lmbda
        w[k] = b
        k = k + 1

    for i1 in range(n-1):
        for i2 in range(i1+1, n):
            for i in range(n):
                x[i, k] = gamma
            x[i1, k] = mu
            x[i2, k] = mu
            w[k] = c
            k = k + 1

    for i1 in range(n-1):
        for i2 in range(i1+1, n):
            for i in range(n):
                x[i, k] = - gamma
            x[i1, k] = - mu
            x[i2, k] = - mu
            w[k] = c
            k = k + 1

    return x, np.atleast_2d(w).T/np.sum(w)


def cn_leg_05_2(n):
    if n < 2:
        print("\n")
        print("CN_LEG_05_2 - Fatal error!")
        print("N must be at least 2.")
        raise ValueError("CN_LEG_05_2 - Fatal error!")

    o = 2 * n**2 + 1
    w = np.zeros(o, dtype=np.float64)
    x = np.zeros((n, o), dtype=np.float64)

    expon = 0
    volume = c1_leg_monomial_integral(expon)
    volume = volume ** n

    b0 = (25 * n * n - 115 * n + 162) * volume / 162.0
    b1 = (70 - 25 * n) * volume / 162.0
    b2 = 25.0 * volume / 324.0

    r = np.sqrt(3.0 / 5.0)

    k = 0

    k += 1
    for i in range(n):
        x[i, k - 1] = 0.0
    w[k - 1] = b0

    for i1 in range(1, n + 1):
        k += 1
        for i in range(n):
            x[i, k - 1] = 0.0
        x[i1 - 1, k - 1] = +r
        w[k - 1] = b1

        k += 1
        for i in range(n):
            x[i, k - 1] = 0.0
        x[i1 - 1, k - 1] = -r
        w[k - 1] = b1

    for i1 in range(1, n):
        for i2 in range(i1 + 1, n + 1):
            k += 1
            for i in range(n):
                x[i, k - 1] = 0.0
            x[i1 - 1, k - 1] = +r
            x[i2 - 1, k - 1] = +r
            w[k - 1] = b2

            k += 1
            for i in range(n):
                x[i, k - 1] = 0.0
            x[i1 - 1, k - 1] = +r
            x[i2 - 1, k - 1] = -r
            w[k - 1] = b2

            k += 1
            for i in range(n):
                x[i, k - 1] = 0.0
            x[i1 - 1, k - 1] = -r
            x[i2 - 1, k - 1] = +r
            w[k - 1] = b2

            k += 1
            for i in range(n):
                x[i, k - 1] = 0.0
            x[i1 - 1, k - 1] = -r
            x[i2 - 1, k - 1] = -r
            w[k - 1] = b2

    return x, np.atleast_2d(w).T/np.sum(w)


def cn_gauss(rdim, degree):
    x, w = np.polynomial.legendre.leggauss(degree)

    X = [x for _ in range(rdim)]

    nodes = np.array(np.meshgrid(*X, indexing='ij')).reshape(rdim, -1)
    weights = np.ones(degree**rdim)/(degree**rdim)

    return nodes, np.atleast_2d(weights).T/np.sum(weights)


def weighted_mean_obj(tab_var, weights):
    rows_sims_no, cols = np.shape(tab_var)
    no_weights, dummy = np.shape(weights)

    if rows_sims_no == no_weights:
        expe = np.zeros((cols, 1))
        outvar = np.zeros((cols, 1))
        for i in range(cols):
            ic(np.dot(tab_var[:, i], weights))
            expe[i, 0] = np.dot(tab_var[:, i], weights)[0]
            outvar[i, 0] = np.dot(tab_var[:, i] ** 2, weights)[0]

        stdDev = np.sqrt(abs(outvar - expe ** 2))
    else:
        expe = 0
        stdDev = 0
        ic('Cols_sims_no != No_weights')

    return list(expe.T[0]), list(stdDev.T[0])


def normal_dist(x, mean, sd):
    prob_density = (np.pi * sd) * np.exp(-0.5 * ((x - mean) / sd) ** 2)
    return prob_density


def text_to_list__(ll):
    """
    Convert text input to list

    Parameters
    ----------
    ll: str
        Input text

    Returns
    -------
    List
    """
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


def animate_width(widget, min_width, standard, enable, option="max"):
    """

    Parameters
    ----------
    widget: QWidget
        QWidget
    min_width: float
        Minimum width
    standard: float
        Default size
    enable:
    option:

    Returns
    -------

    """
    global animate
    if enable:
        # GET WIDTH
        width = widget.width()
        # SET MAX WIDTH
        if width > 0:
            widthCollapsed = min_width
            widget.setMinimumWidth(0)
        else:
            widthCollapsed = standard
            # widget.setMinimumWidth(standard)

        # ANIMATION
        if option == 'max':
            animate = QPropertyAnimation(widget, b"maximumWidth")
        else:
            animate = QPropertyAnimation(widget, b"minimumWidth")

        animate.setDuration(AN_DURATION)
        animate.setStartValue(width)
        animate.setEndValue(widthCollapsed)
        animate.setEasingCurve(QEasingCurve.InOutQuart)
        animate.start()


def animate_height(widget, min_height, standard, enable, option="max"):
    """

    Parameters
    ----------
    widget: QWidget
        QWidget
    min_height: float
        Minimum width
    standard: float
        Default size
    enable:
    option:

    Returns
    -------

    """
    if enable:
        # GET WIDTH
        height = widget.height()

        # SET MAX WIDTH
        if height > 0:
            heightCollapsed = min_height
            widget.setMinimumHeight(0)
        else:
            heightCollapsed = standard
            # ui.w_Shape_Parameters.setMinimumSize(0, 250)

        # ANIMATION
        if option == 'max':
            animation = QPropertyAnimation(widget, b"maximumHeight")
        else:
            animation = QPropertyAnimation(widget, b"minimumHeight")
        animation.setDuration(AN_DURATION)
        animation.setStartValue(height)
        animation.setEndValue(heightCollapsed)
        animation.setEasingCurve(QEasingCurve.InOutQuart)
        animation.start()


def f2b_slashes(path):
    """
    Replaces forward slashes with backward slashes for windows OS

    Parameters
    ----------
    path: str
        Directory path

    Returns
    -------

    """

    if os.name == 'nt':
        path = path.replace("/", "\\")
    else:
        path = path.replace('\\', '/')
    return path


def button_clicked(button):
    """
    Return text on button clicked

    Parameters
    ----------
    button: QPushButton
        PyQt object

    Returns
    -------
    String containing text on PyQt object

    """
    return button.text()


def text_to_list(txt):
    """
    Converts text entered in QLineEdit to Python list
    Parameters
    ----------
    txt: str
        String from QLineEdit

    Returns
    -------

    """
    if "range" in txt:
        txt = txt.replace('range', '')
        ll = ast.literal_eval(txt)
        return range(ll[0], ll[1], ll[2])
    elif 'linspace' in txt:
        ll = eval(f'np.{txt}')
        return ll
    elif txt == '':
        return [1]
    else:
        ll = ast.literal_eval(txt)
        if isinstance(ll, int) or isinstance(ll, float):
            return [ll]
        else:
            return list(ll)


def get_geometric_parameters(frame_control, code, scales=None):
    if scales is None:
        scales = [1]

    shape_space = {}
    for scale in scales:
        loaded_shape_space = copy.deepcopy(
            frame_control.loaded_shape_space)  # have to copy it otherwise the dict is updated
        if frame_control.ui.cb_Shape_Entry_Mode.currentIndex() == 0:
            try:
                # get selected keys
                frame_control.selected_keys = frame_control.ui.cb_Shape_Space_Keys.currentText()
                # ic("selected keys: ", frame_control.ui.cb_Shape_Space_Keys.currentText())
                # print("Selected keys: ", frame_control.selected_keys, type(frame_control.selected_keys[0]))

                # check keys of shape space if results already exist
                to_all = None
                for key, val in loaded_shape_space.items():
                    # process for only keys selected in combobox
                    if frame_control.ui.cb_Shape_Space_Keys.currentText() == "" \
                            or frame_control.ui.cb_Shape_Space_Keys.currentText() == "All":
                        pass
                    else:
                        if key not in frame_control.selected_keys.split(', '):
                            continue

                    if not to_all:
                        ans = frame_control.prompt(key)
                        if ans == 'Yes':
                            if scale == 1 or scale == 0:
                                shape_space[key] = val
                            else:
                                shape_space[f"{key}_s{scale}"] = scale_cavity_geometry(val, scale)

                        if ans == 'No':
                            continue

                        if ans == 'YesToAll':
                            if scale == 1 or scale == 0:
                                shape_space[key] = val
                            else:
                                shape_space[f"{key}_s{scale}"] = scale_cavity_geometry(val, scale)
                            to_all = 'YesToAll'

                        if ans == 'NoToAll':
                            to_all = 'NoToAll'

                        if ans == "Does not exist":
                            if scale == 1 or scale == 0:
                                shape_space[key] = val
                            else:
                                shape_space[f"{key}_s{scale}"] = scale_cavity_geometry(val, scale)
                            to_all = None
                    else:
                        if to_all == 'YesToAll':
                            if scale == 1 or scale == 0:
                                shape_space[key] = val
                            else:
                                shape_space[f"{key}_s{scale}"] = scale_cavity_geometry(val, scale)
                        else:
                            path = f'{frame_control.main_control.projectDir}/SimulationData/{code}/{key}'
                            if os.path.exists(path):
                                continue
                            else:
                                if scale == 1 or scale == 0:
                                    shape_space[key] = val
                                else:
                                    shape_space[f"{key}_s{scale}"] = scale_cavity_geometry(val, scale)

            except Exception as e:
                print(f"File not found, check path:: {e}")
        else:
            if frame_control.ui.cb_Inner_Cell.checkState() == 2:
                # Middle Ellipse data
                A_i_space = text_to_list(frame_control.ui.le_A_i.text())
                B_i_space = text_to_list(frame_control.ui.le_B_i.text())
                a_i_space = text_to_list(frame_control.ui.le_a_i.text())
                b_i_space = text_to_list(frame_control.ui.le_b_i.text())
                Ri_i_space = text_to_list(frame_control.ui.le_Ri_i.text())
                L_i_space = text_to_list(frame_control.ui.le_L_i.text())
                Req_i_space = text_to_list(frame_control.ui.le_Req_i.text())
                alpha_i_space = text_to_list(frame_control.ui.le_Alpha.text())

                inner_cell_space = [A_i_space, B_i_space, a_i_space, b_i_space, Ri_i_space, L_i_space, Req_i_space,
                                    alpha_i_space]
            else:
                inner_cell_space = [[0], [0], [0], [0], [0], [0], [0], [0]]

            if frame_control.ui.cb_Outer_Cell_L.checkState() == 2:
                # Middle Ellipse data
                A_ol_space = text_to_list(frame_control.ui.le_A_ol.text())
                B_ol_space = text_to_list(frame_control.ui.le_B_ol.text())
                a_ol_space = text_to_list(frame_control.ui.le_a_ol.text())
                b_ol_space = text_to_list(frame_control.ui.le_b_ol.text())
                Ri_ol_space = text_to_list(frame_control.ui.le_Ri_ol.text())
                L_ol_space = text_to_list(frame_control.ui.le_L_ol.text())
                Req_ol_space = text_to_list(frame_control.ui.le_Req_ol.text())
                alpha_ol_space = text_to_list(frame_control.ui.le_Alpha_ol.text())

                outer_cell_L_space = [A_ol_space, B_ol_space, a_ol_space, b_ol_space, Ri_ol_space, L_ol_space,
                                      Req_ol_space, alpha_ol_space]
            else:
                outer_cell_L_space = inner_cell_space

            if frame_control.ui.cb_Outer_Cell_R.checkState() == 2:
                # Middle Ellipse data
                A_or_space = text_to_list(frame_control.ui.le_A_or.text())
                B_or_space = text_to_list(frame_control.ui.le_B_or.text())
                a_or_space = text_to_list(frame_control.ui.le_a_or.text())
                b_or_space = text_to_list(frame_control.ui.le_b_or.text())
                Ri_or_space = text_to_list(frame_control.ui.le_Ri_or.text())
                L_or_space = text_to_list(frame_control.ui.le_L_or.text())
                Req_or_space = text_to_list(frame_control.ui.le_Req_or.text())
                alpha_or_space = text_to_list(frame_control.ui.le_Alpha_or.text())

                outer_cell_R_space = [A_or_space, B_or_space, a_or_space, b_or_space, Ri_or_space, L_or_space,
                                      Req_or_space, alpha_or_space]
            else:
                outer_cell_R_space = inner_cell_space
            count = 0

            for A_i in inner_cell_space[0]:
                for B_i in inner_cell_space[1]:
                    for a_i in inner_cell_space[2]:
                        for b_i in inner_cell_space[3]:
                            for Ri_i in inner_cell_space[4]:
                                for L_i in inner_cell_space[5]:
                                    for Req_i in inner_cell_space[6]:
                                        if outer_cell_L_space == inner_cell_space:
                                            inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, 0]
                                            outer_cell_L = inner_cell

                                            if frame_control.ui.cb_LBP.checkState() == 2 and frame_control.ui.cb_RBP.checkState() == 2:
                                                shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L,
                                                                      'OC_R': outer_cell_L, 'BP': 'both', 'FREQ': None}
                                            elif frame_control.ui.cb_LBP.checkState() == 2 and frame_control.ui.cb_RBP.checkState() == 0:
                                                shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L,
                                                                      'OC_R': outer_cell_L, 'BP': 'left', 'FREQ': None}
                                            elif frame_control.ui.cb_LBP.checkState() == 0 and frame_control.ui.cb_RBP.checkState() == 2:
                                                shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L,
                                                                      'OC_R': outer_cell_L, 'BP': 'right', 'FREQ': None}
                                            else:
                                                shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L,
                                                                      'OC_R': outer_cell_L, 'BP': 'none', 'FREQ': None}

                                            count += 1
                                        else:
                                            for A_ol in outer_cell_L_space[0]:
                                                for B_ol in outer_cell_L_space[1]:
                                                    for a_ol in outer_cell_L_space[2]:
                                                        for b_ol in outer_cell_L_space[3]:
                                                            for Ri_ol in outer_cell_L_space[4]:
                                                                for L_ol in outer_cell_L_space[5]:
                                                                    # for Req_ol in outer_cell_L_space[6]:
                                                                    if outer_cell_L_space == outer_cell_R_space:
                                                                        inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i,
                                                                                      Req_i, 0]
                                                                        outer_cell_L = [A_ol, B_ol, a_ol, b_ol, Ri_ol,
                                                                                        L_ol, Req_i, 0]
                                                                        outer_cell_R = outer_cell_L
                                                                        if frame_control.ui.cb_LBP.checkState() == 2 and frame_control.ui.cb_RBP.checkState() == 0:
                                                                            shape_space[count] = {'IC': inner_cell,
                                                                                                  'OC': outer_cell_L,
                                                                                                  'OC_R': outer_cell_R,
                                                                                                  'BP': 'left',
                                                                                                  'FREQ': None}
                                                                        elif frame_control.ui.cb_LBP.checkState() == 0 and frame_control.ui.cb_RBP.checkState() == 2:
                                                                            shape_space[count] = {'IC': inner_cell,
                                                                                                  'OC': outer_cell_L,
                                                                                                  'OC_R': outer_cell_R,
                                                                                                  'BP': 'right',
                                                                                                  'FREQ': None}
                                                                        elif frame_control.ui.cb_LBP.checkState() == 2 and frame_control.ui.cb_RBP.checkState() == 2:
                                                                            shape_space[count] = {'IC': inner_cell,
                                                                                                  'OC': outer_cell_L,
                                                                                                  'OC_R': outer_cell_R,
                                                                                                  'BP': 'both',
                                                                                                  'FREQ': None}
                                                                        else:
                                                                            shape_space[count] = {'IC': inner_cell,
                                                                                                  'OC': outer_cell_L,
                                                                                                  'OC_R': outer_cell_R,
                                                                                                  'BP': 'none',
                                                                                                  'FREQ': None}

                                                                        count += 1
                                                                    else:
                                                                        for A_or in outer_cell_R_space[0]:
                                                                            for B_or in outer_cell_R_space[1]:
                                                                                for a_or in outer_cell_R_space[2]:
                                                                                    for b_or in outer_cell_R_space[3]:
                                                                                        for Ri_or in outer_cell_R_space[
                                                                                            4]:
                                                                                            for L_or in \
                                                                                                    outer_cell_R_space[
                                                                                                        5]:
                                                                                                # for Req_or in outer_cell_R_space[6]:
                                                                                                inner_cell = [A_i, B_i,
                                                                                                              a_i, b_i,
                                                                                                              Ri_i, L_i,
                                                                                                              Req_i, 0]
                                                                                                outer_cell_L = [A_ol,
                                                                                                                B_ol,
                                                                                                                a_ol,
                                                                                                                b_ol,
                                                                                                                Ri_ol,
                                                                                                                L_ol,
                                                                                                                Req_i,
                                                                                                                0]
                                                                                                outer_cell_R = [A_or,
                                                                                                                B_or,
                                                                                                                a_or,
                                                                                                                b_or,
                                                                                                                Ri_or,
                                                                                                                L_or,
                                                                                                                Req_i,
                                                                                                                0]
                                                                                                if frame_control.ui.cb_LBP.checkState() == 2 and frame_control.ui.cb_RBP.checkState() == 0:
                                                                                                    shape_space[
                                                                                                        count] = {
                                                                                                        'IC': inner_cell,
                                                                                                        'OC': outer_cell_L,
                                                                                                        'OC_R': outer_cell_R,
                                                                                                        'BP': 'left',
                                                                                                        'FREQ': None}
                                                                                                elif frame_control.ui.cb_LBP.checkState() == 0 and frame_control.ui.cb_RBP.checkState() == 2:
                                                                                                    shape_space[
                                                                                                        count] = {
                                                                                                        'IC': inner_cell,
                                                                                                        'OC': outer_cell_L,
                                                                                                        'OC_R': outer_cell_R,
                                                                                                        'BP': 'right',
                                                                                                        'FREQ': None}
                                                                                                elif frame_control.ui.cb_LBP.checkState() == 2 and frame_control.ui.cb_RBP.checkState() == 2:
                                                                                                    shape_space[
                                                                                                        count] = {
                                                                                                        'IC': inner_cell,
                                                                                                        'OC': outer_cell_L,
                                                                                                        'OC_R': outer_cell_R,
                                                                                                        'BP': 'both',
                                                                                                        'FREQ': None}
                                                                                                else:
                                                                                                    shape_space[
                                                                                                        count] = {
                                                                                                        'IC': inner_cell,
                                                                                                        'OC': outer_cell_L,
                                                                                                        'OC_R': outer_cell_R,
                                                                                                        'BP': 'none',
                                                                                                        'FREQ': None}

                                                                                                count += 1

            # scale geometry
            if scale != 1 or scale != 0:
                scaled_shape_space = {}
                for key, val in shape_space.items():
                    scaled_shape_space[f"{key}_s{scale}"] = scale_cavity_geometry(val, scale)

                shape_space = scaled_shape_space

            return shape_space

    return shape_space


def scale_cavity_geometry(shape, scale):
    """
    Scale cavity geometry parameters
    Parameters
    ----------
    shape: dict
        Python dictionary of cavity geometry parameters
    scale: float
        Desired scale

    Returns
    -------

    """
    shape['IC'] = list(np.array(shape['IC']) * scale)
    shape['OC'] = list(np.array(shape['OC']) * scale)
    shape['FREQ'] = shape['FREQ'] / scale
    try:
        shape['OC_R'] = list(np.array(shape['OC_R']) * scale)
    except KeyError:
        pass

    return shape


def validating(le, default='1'):
    """
    Validate QLineEdit meant to accept only numbers
    Parameters
    ----------
    le: QLineEdit
        QLineEdit object
    default: str
        Set the value that the QLineEdit.text() defaults to

    Returns
    -------

    """
    validation_rule = QDoubleValidator(0, 1e12, 10)
    if validation_rule.validate(le.text(), 20)[0] == QValidator.Acceptable:
        le.setFocus()
    else:
        le.setText(default)


# def mid_only():
#
#     return inner_cell, inner_cell, inner_cell
#
#
# def mid_left():
#
#     return inner_cell, outer_cell_left, outer_cell_left
#
#
# def mid_left_right():
#
#     return inner_cell, outer_cell_left, outer_cell_right


def open_file(le, start_folder=''):
    """
    Load shape space from json file to Python dictionary, opens a file dialog
    Parameters
    ----------
    le: QLineEdit
        QLineEdit object to get filepath from
    start_folder: str, path
        Opens a file dialog for json file selection

    Returns
    -------

    """

    if os.path.exists(start_folder):
        pass
    else:
        os.mkdir(start_folder)

    filename, _ = QFileDialog.getOpenFileName(None, "Open File", str(start_folder), "Json Files (*.json)")
    try:
        le.setText(filename)

    except Exception as e:
        print('Failed to open file:: ', e)


def load_shape_space(frame_control, le, cb):
    """
    Load shape space from json file to Python dictionary
    Parameters
    ----------
    frame_control: GeometryInputControl
        Frame control object
    le: QLineEdit
        QLineEdit object to get filepath from
    cb: QCheckableComboBox
        QCheckableComboBox to save geometry keys in the json file to for selection later

    Returns
    -------

    """
    filename = le.text()
    if os.path.exists(filename):
        try:
            le.setText(filename)
            with open(filename, 'r') as file:
                dd = json.load(file)

            if dd:
                # clear combobox
                cb.clear()
                cb.addItem('All')

                # populate checkboxes with key
                for col in dd.keys():
                    cb.addItem(fr'{col}')

            frame_control.loaded_shape_space = dd

        except Exception as e:
            print('Failed to open file:: ', e)


def uq(key, shape, qois, n_cells, n_modules, n_modes, f_shift, bc, parentDir, projectDir):
    """
    Uncertainty quantification (work in progress)
    Parameters
    ----------
    key
    shape:

    qois:
        Quantities of interest
    n_cells: int
        Number of cavity cells
    n_modules: int
        Number of cavity modules
    n_modes: int
        Number of modes calculated for
    f_shift:
        Frequency shift (eigenmode solver)
    bc: int
        Boundary condition
    parentDir: str, path
        GUI directory
    projectDir: str, path
        Project directory

    Returns
    -------

    """
    err = False
    result_dict_slans = {}
    slans_obj_list = qois
    for o in qois:
        result_dict_slans[o] = {'expe': [], 'stdDev': []}

    # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
    p_true = shape['IC'][0:5]

    rdim = len(p_true)  # How many variabels will be considered as random in our case 5
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
    p_init = np.zeros(np.shape(p_true))

    no_parm, no_sims = np.shape(nodes)
    delta = 0.005  # or 0.1

    Ttab_val_f = []

    sub_dir = fr'{key}'  # the simulation runs at the quadrature points are saved to the key of mean value run
    for i in range(no_sims):
        skip = False
        p_init[0] = p_true[0] * (1 + delta * nodes[0, i])
        p_init[1] = p_true[1] * (1 + delta * nodes[1, i])
        p_init[2] = p_true[2] * (1 + delta * nodes[2, i])
        p_init[3] = p_true[3] * (1 + delta * nodes[3, i])
        p_init[4] = p_true[4] * (1 + delta * nodes[4, i])

        par_mid = np.append(p_init, shape['IC'][5:]).tolist()
        par_end = par_mid

        # perform checks on geometry
        ok = perform_geometry_checks(par_mid, par_end)
        if not ok:
            err = True
            break
        fid = fr'{key}_Q{i}'

        # check if folder exists and skip if it does
        if os.path.exists(fr'{projectDir}\SimulationData\SLANS\{key}\{fid}'):
            skip = True

        # skip analysis if folder already exists.
        if not skip:
            #  run model using SLANS or CST
            # # create folders for all keys
            slans_geom.createFolder(fid, projectDir, subdir=sub_dir)
            try:
                slans_geom.cavity(n_cells, n_modules, par_mid, par_end, par_end,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
            except KeyError:
                slans_geom.cavity(n_cells, n_modules, par_mid, par_end, par_end,
                                  n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)

        filename = fr'{projectDir}\SimulationData\SLANS\{key}\{fid}\cavity_{bc}.svl'
        if os.path.exists(filename):
            params = fr.svl_reader(filename)
            norm_length = 2 * n_cells * shape['IC'][5]

            qois_result = get_qoi_value(params, slans_obj_list, n_cells, norm_length)
            print(qois_result)
            # sometimes some degenerate shapes are still generated and the solver returns zero
            # for the objective functions, such shapes are considered invalid
            for objr in qois_result:
                if objr == 0:
                    # skip key
                    err = True
                    break

            tab_val_f = qois_result

            Ttab_val_f.append(tab_val_f)
        else:
            err = True

    # # add original point
    # filename = fr'{projectDir}\SimulationData\SLANS\{key}\cavity_33.svl'
    # params = fr.svl_reader(filename)
    # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
    # tab_val_f = obj_result
    # Ttab_val_f.append(tab_val_f)

    # import matplotlib.pyplot as plt
    if not err:
        v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights)

        # append results to dict
        for i, o in enumerate(slans_obj_list):
            result_dict_slans[o]['expe'].append(v_expe_fobj[i])
            result_dict_slans[o]['stdDev'].append(v_stdDev_fobj[i])

            # pdf = normal_dist(np.sort(np.array(Ttab_val_f).T[i]), v_expe_fobj[i], v_stdDev_fobj[i])
            # plt.plot(np.sort(np.array(Ttab_val_f).T[i]), pdf)

        # plt.show()

        with open(fr"{projectDir}\SimulationData\SLANS\{key}\uq.json", 'w') as file:
            file.write(json.dumps(result_dict_slans, indent=4, separators=(',', ': ')))
    else:
        print(fr"There was a problem running UQ analysis for {key}")


def write_cavity_for_custom_eig_solver(file_path, n_cell, mid_cell, end_cell_left=None, end_cell_right=None,
                                       beampipe='none', plot=False):
    """
    Write cavity geometry to be used by Multipac for multipacting analysis
    Parameters
    ----------
    file_path: str
        File path to write geometry to
    n_cell: int
        Number of cavity cells
    mid_cell: list, ndarray
        Array of cavity middle cells' geometric parameters
    end_cell_left: list, ndarray
        Array of cavity left end cell's geometric parameters
    end_cell_right: list, ndarray
        Array of cavity left end cell's geometric parameters
    beampipe: str {"left", "right", "both", "none"}
        Specify if beam pipe is on one or both ends or at no end at all
    plot: bool
        If True, the cavity geometry is plotted for viewing

    Returns
    -------

    """
    if plot:
        plt.rcParams["figure.figsize"] = (12, 2)

    if end_cell_left is None:
        end_cell_left = mid_cell

    if end_cell_right is None:
        if end_cell_left is None:
            end_cell_right = mid_cell
        else:
            end_cell_right = end_cell_left

    A_m, B_m, a_m, b_m, Ri_m, L_m, Req = np.array(mid_cell[:7])*1e-3
    A_el, B_el, a_el, b_el, Ri_el, L_el, Req = np.array(end_cell_left[:7])*1e-3
    A_er, B_er, a_er, b_er, Ri_er, L_er, Req = np.array(end_cell_right[:7])*1e-3

    step = 0.0005  # step in boundary points in mm

    if beampipe.lower() == 'both':
        L_bp_l = 4 * L_m
        L_bp_r = 4 * L_m
    elif beampipe.lower() == 'none':
        L_bp_l = 0.0  # 4 * L_m  #
        L_bp_r = 0.0  # 4 * L_m  #
    elif beampipe.lower() == 'left':
        L_bp_l = 4 * L_m
        L_bp_r = 0.0
    elif beampipe.lower() == 'right':
        L_bp_l = 0.00
        L_bp_r = 4 * L_m
    else:
        L_bp_l = 0.0  # 4 * L_m  #
        L_bp_r = 0.0  # 4 * L_m  #

    # calculate shift
    shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er) / 2

    # calculate angles outside loop
    # CALCULATE x1_el, y1_el, x2_el, y2_el

    df = tangent_coords(A_el, B_el, a_el, b_el, Ri_el, L_el, Req, L_bp_l)
    x1el, y1el, x2el, y2el = df[0]

    # CALCULATE x1, y1, x2, y2
    df = tangent_coords(A_m, B_m, a_m, b_m, Ri_m, L_m, Req, L_bp_l)
    x1, y1, x2, y2 = df[0]

    # CALCULATE x1_er, y1_er, x2_er, y2_er
    df = tangent_coords(A_er, B_er, a_er, b_er, Ri_er, L_er, Req, L_bp_r)
    x1er, y1er, x2er, y2er = df[0]

    with open(file_path, 'w') as fil:
        fil.write("   2.0000000e-03   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")
        fil.write("   1.25000000e-02   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")  # a point inside the structure
        fil.write("  -3.1415927e+00  -2.7182818e+00   0.0000000e+00   0.0000000e+00\n")  # a point outside the structure

        # SHIFT POINT TO START POINT
        start_point = [-shift, 0]
        fil.write(f"  {start_point[1]:.16E}  {start_point[0]:.16E}   3.0000000e+00   0.0000000e+00\n")

        lineTo(start_point, [-shift, Ri_el], step, plot)
        pt = [-shift, Ri_el]
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        # ADD BEAM PIPE LENGTH
        lineTo(pt, [L_bp_l - shift, Ri_el], step, plot)
        pt = [L_bp_l - shift, Ri_el]
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        for n in range(1, n_cell + 1):
            if n == 1:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el], plot)
                pt = [-shift + x1el, y1el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2el, y2el], step, plot)
                pt = [-shift + x2el, y2el]
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_el + L_bp_l - shift, Req - B_el, A_el, B_el, step, pt, [L_bp_l + L_el - shift, Req],
                            plot)
                pt = [L_bp_l + L_el - shift, Req]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                if n_cell == 1:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                                [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, Req], plot)
                    pt = [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step, plot)
                    pt = [L_el + L_er - x1er + + L_bp_l + L_bp_r - shift, y1er]
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                [L_bp_l + L_el + L_er - shift, y1er], plot)

                    pt = [L_bp_l + L_el + L_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_er)
                else:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                                [L_el + L_m - x2 + 2 * L_bp_l - shift, Req], plot)
                    pt = [L_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_el + L_m - x1 + 2 * L_bp_l - shift, y1], step, plot)
                    pt = [L_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                [L_bp_l + L_el + L_m - shift, y1], plot)
                    pt = [L_bp_l + L_el + L_m - shift, Ri_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_m)

            elif n > 1 and n != n_cell:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req], plot)
                pt = [L_bp_l + L_m - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                            [L_m + L_m - x2 + 2 * L_bp_l - shift, Req], plot)
                pt = [L_m + L_m - x2 + 2 * L_bp_l - shift, y2]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                lineTo(pt, [L_m + L_m - x1 + 2 * L_bp_l - shift, y1], step, plot)
                pt = [L_m + L_m - x1 + 2 * L_bp_l - shift, y1]
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                            [L_bp_l + L_m + L_m - shift, y1], plot)
                pt = [L_bp_l + L_m + L_m - shift, Ri_m]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # calculate new shift
                shift = shift - 2 * L_m
            else:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req], plot)
                pt = [L_bp_l + L_m - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                            [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, Req], plot)
                pt = [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                lineTo(pt, [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step, plot)
                pt = [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                            [L_bp_l + L_m + L_er - shift, y1er], plot)
                pt = [L_bp_l + L_m + L_er - shift, Ri_er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        # BEAM PIPE
        # reset shift
        shift = (L_bp_r + L_bp_l + (n_cell - 1) * 2 * L_m + L_el + L_er) / 2
        lineTo(pt, [L_bp_r + L_bp_l + 2 * (n_cell - 1) * L_m + L_el + L_er - shift, Ri_er], step, plot)
        pt = [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, Ri_er]
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   3.0000000e+00   0.0000000e+00\n")

        # END PATH
        lineTo(pt, [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0], step,
               plot)  # to add beam pipe to right
        pt = [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0]
        # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
        # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   0.0000000e+00   0.0000000e+00\n")

        # CLOSE PATH
        lineTo(pt, start_point, step, plot)
        fil.write(f"  {start_point[1]:.16E}  {start_point[0]:.16E}   0.0000000e+00   0.0000000e+00\n")

    if plot:
        plt.tight_layout()
        plt.show()

        plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]


def linspace(start, stop, step=1.):
    """
    Like np.linspace but uses step instead of num
    This is inclusive to stop, so if start=1, stop=3, step=0.5
    Output is: array([1., 1.5, 2., 2.5, 3.])
    """
    if start < stop:
        ll = np.linspace(start, stop, int(np.ceil((stop - start) / abs(step) + 1)))
        if stop not in ll:
            ll = np.append(ll, stop)

        return ll
    else:
        ll = np.linspace(stop, start, int(np.ceil((start - stop) / abs(step) + 1)))
        if start not in ll:
            ll = np.append(ll, start)
        return ll


def lineTo(prevPt, nextPt, step, plot=False):
    if prevPt[0] == nextPt[0]:
        # vertical line
        # check id nextPt is greater
        if prevPt[1] < nextPt[1]:
            py = linspace(prevPt[1], nextPt[1], step)
        else:
            py = linspace(nextPt[1], prevPt[1], step)
            py = py[::-1]
        px = np.ones(len(py)) * prevPt[0]

    elif prevPt[1] == nextPt[1]:
        # horizontal line
        if prevPt[0] < nextPt[1]:
            px = linspace(prevPt[0], nextPt[0], step)
        else:
            px = linspace(nextPt[0], prevPt[0], step)

        py = np.ones(len(px)) * prevPt[1]
    else:
        # calculate angle to get appropriate step size for x and y
        ang = np.arctan((nextPt[1] - prevPt[1]) / (nextPt[0] - prevPt[0]))
        if prevPt[0] < nextPt[0] and prevPt[1] < nextPt[1]:
            px = linspace(prevPt[0], nextPt[0], step * np.cos(ang))
            py = linspace(prevPt[1], nextPt[1], step * np.sin(ang))
        elif prevPt[0] > nextPt[0] and prevPt[1] < nextPt[1]:
            px = linspace(nextPt[0], prevPt[0], step * np.cos(ang))
            px = px[::-1]
            py = linspace(prevPt[1], nextPt[1], step * np.sin(ang))
        elif prevPt[0] < nextPt[0] and prevPt[1] > nextPt[1]:
            px = linspace(prevPt[0], nextPt[0], step * np.cos(ang))
            py = linspace(nextPt[1], prevPt[1], step * np.sin(ang))
            py = py[::-1]
        else:
            px = linspace(nextPt[0], prevPt[0], step * np.cos(ang))
            px = px[::-1]
            py = linspace(nextPt[1], prevPt[1], step * np.sin(ang))
            py = py[::-1]
    if plot:
        plt.plot(px, py) #, marker='x'

    return np.array([px, py]).T


def arcTo(x_center, y_center, a, b, step, start, end, plot=False):
    u = x_center  # <- x-position of the center
    v = y_center  # <- y-position of the center
    a = a  # <- radius on the x-axis
    b = b  # <- radius on the y-axis
    C = np.pi*(a+b)  # <- approximate perimeter of ellipse

    t = np.arange(0, 2 * np.pi, np.pi / int(np.ceil(C/step)))

    x = u + a * np.cos(t)
    y = v + b * np.sin(t)
    pts = np.column_stack((x, y))
    inidx = np.all(np.logical_and(np.array(start) < pts, pts < np.array(end)), axis=1)
    inbox = pts[inidx]
    inbox = inbox[inbox[:, 0].argsort()]

    if plot:
        plt.plot(inbox[:, 0], inbox[:, 1])#, marker='x'

    return inbox


def plot_cavity_geometry(plot, IC, OC, OC_R, BP, n_cell, bc, scale=1):
    """
    Plot cavity geometry for display in w_GeometryView
    Parameters
    ----------
    plot: Plot object
        Matplotlib plot object from plotter
    IC: list, ndarray
        Inner Cell geometric parameters list
    OC: list, ndarray
        Left outer Cell geometric parameters list
    OC_R: list, ndarray
        Right outer Cell geometric parameters list
    BP: str {"left", "right", "both", "none"}
        Specify if beam pipe is on one or both ends or at no end at all
    n_cell: int
        Number of cavity cells
    bc: list, ndarray
        List specifying the left and right boundary conditions' color
    scale: float
        Scale of the cavity geometry

    Returns
    -------

    """
    fig = plot.fig
    ax = plot.ax

    # 21578127116
    A_m, B_m, a_m, b_m, Ri_m, L_m, Req = np.array(IC)[:7] * scale * 1e-3
    A_el, B_el, a_el, b_el, Ri_el, L_el, Req = np.array(OC)[:7] * scale * 1e-3
    A_er, B_er, a_er, b_er, Ri_er, L_er, Req = np.array(OC_R)[:7] * scale * 1e-3

    if BP.lower() == 'both':
        L_bp_l = 4 * L_m
        L_bp_r = 4 * L_m

    elif BP.lower() == 'left':
        L_bp_l = 4 * L_m
        L_bp_r = 0.000

    elif BP.lower() == 'right':
        L_bp_l = 0.000
        L_bp_r = 4 * L_m
    else:
        L_bp_l = 0.000
        L_bp_r = 0.000

    step = 0.0005

    # calculate shift
    shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er) / 2
    # shift = 0
    # shift = L_m  # for end cell

    # calculate angles outside loop
    # CALCULATE x1_el, y1_el, x2_el, y2_el
    data = ([0 + L_bp_l, Ri_el + b_el, L_el + L_bp_l, Req - B_el],
            [a_el, b_el, A_el, B_el])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])

    x1el, y1el, x2el, y2el = fsolve(ellipse_tangent, np.array(
        [a_el + L_bp_l, Ri_el + 0.85 * b_el, L_el - A_el + L_bp_l, Req - 0.85 * B_el]),
                                    args=data,
                                    xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req-0.7*B_m] initial guess

    # CALCULATE x1, y1, x2, y2
    data = ([0 + L_bp_l, Ri_m + b_m, L_m + L_bp_l, Req - B_m],
            [a_m, b_m, A_m, B_m])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    x1, y1, x2, y2 = fsolve(ellipse_tangent,
                            np.array([a_m + L_bp_l, Ri_m + 0.85 * b_m, L_m - A_m + L_bp_l, Req - 0.85 * B_m]),
                            args=data, xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req-0.7*B_m] initial guess

    # CALCULATE x1_er, y1_er, x2_er, y2_er
    data = ([0 + L_bp_r, Ri_er + b_er, L_er + L_bp_r, Req - B_er],
            [a_er, b_er, A_er, B_er])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    x1er, y1er, x2er, y2er = fsolve(ellipse_tangent, np.array(
        [a_er + L_bp_r, Ri_er + 0.85 * b_er, L_er - A_er + L_bp_r, Req - 0.85 * B_er]),
                                    args=data,
                                    xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req-0.7*B_m] initial guess
    geo = []

    # SHIFT POINT TO START POINT
    start_point = [-shift, 0]
    geo.append([start_point[1], start_point[0]])

    lineTo(start_point, [-shift, Ri_el], step)
    pt = [-shift, Ri_el]
    geo.append([pt[1], pt[0]])

    # draw left boundary condition
    ax.plot([-shift, -shift], [-Ri_el, Ri_el],
            [-shift - 0.2 * L_m, -shift - 0.2 * L_m], [-0.5 * Ri_el, 0.5 * Ri_el],
            [-shift - 0.4 * L_m, -shift - 0.4 * L_m], [-0.1 * Ri_el, 0.1 * Ri_el], c=bc[0], lw=4, zorder=100)

    # ADD BEAM PIPE LENGTH
    if L_bp_l != 0:
        lineTo(pt, [L_bp_l - shift, Ri_el], step)
        pt = [L_bp_l - shift, Ri_el]

        geo.append([pt[1], pt[0]])

    for n in range(1, n_cell + 1):
        if n == 1:
            # DRAW ARC:
            pts = arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el])
            pt = [-shift + x1el, y1el]
            for pp in pts:
                geo.append([pp[1], pp[0]])
            geo.append([pt[1], pt[0]])

            # DRAW LINE CONNECTING ARCS
            lineTo(pt, [-shift + x2el, y2el], step)
            pt = [-shift + x2el, y2el]
            geo.append([pt[1], pt[0]])

            # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
            pts = arcTo(L_el + L_bp_l - shift, Req - B_el, A_el, B_el, step, pt, [L_bp_l + L_el - shift, Req])
            pt = [L_bp_l + L_el - shift, Req]
            for pp in pts:
                geo.append([pp[1], pp[0]])
            geo.append([pt[1], pt[0]])

            if n_cell == 1:
                if L_bp_r > 1:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                                [L_el + L_er - x2er + + L_bp_l + L_bp_r - shift, Req])
                    pt = [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            geo.append([pp[1], pp[0]])
                    geo.append([pt[1], pt[0]])

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step)
                    pt = [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                    geo.append([pt[1], pt[0]])

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                [L_bp_l + L_el + L_er - shift, y1er])

                    pt = [L_bp_l + L_el + L_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            geo.append([pp[1], pp[0]])

                    geo.append([pt[1], pt[0]])

                    # calculate new shift
                    shift = shift - (L_el + L_er)
                else:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                                [L_m + L_m - x2 + 2 * L_bp_l - shift, Req])
                    pt = [L_m + L_m - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            geo.append([pp[1], pp[0]])

                    geo.append([pt[1], pt[0]])

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_m + L_m - x1 + 2 * L_bp_l - shift, y1], step)
                    pt = [L_m + L_m - x1 + 2 * L_bp_l - shift, y1]
                    geo.append([pt[1], pt[0]])

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_m + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                [L_bp_l + L_m + L_m - shift, y1])
                    pt = [L_bp_l + L_m + L_m - shift, Ri_m]

                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            geo.append([pp[1], pp[0]])
                    geo.append([pt[1], pt[0]])

                    # calculate new shift
                    shift = shift - L_m

            else:
                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_el + L_bp_l - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                            [L_el + L_m - x2 + 2 * L_bp_l - shift, Req])
                pt = [L_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        geo.append([pp[1], pp[0]])
                geo.append([pt[1], pt[0]])

                # STRAIGHT LINE TO NEXT POINT
                lineTo(pt, [L_el + L_m - x1 + 2 * L_bp_l - shift, y1], step)
                pt = [L_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                geo.append([pt[1], pt[0]])

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                            [L_bp_l + L_el + L_m - shift, y1])
                pt = [L_bp_l + L_el + L_m - shift, Ri_m]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        geo.append([pp[1], pp[0]])
                geo.append([pt[1], pt[0]])

                # calculate new shift
                shift = shift - (L_el + L_m)
                # ic(shift)

        elif n > 1 and n != n_cell:
            # DRAW ARC:
            pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
            pt = [-shift + x1, y1]
            for pp in pts:
                if (np.around(pp, 12) != np.around(pt, 12)).all():
                    geo.append([pp[1], pp[0]])
            geo.append([pt[1], pt[0]])

            # DRAW LINE CONNECTING ARCS
            lineTo(pt, [-shift + x2, y2], step)
            pt = [-shift + x2, y2]
            geo.append([pt[1], pt[0]])

            # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
            pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req])
            pt = [L_bp_l + L_m - shift, Req]
            for pp in pts:
                if (np.around(pp, 12) != np.around(pt, 12)).all():
                    geo.append([pp[1], pp[0]])

            geo.append([pt[1], pt[0]])

            # EQUATOR ARC TO NEXT POINT
            # half of bounding box is required,
            # start is the lower coordinate of the bounding box and end is the upper
            pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                        [L_m + L_m - x2 + 2 * L_bp_l - shift, Req])
            pt = [L_m + L_m - x2 + 2 * L_bp_l - shift, y2]
            for pp in pts:
                if (np.around(pp, 12) != np.around(pt, 12)).all():
                    geo.append([pp[1], pp[0]])

            geo.append([pt[1], pt[0]])

            # STRAIGHT LINE TO NEXT POINT
            lineTo(pt, [L_m + L_m - x1 + 2 * L_bp_l - shift, y1], step)
            pt = [L_m + L_m - x1 + 2 * L_bp_l - shift, y1]
            geo.append([pt[1], pt[0]])

            # ARC
            # half of bounding box is required,
            # start is the lower coordinate of the bounding box and end is the upper
            pts = arcTo(L_m + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                        [L_bp_l + L_m + L_m - shift, y1])
            pt = [L_bp_l + L_m + L_m - shift, Ri_m]

            for pp in pts:
                if (np.around(pp, 12) != np.around(pt, 12)).all():
                    geo.append([pp[1], pp[0]])
            geo.append([pt[1], pt[0]])

            # calculate new shift
            shift = shift - 2 * L_m
        else:
            # DRAW ARC:
            pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
            pt = [-shift + x1, y1]
            for pp in pts:
                if (np.around(pp, 12) != np.around(pt, 12)).all():
                    geo.append([pp[1], pp[0]])
            geo.append([pt[1], pt[0]])

            # DRAW LINE CONNECTING ARCS
            lineTo(pt, [-shift + x2, y2], step)
            pt = [-shift + x2, y2]
            geo.append([pt[1], pt[0]])

            # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
            pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req])
            pt = [L_bp_l + L_m - shift, Req]
            for pp in pts:
                if (np.around(pp, 12) != np.around(pt, 12)).all():
                    geo.append([pp[1], pp[0]])
            geo.append([pt[1], pt[0]])

            # EQUATOR ARC TO NEXT POINT
            # half of bounding box is required,
            # start is the lower coordinate of the bounding box and end is the upper
            pts = arcTo(L_m + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                        [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, Req])
            pt = [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
            for pp in pts:
                if (np.around(pp, 12) != np.around(pt, 12)).all():
                    geo.append([pp[1], pp[0]])
            geo.append([pt[1], pt[0]])

            # STRAIGHT LINE TO NEXT POINT
            lineTo(pt, [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step)
            pt = [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
            geo.append([pt[1], pt[0]])

            # ARC
            # half of bounding box is required,
            # start is the lower coordinate of the bounding box and end is the upper
            pts = arcTo(L_m + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                        [L_bp_l + L_m + L_er - shift, y1er])
            pt = [L_bp_l + L_m + L_er - shift, Ri_er]
            for pp in pts:
                if (np.around(pp, 12) != np.around(pt, 12)).all():
                    geo.append([pp[1], pp[0]])
            geo.append([pt[1], pt[0]])

    # BEAM PIPE
    # reset shift
    shift = (L_bp_r + L_bp_l + (n_cell - 1) * 2 * L_m + L_el + L_er) / 2

    if L_bp_r > 0:  # if there's a problem, check here.
        lineTo(pt, [L_bp_r + L_bp_l + 2 * (n_cell - 1) * L_m + L_el + L_er - shift, Ri_er], step)
        pt = [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, Ri_er]
        geo.append([pt[1], pt[0]])

    # END PATH
    lineTo(pt, [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0], step)  # to add beam pipe to right
    pt = [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0]
    # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
    # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
    geo.append([pt[1], pt[0]])

    # draw rightt boundary condition
    ax.plot([shift, shift], [-Ri_er, Ri_er],
            [shift + 0.2 * L_m, shift + 0.2 * L_m], [-0.5 * Ri_er, 0.5 * Ri_er],
            [shift + 0.4 * L_m, shift + 0.4 * L_m], [-0.1 * Ri_er, 0.1 * Ri_er], c=bc[1], lw=4, zorder=100)

    # CLOSE PATH
    # lineTo(pt, start_point, step)
    # geo.append([start_point[1], start_point[0]])
    geo = np.array(geo)

    top = ax.plot(geo[:, 1], geo[:, 0], lw=4)
    bottom = ax.plot(geo[:, 1], -geo[:, 0], lw=4, c=top[0].get_color())
    fig.canvas.draw()


def orig_writeCavityForMultipac(file_path, n_cell, mid_cell, end_cell_left=None, end_cell_right=None, beampipe='none',
                           plot=False, unit=1e-3, scale=1):
    """
    Write cavity geometry to be used by Multipac for multipacting analysis
    Parameters
    ----------
    file_path: str
        File path to write geometry to
    n_cell: int
        Number of cavity cells
    mid_cell: list, ndarray
        Array of cavity middle cells' geometric parameters
    end_cell_left: list, ndarray
        Array of cavity left end cell's geometric parameters
    end_cell_right: list, ndarray
        Array of cavity left end cell's geometric parameters
    beampipe: str {"left", "right", "both", "none"}
        Specify if beam pipe is on one or both ends or at no end at all
    plot: bool
        If True, the cavity geometry is plotted for viewing

    Returns
    -------

    """
    if plot:
        plt.rcParams["figure.figsize"] = (12, 2)

    if end_cell_left is None:
        end_cell_left = mid_cell

    if end_cell_right is None:
        if end_cell_left is None:
            end_cell_right = mid_cell
        else:
            end_cell_right = end_cell_left

    A_m, B_m, a_m, b_m, Ri_m, L_m, Req = np.array(mid_cell[:7])*unit*scale
    A_el, B_el, a_el, b_el, Ri_el, L_el, Req = np.array(end_cell_left[:7])*unit*scale
    A_er, B_er, a_er, b_er, Ri_er, L_er, Req = np.array(end_cell_right[:7])*unit*scale

    step = 0.001

    if beampipe.lower() == 'both':
        L_bp_l = 4 * L_m
        L_bp_r = 4 * L_m
    elif beampipe.lower() == 'none':
        L_bp_l = 0.000  # 4 * L_m  #
        L_bp_r = 0.000  # 4 * L_m  #
    elif beampipe.lower() == 'left':
        L_bp_l = 4 * L_m
        L_bp_r = 0.000
    elif beampipe.lower() == 'right':
        L_bp_l = 0.000
        L_bp_r = 4 * L_m
    else:
        L_bp_l = 0.000  # 4 * L_m  #
        L_bp_r = 0.000  # 4 * L_m  #

    # calculate shift
    shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er) / 2

    # calculate angles outside loop
    # CALCULATE x1_el, y1_el, x2_el, y2_el

    df = tangent_coords(A_el, B_el, a_el, b_el, Ri_el, L_el, Req, L_bp_l)
    x1el, y1el, x2el, y2el = df[0]

    # CALCULATE x1, y1, x2, y2
    df = tangent_coords(A_m, B_m, a_m, b_m, Ri_m, L_m, Req, L_bp_l)
    x1, y1, x2, y2 = df[0]

    # CALCULATE x1_er, y1_er, x2_er, y2_er
    df = tangent_coords(A_er, B_er, a_er, b_er, Ri_er, L_er, Req, L_bp_r)
    x1er, y1er, x2er, y2er = df[0]

    with open(file_path, 'w') as fil:
        fil.write("   2.0000000e-03   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")
        fil.write("   1.25000000e-02   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")  # a point inside the structure
        fil.write("  -3.1415927e+00  -2.7182818e+00   0.0000000e+00   0.0000000e+00\n")  # a point outside the structure

        # SHIFT POINT TO START POINT
        start_point = [-shift, 0]
        fil.write(f"  {start_point[1]:.16E}  {start_point[0]:.16E}   3.0000000e+00   0.0000000e+00\n")

        pts = lineTo(start_point, [-shift, Ri_el], step, plot)
        pt = [-shift, Ri_el]
        # for pp in pts:
        #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        # ADD BEAM PIPE LENGTH
        pts = lineTo(pt, [L_bp_l - shift, Ri_el], step, plot)
        pt = [L_bp_l - shift, Ri_el]
        for pp in pts:
            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        for n in range(1, n_cell + 1):
            if n == 1:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el], plot)
                pt = [-shift + x1el, y1el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2el, y2el], step, plot)
                pt = [-shift + x2el, y2el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_el + L_bp_l - shift, Req - B_el, A_el, B_el, step, pt, [L_bp_l + L_el - shift, Req],
                            plot)
                pt = [L_bp_l + L_el - shift, Req]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                if n_cell == 1:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                                [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, Req], plot)
                    pt = [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    pts = lineTo(pt, [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step, plot)
                    pt = [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                [L_bp_l + L_el + L_er - shift, y1er], plot)

                    pt = [L_bp_l + L_el + L_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_er)
                else:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                                [L_el + L_m - x2 + 2 * L_bp_l - shift, Req], plot)
                    pt = [L_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    pts = lineTo(pt, [L_el + L_m - x1 + 2 * L_bp_l - shift, y1], step, plot)
                    pt = [L_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                [L_bp_l + L_el + L_m - shift, y1], plot)
                    pt = [L_bp_l + L_el + L_m - shift, Ri_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_m)

            elif n > 1 and n != n_cell:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req], plot)
                pt = [L_bp_l + L_m - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                            [L_m + L_m - x2 + 2 * L_bp_l - shift, Req], plot)
                pt = [L_m + L_m - x2 + 2 * L_bp_l - shift, y2]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                pts = lineTo(pt, [L_m + L_m - x1 + 2 * L_bp_l - shift, y1], step, plot)
                pt = [L_m + L_m - x1 + 2 * L_bp_l - shift, y1]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                            [L_bp_l + L_m + L_m - shift, y1], plot)
                pt = [L_bp_l + L_m + L_m - shift, Ri_m]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # calculate new shift
                shift = shift - 2 * L_m
            else:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req], plot)
                pt = [L_bp_l + L_m - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                            [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, Req], plot)
                pt = [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                pts = lineTo(pt, [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step, plot)
                pt = [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                            [L_bp_l + L_m + L_er - shift, y1er], plot)
                pt = [L_bp_l + L_m + L_er - shift, Ri_er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        # BEAM PIPE
        # reset shift
        shift = (L_bp_r + L_bp_l + (n_cell - 1) * 2 * L_m + L_el + L_er) / 2
        pts = lineTo(pt, [L_bp_r + L_bp_l + 2 * (n_cell - 1) * L_m + L_el + L_er - shift, Ri_er], step, plot)
        pt = [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, Ri_er]
        for pp in pts:
            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   3.0000000e+00   0.0000000e+00\n")

        # END PATH
        pts = lineTo(pt, [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0], step,
               plot)  # to add beam pipe to right
        pt = [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0]
        # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
        # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
        # for pp in pts:
        #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   0.0000000e+00   0.0000000e+00\n")

        # CLOSE PATH
        pts = lineTo(pt, start_point, step, plot)
        # for pp in pts:
        #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {start_point[1]:.16E}  {start_point[0]:.16E}   0.0000000e+00   0.0000000e+00\n")

    if plot:
        plt.tight_layout()
        plt.show()

        plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]


def writeCavityForMultipac(file_path, n_cell, mid_cell, end_cell_left=None, end_cell_right=None, beampipe='none',
                           plot=False, unit=1e-3, scale=1):
    """
    Write cavity geometry to be used by Multipac for multipacting analysis
    Parameters
    ----------
    file_path: str
        File path to write geometry to
    n_cell: int
        Number of cavity cells
    mid_cell: list, ndarray
        Array of cavity middle cells' geometric parameters
    end_cell_left: list, ndarray
        Array of cavity left end cell's geometric parameters
    end_cell_right: list, ndarray
        Array of cavity left end cell's geometric parameters
    beampipe: str {"left", "right", "both", "none"}
        Specify if beam pipe is on one or both ends or at no end at all
    plot: bool
        If True, the cavity geometry is plotted for viewing

    Returns
    -------

    """
    if plot:
        plt.rcParams["figure.figsize"] = (12, 2)

    if end_cell_left is None:
        end_cell_left = mid_cell

    if end_cell_right is None:
        if end_cell_left is None:
            end_cell_right = mid_cell
        else:
            end_cell_right = end_cell_left

    us = unit*scale
    A_m, B_m, a_m, b_m, Ri_m, L_m, Req = np.array(mid_cell[:7])*us
    A_el, B_el, a_el, b_el, Ri_el, L_el, Req = np.array(end_cell_left[:7])*us
    A_er, B_er, a_er, b_er, Ri_er, L_er, Req = np.array(end_cell_right[:7])*us

    step = 0.005

    if beampipe.lower() == 'both':
        L_bp_l = 4 * L_m
        L_bp_r = 4 * L_m
    elif beampipe.lower() == 'none':
        L_bp_l = 0.000  # 4 * L_m  #
        L_bp_r = 0.000  # 4 * L_m  #
    elif beampipe.lower() == 'left':
        L_bp_l = 4 * L_m
        L_bp_r = 0.000
    elif beampipe.lower() == 'right':
        L_bp_l = 0.000
        L_bp_r = 4 * L_m
    else:
        L_bp_l = 0.000  # 4 * L_m  #
        L_bp_r = 0.000  # 4 * L_m  #

    # calculate shift
    shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er) / 2

    # calculate angles outside loop
    # CALCULATE x1_el, y1_el, x2_el, y2_el

    df = tangent_coords(A_el, B_el, a_el, b_el, Ri_el, L_el, Req, L_bp_l)
    x1el, y1el, x2el, y2el = df[0]

    # CALCULATE x1, y1, x2, y2
    df = tangent_coords(A_m, B_m, a_m, b_m, Ri_m, L_m, Req, L_bp_l)
    x1, y1, x2, y2 = df[0]

    # CALCULATE x1_er, y1_er, x2_er, y2_er
    df = tangent_coords(A_er, B_er, a_er, b_er, Ri_er, L_er, Req, L_bp_r)
    x1er, y1er, x2er, y2er = df[0]

    with open(file_path, 'w') as fil:
        fil.write("   2.0000000e-03   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")
        fil.write("   1.25000000e-02   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")  # a point inside the structure
        fil.write("  -3.1415927e+00  -2.7182818e+00   0.0000000e+00   0.0000000e+00\n")  # a point outside the structure

        # SHIFT POINT TO START POINT
        start_point = [-shift, 0]
        fil.write(f"  {start_point[1]:.16E}  {start_point[0]:.16E}   3.0000000e+00   0.0000000e+00\n")

        pts = lineTo(start_point, [-shift, Ri_el], step, plot)
        pt = [-shift, Ri_el]
        # for pp in pts:
        #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        # ADD BEAM PIPE LENGTH
        pts = lineTo(pt, [L_bp_l - shift, Ri_el], step, plot)
        pt = [L_bp_l - shift, Ri_el]
        for pp in pts:
            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        for n in range(1, n_cell + 1):
            if n == 1:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el], plot)
                pt = [-shift + x1el, y1el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2el, y2el], step, plot)
                pt = [-shift + x2el, y2el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_el + L_bp_l - shift, Req - B_el, A_el, B_el, step, pt, [L_bp_l + L_el - shift, Req],
                            plot)
                pt = [L_bp_l + L_el - shift, Req]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                if n_cell == 1:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                                [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, Req], plot)
                    pt = [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    pts = lineTo(pt, [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step, plot)
                    pt = [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                [L_bp_l + L_el + L_er - shift, y1er], plot)

                    pt = [L_bp_l + L_el + L_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")

                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_er)
                else:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                                [L_el + L_m - x2 + 2 * L_bp_l - shift, Req], plot)
                    pt = [L_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    pts = lineTo(pt, [L_el + L_m - x1 + 2 * L_bp_l - shift, y1], step, plot)
                    pt = [L_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                [L_bp_l + L_el + L_m - shift, y1], plot)
                    pt = [L_bp_l + L_el + L_m - shift, Ri_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_m)

            elif n > 1 and n != n_cell:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req], plot)
                pt = [L_bp_l + L_m - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                            [L_m + L_m - x2 + 2 * L_bp_l - shift, Req], plot)
                pt = [L_m + L_m - x2 + 2 * L_bp_l - shift, y2]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                pts = lineTo(pt, [L_m + L_m - x1 + 2 * L_bp_l - shift, y1], step, plot)
                pt = [L_m + L_m - x1 + 2 * L_bp_l - shift, y1]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                            [L_bp_l + L_m + L_m - shift, y1], plot)
                pt = [L_bp_l + L_m + L_m - shift, Ri_m]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # calculate new shift
                shift = shift - 2 * L_m
            else:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req], plot)
                pt = [L_bp_l + L_m - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                            [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, Req], plot)
                pt = [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                pts = lineTo(pt, [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step, plot)
                pt = [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                            [L_bp_l + L_m + L_er - shift, y1er], plot)
                pt = [L_bp_l + L_m + L_er - shift, Ri_er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

        # BEAM PIPE
        # reset shift
        shift = (L_bp_r + L_bp_l + (n_cell - 1) * 2 * L_m + L_el + L_er) / 2
        pts = lineTo(pt, [L_bp_r + L_bp_l + 2 * (n_cell - 1) * L_m + L_el + L_er - shift, Ri_er], step, plot)
        pt = [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, Ri_er]
        for pp in pts:
            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   3.0000000e+00   0.0000000e+00\n")

        # END PATH
        pts = lineTo(pt, [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0], step,
               plot)  # to add beam pipe to right
        pt = [2 * (n_cell - 1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0]
        # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
        # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
        # for pp in pts:
        #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   0.0000000e+00   0.0000000e+00\n")

        # CLOSE PATH
        pts = lineTo(pt, start_point, step, plot)
        # for pp in pts:
        #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {start_point[1]:.16E}  {start_point[0]:.16E}   0.0000000e+00   0.0000000e+00\n")

    if plot:
        plt.tight_layout()
        plt.show()

        plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]


def writeCavityForMultipac_multicell(file_path, n_cell, mid_cell, end_cell_left=None, end_cell_right=None, beampipe='none',
                           plot=False, unit=1e-3, scale=1):
    """
    Write cavity geometry to be used by Multipac for multipacting analysis
    Parameters
    ----------
    file_path: str
        File path to write geometry to
    n_cell: int
        Number of cavity cells
    mid_cell: list, ndarray
        Array of cavity middle cells' geometric parameters
    end_cell_left: list, ndarray
        Array of cavity left end cell's geometric parameters
    end_cell_right: list, ndarray
        Array of cavity left end cell's geometric parameters
    beampipe: str {"left", "right", "both", "none"}
        Specify if beam pipe is on one or both ends or at no end at all
    plot: bool
        If True, the cavity geometry is plotted for viewing

    Returns
    -------

    """
    if plot:
        plt.rcParams["figure.figsize"] = (12, 2)

    if end_cell_left is None:
        end_cell_left = mid_cell

    if end_cell_right is None:
        if end_cell_left is None:
            end_cell_right = mid_cell
        else:
            end_cell_right = end_cell_left

    us = unit*scale
    A_m_array, B_m_array, a_m_array, b_m_array, Ri_m_array, L_m_array, Req_array = np.array(mid_cell[:7])*us
    A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = np.array(end_cell_left[:7])*us
    A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = np.array(end_cell_right[:7])*us
    # ic(np.array(mid_cell[:7]), np.array(mid_cell[:7]).shape)

    step = 0.005

    if beampipe.lower() == 'both':
        L_bp_l = 4 * L_el
        L_bp_r = 4 * L_el
    elif beampipe.lower() == 'none':
        L_bp_l = 0.000  # 4 * L_m  #
        L_bp_r = 0.000  # 4 * L_m  #
    elif beampipe.lower() == 'left':
        L_bp_l = 4 * L_el
        L_bp_r = 0.000
    elif beampipe.lower() == 'right':
        L_bp_l = 0.000
        L_bp_r = 4 * L_el
    else:
        L_bp_l = 0.000  # 4 * L_m  #
        L_bp_r = 0.000  # 4 * L_m  #

    # calculate shift
    if n_cell == 1:
        shift = (L_bp_r + L_bp_l + L_el + L_er) / 2
    else:
        shift = (L_bp_r + L_bp_l + L_el + np.sum(L_m_array) + L_er) / 2

    # calculate angles outside loop
    # CALCULATE x1_el, y1_el, x2_el, y2_el

    # enforce welding by average equator radius between adjacent cells
    if n_cell > 1:
        # ic(Req_array)
        Req_el = (Req_el + Req_array[0][0])/2
    else:
        Req_el = (Req_el + Req_er)/2

    # print(A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, L_bp_l)
    df = tangent_coords(A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, L_bp_l)
    x1el, y1el, x2el, y2el = df[0]

    # # CALCULATE x1, y1, x2, y2
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = A_m_array[0], B_m_array[0], a_m_array[0], b_m_array[0], Ri_m_array[0], L_m_array[0], Req_array[0]
    # df = tangent_coords(A_m, B_m, a_m, b_m, Ri_m, L_m, Req_el, L_bp_l)
    # x1, y1, x2, y2 = df[0]

    # # CALCULATE x1_er, y1_er, x2_er, y2_er
    # df = tangent_coords(A_er, B_er, a_er, b_er, Ri_er, L_er, Req_m, L_bp_r)
    # x1er, y1er, x2er, y2er = df[0]

    with open(file_path, 'w') as fil:
        fil.write("   2.0000000e-03   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")
        fil.write("   1.25000000e-02   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")  # a point inside the structure
        fil.write("  -3.1415927e+00  -2.7182818e+00   0.0000000e+00   0.0000000e+00\n")  # a point outside the structure

        # SHIFT POINT TO START POINT
        start_point = [-shift, 0]
        fil.write(f"  {start_point[1]:.16E}  {start_point[0]:.16E}   3.0000000e+00   0.0000000e+00\n")

        pts = lineTo(start_point, [-shift, Ri_el], step, plot)
        pt = [-shift, Ri_el]
        # for pp in pts:
        #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        # ADD BEAM PIPE LENGTH
        pts = lineTo(pt, [L_bp_l - shift, Ri_el], step, plot)
        pt = [L_bp_l - shift, Ri_el]
        for pp in pts:
            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        for n in range(1, n_cell + 1):
            if n == 1:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el], plot)
                pt = [-shift + x1el, y1el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2el, y2el], step, plot)
                pt = [-shift + x2el, y2el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_el + L_bp_l - shift, Req_el - B_el, A_el, B_el, step, pt, [L_bp_l + L_el - shift, Req_el],
                            plot)
                pt = [L_bp_l + L_el - shift, Req_el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                if n_cell == 1:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper

                    # change parameters to right end cell parameters
                    # CALCULATE x1_er, y1_er, x2_er, y2_er
                    df = tangent_coords(A_er, B_er, a_er, b_er, Ri_er, L_er, Req_el, L_bp_r)
                    x1er, y1er, x2er, y2er = df[0]

                    pts = arcTo(L_el + L_bp_l - shift, Req_el - B_er, A_er, B_er, step, [pt[0], Req_el - B_er],
                                [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, Req_el], plot)
                    pt = [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    pts = lineTo(pt, [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step, plot)
                    pt = [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                [L_bp_l + L_el + L_er - shift, y1er], plot)

                    pt = [L_bp_l + L_el + L_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")

                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_er)
                else:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper

                    # change midcell parameters
                    A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = A_m_array[n-1][0], B_m_array[n-1][0], a_m_array[n-1][0], \
                                                     b_m_array[n-1][0], Ri_m_array[n-1][0], L_m_array[n-1][0], Req_array[n-1][0]

                    # enforce welding by average iris radius between adjacent cells
                    Ri_m = (Ri_m + Ri_m_array[n-1][1]) / 2

                    # CALCULATE x1, y1, x2, y2
                    df = tangent_coords(A_m, B_m, a_m, b_m, Ri_m, L_m, Req_el, L_bp_l)
                    x1, y1, x2, y2 = df[0]

                    pts = arcTo(L_el + L_bp_l - shift, Req_el - B_m, A_m, B_m, step, [pt[0], Req_el - B_m],
                                [L_el + L_m - x2 + 2 * L_bp_l - shift, Req_el], plot)
                    pt = [L_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    pts = lineTo(pt, [L_el + L_m - x1 + 2 * L_bp_l - shift, y1], step, plot)
                    pt = [L_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                [L_bp_l + L_el + L_m - shift, y1], plot)
                    pt = [L_bp_l + L_el + L_m - shift, Ri_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_m_array[n-1][0])

            elif n > 1 and n != n_cell:
                # ic(n)
                # change mid cell parameters
                A_m, B_m, a_m, b_m, _, L_m, Req_m = A_m_array[n-2][1], B_m_array[n-2][1], a_m_array[n-2][1], \
                                                 b_m_array[n-2][1], Ri_m_array[n-2][1], L_m_array[n-2][1], Req_array[n-2][1]

                # enforce welding by average equator radius between adjacent cells
                Req_m = (Req_m + Req_array[n-1][0]) / 2

                # CALCULATE x1, y1, x2, y2
                df = tangent_coords(A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, L_bp_l)
                x1, y1, x2, y2 = df[0]

                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req_m], plot)
                pt = [L_bp_l + L_m - shift, Req_m]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # change midcell parameters
                A_m, B_m, a_m, b_m, Ri_m, L_m, _ = A_m_array[n-1][0], B_m_array[n-1][0], \
                                                       a_m_array[n-1][0], \
                                                       b_m_array[n-1][0], Ri_m_array[n-1][0], \
                                                       L_m_array[n-1][0], Req_array[n-1][0]

                # enforce welding by average iris radius between adjacent cells
                Ri_m = (Ri_m + Ri_m_array[n-1][1]) / 2

                # CALCULATE x1, y1, x2, y2
                df = tangent_coords(A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, L_bp_l)
                x1, y1, x2, y2 = df[0]

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m_array[n-2][1] + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, [pt[0], Req_m - B_m],
                            [L_m_array[n-2][1] + L_m - x2 + 2 * L_bp_l - shift, Req_m], plot)
                pt = [L_m_array[n-2][1] + L_m - x2 + 2 * L_bp_l - shift, y2]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                pts = lineTo(pt, [L_m_array[n-2][1] + L_m - x1 + 2 * L_bp_l - shift, y1], step, plot)
                pt = [L_m_array[n-2][1] + L_m - x1 + 2 * L_bp_l - shift, y1]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m_array[n-2][1] + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                            [L_bp_l + L_m_array[n-2][1] + L_m - shift, y1], plot)
                pt = [L_bp_l + L_m_array[n-2][1] + L_m - shift, Ri_m]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # calculate new shift
                shift = shift - (L_m_array[n-2][1] + L_m_array[n-1][0])
            else:
                # change midcell parameters
                A_m, B_m, a_m, b_m, _, L_m, Req_m = A_m_array[n-2][1], B_m_array[n-2][1], a_m_array[n-2][1], \
                                                 b_m_array[n-2][1], Ri_m_array[n-2][1], L_m_array[n-2][1], Req_array[n-2][1]

                # enforce welding by average equator radius between adjacent cells
                Req_m = (Req_m + Req_er) / 2

                # CALCULATE x1, y1, x2, y2
                df = tangent_coords(A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, L_bp_l)
                x1, y1, x2, y2 = df[0]

                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # change parameters to right end cell parameters
                # CALCULATE x1_er, y1_er, x2_er, y2_er
                df = tangent_coords(A_er, B_er, a_er, b_er, Ri_er, L_er, Req_m, L_bp_r)
                x1er, y1er, x2er, y2er = df[0]

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req_m], plot)
                pt = [L_bp_l + L_m - shift, Req_m]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l - shift, Req_m - B_er, A_er, B_er, step, [pt[0], Req_m - B_er],
                            [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, Req_m], plot)
                pt = [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                pts = lineTo(pt, [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step, plot)
                pt = [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                            [L_bp_l + L_m + L_er - shift, y1er], plot)
                pt = [L_bp_l + L_m + L_er - shift, Ri_er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

        # BEAM PIPE
        # reset shift
        if n_cell == 1:
            shift = (L_bp_r + L_bp_l + L_el + L_er) / 2
            pts = lineTo(pt, [L_bp_r + L_bp_l + L_el + L_er - shift, Ri_er], step, plot)
            pt = [L_el + L_er + L_bp_l + L_bp_r - shift, Ri_er]
            for pp in pts:
                fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
            fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   3.0000000e+00   0.0000000e+00\n")

            # END PATH
            pts = lineTo(pt, [L_el + L_er + L_bp_l + L_bp_r - shift, 0], step,
                   plot)  # to add beam pipe to right
            pt = [L_el + L_er + L_bp_l + L_bp_r - shift, 0]
            # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
            # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
            # for pp in pts:
            #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
            fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   0.0000000e+00   0.0000000e+00\n")

            # CLOSE PATH
            pts = lineTo(pt, start_point, step, plot)
            # for pp in pts:
            #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
            fil.write(f"  {start_point[1]:.16E}  {start_point[0]:.16E}   0.0000000e+00   0.0000000e+00\n")
        else:
            shift = (L_bp_r + L_bp_l + np.sum(L_m_array) + L_el + L_er) / 2
            pts = lineTo(pt, [L_bp_r + L_bp_l + np.sum(L_m_array) + L_el + L_er - shift, Ri_er], step, plot)
            pt = [np.sum(L_m_array) + L_el + L_er + L_bp_l + L_bp_r - shift, Ri_er]
            for pp in pts:
                fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
            fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   3.0000000e+00   0.0000000e+00\n")

            # END PATH
            pts = lineTo(pt, [np.sum(L_m_array) + L_el + L_er + L_bp_l + L_bp_r - shift, 0], step,
                   plot)  # to add beam pipe to right
            pt = [np.sum(L_m_array) + L_el + L_er + L_bp_l + L_bp_r - shift, 0]
            # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
            # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
            # for pp in pts:
            #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
            fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   0.0000000e+00   0.0000000e+00\n")

            # CLOSE PATH
            pts = lineTo(pt, start_point, step, plot)
            # for pp in pts:
            #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
            fil.write(f"  {start_point[1]:.16E}  {start_point[0]:.16E}   0.0000000e+00   0.0000000e+00\n")

    if plot:
        plt.tight_layout()
        plt.show()

        plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]


def orig_writeCavityForMultipac_flat_top(file_path, n_cell, mid_cell, end_cell_left=None, end_cell_right=None, beampipe='none',
                           plot=False, unit=1e-3, scale=1):
    """
    Write cavity geometry to be used by Multipac for multipacting analysis
    Parameters
    ----------
    file_path: str
        File path to write geometry to
    n_cell: int
        Number of cavity cells
    mid_cell: list, ndarray
        Array of cavity middle cells' geometric parameters
    end_cell_left: list, ndarray
        Array of cavity left end cell's geometric parameters
    end_cell_right: list, ndarray
        Array of cavity left end cell's geometric parameters
    beampipe: str {"left", "right", "both", "none"}
        Specify if beam pipe is on one or both ends or at no end at all
    plot: bool
        If True, the cavity geometry is plotted for viewing

    Returns
    -------

    """

    if plot:
        plt.rcParams["figure.figsize"] = (12, 2)

    if end_cell_left is None:
        end_cell_left = mid_cell

    if end_cell_right is None:
        if end_cell_left is None:
            end_cell_right = mid_cell
        else:
            end_cell_right = end_cell_left

    A_m, B_m, a_m, b_m, Ri_m, L_m, Req, lft = np.array(mid_cell[:8])*unit*scale
    A_el, B_el, a_el, b_el, Ri_el, L_el, Req, lft_el = np.array(end_cell_left[:8])*unit*scale
    A_er, B_er, a_er, b_er, Ri_er, L_er, Req, lft_er = np.array(end_cell_right[:8])*unit*scale

    step = 0.001

    if beampipe.lower() == 'both':
        L_bp_l = 4 * L_m
        L_bp_r = 4 * L_m
    elif beampipe.lower() == 'none':
        L_bp_l = 0.000  # 4 * L_m  #
        L_bp_r = 0.000  # 4 * L_m  #
    elif beampipe.lower() == 'left':
        L_bp_l = 4 * L_m
        L_bp_r = 0.000
    elif beampipe.lower() == 'right':
        L_bp_l = 0.000
        L_bp_r = 4 * L_m
    else:
        L_bp_l = 0.000  # 4 * L_m  #
        L_bp_r = 0.000  # 4 * L_m  #

    # calculate shift
    shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er + (n_cell-2)*lft + lft_el + lft_er) / 2

    # calculate angles outside loop
    # CALCULATE x1_el, y1_el, x2_el, y2_el

    df = tangent_coords(A_el, B_el, a_el, b_el, Ri_el, L_el, Req, L_bp_l)
    x1el, y1el, x2el, y2el = df[0]

    # CALCULATE x1, y1, x2, y2
    df = tangent_coords(A_m, B_m, a_m, b_m, Ri_m, L_m, Req, L_bp_l)
    x1, y1, x2, y2 = df[0]

    # CALCULATE x1_er, y1_er, x2_er, y2_er
    df = tangent_coords(A_er, B_er, a_er, b_er, Ri_er, L_er, Req, L_bp_r)
    x1er, y1er, x2er, y2er = df[0]

    with open(file_path, 'w') as fil:
        fil.write("   2.0000000e-03   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")
        fil.write("   1.25000000e-02   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")  # a point inside the structure
        fil.write("  -3.1415927e+00  -2.7182818e+00   0.0000000e+00   0.0000000e+00\n")  # a point outside the structure

        # SHIFT POINT TO START POINT
        start_point = [-shift, 0]
        fil.write(f"  {start_point[1]:.16E}  {start_point[0]:.16E}   3.0000000e+00   0.0000000e+00\n")

        pts = lineTo(start_point, [-shift, Ri_el], step, plot)
        pt = [-shift, Ri_el]
        # for pp in pts:
        #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        # ADD BEAM PIPE LENGTH
        if L_bp_l != 0:
            pts = lineTo(pt, [L_bp_l - shift, Ri_el], step, plot)
            pt = [L_bp_l - shift, Ri_el]
            for pp in pts:
                fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
            fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        for n in range(1, n_cell + 1):
            if n == 1:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el], plot)
                pt = [-shift + x1el, y1el]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2el, y2el], step, plot)
                pt = [-shift + x2el, y2el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_el + L_bp_l - shift, Req - B_el, A_el, B_el, step, pt, [L_bp_l + L_el - shift, Req], plot)
                pt = [L_bp_l + L_el - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # flat top
                pts = lineTo(pt, [L_bp_l + L_el + lft_el - shift, Req], step, plot)
                pt = [L_bp_l + L_el + lft_el - shift, Req]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                if n_cell == 1:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l + lft_el - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                                [L_el + lft_el + L_er - x2er + 2 * L_bp_l - shift, Req], plot)
                    pt = [L_el + lft_el + L_er - x2er + 2 * L_bp_l - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    pts = lineTo(pt, [L_el + lft_el + L_er - x1er + 2 * L_bp_l - shift, y1er], step, plot)
                    pt = [L_el + lft_el + L_er - x1er + 2 * L_bp_l - shift, y1er]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + lft_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                [L_bp_l + L_el + lft_el + L_er - shift, y1er], plot)

                    pt = [L_bp_l + lft_el + L_el + L_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_er + lft_el)
                else:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l + lft_el - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                                [L_el + lft_el + L_m - x2 + 2 * L_bp_l - shift, Req], plot)
                    pt = [L_el + lft_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    pts = lineTo(pt, [L_el + lft_el + L_m - x1 + 2 * L_bp_l - shift, y1], step, plot)
                    pt = [L_el + lft_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + lft_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                [L_bp_l + L_el + lft_el + L_m - shift, y1], plot)
                    pt = [L_bp_l + L_el + lft_el + L_m - shift, Ri_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_m + lft_el)
                    # ic(shift)

            elif n > 1 and n != n_cell:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req], plot)
                pt = [L_bp_l + L_m - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # flat top
                pts = lineTo(pt, [L_bp_l + L_m + lft - shift, Req], step, plot)
                pt = [L_bp_l + L_el + lft - shift, Req]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l + lft - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                            [L_m + L_m + lft - x2 + 2 * L_bp_l - shift, Req], plot)
                pt = [L_m + L_m + lft - x2 + 2 * L_bp_l - shift, y2]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                pts = lineTo(pt, [L_m + L_m + lft - x1 + 2 * L_bp_l - shift, y1], step, plot)
                pt = [L_m + L_m + lft - x1 + 2 * L_bp_l - shift, y1]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_m + lft + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                            [L_bp_l + L_m + L_m + lft - shift, y1], plot)
                pt = [L_bp_l + L_m + L_m + lft - shift, Ri_m]
                ic(pt)
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # calculate new shift
                shift = shift - 2*L_m - (lft_el + lft)
            else:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req], plot)
                pt = [L_bp_l + L_m - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # flat top
                pts = lineTo(pt, [L_bp_l + L_m + lft_er - shift, Req], step, plot)
                pt = [L_bp_l + L_m + lft_er - shift, Req, Req]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + lft_er + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                            [L_m + L_er + lft_er - x2er + L_bp_l + L_bp_r - shift, Req], plot)
                pt = [L_m + L_er + lft_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                pts = lineTo(pt, [L_m + L_er + lft_er - x1er + L_bp_l + L_bp_r - shift, y1er], step, plot)
                pt = [L_m + L_er + lft_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_er + lft_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                            [L_bp_l + L_m + L_er + lft_er - shift, y1er], plot)
                pt = [L_bp_l + L_m + L_er + lft_er - shift, Ri_er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        # BEAM PIPE
        # reset shift

        shift = (L_bp_r + L_bp_l + L_el + lft_el + (n_cell - 1) * 2 * L_m + (n_cell - 2)*lft + L_er + lft_er) / 2
        pts = lineTo(pt, [L_bp_r + L_bp_l + 2 * (n_cell-1) * L_m + (n_cell-2)*lft + lft_el + lft_er + L_el + L_er - shift, Ri_er], step, plot)

        if L_bp_r != 0:
            pt = [2 * (n_cell-1) * L_m + L_el + L_er + L_bp_l + L_bp_r + (n_cell-2)*lft + lft_el + lft_er - shift, Ri_er]
            for pp in pts:
                fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
            fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        # END PATH
        pts = lineTo(pt, [2 * (n_cell-1) * L_m + L_el + L_er + (n_cell-2)*lft + lft_el + lft_er + L_bp_l + L_bp_r - shift, 0], step, plot)  # to add beam pipe to right
        pt = [2 * (n_cell-1) * L_m + L_el + L_er + (n_cell-2)*lft + lft_el + lft_er + L_bp_l + L_bp_r - shift, 0]
        # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
        # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
        # for pp in pts:
        #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   0.0000000e+00   0.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   0.0000000e+00   1.0000000e+00\n")

        # CLOSE PATH
        pts = lineTo(pt, start_point, step, plot)
        # for pp in pts:
        #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   0.0000000e+00   0.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   0.0000000e+00   0.0000000e+00\n")

    print(' plot state', plot)
    if plot:
        print(' in here plotting')
        plt.show()


def writeCavityForMultipac_flat_top(file_path, n_cell, mid_cell, end_cell_left=None, end_cell_right=None, beampipe='none',
                           plot=False, unit=1e-3, scale=1):
    """
    Write cavity geometry to be used by Multipac for multipacting analysis
    Parameters
    ----------
    file_path: str
        File path to write geometry to
    n_cell: int
        Number of cavity cells
    mid_cell: list, ndarray
        Array of cavity middle cells' geometric parameters
    end_cell_left: list, ndarray
        Array of cavity left end cell's geometric parameters
    end_cell_right: list, ndarray
        Array of cavity left end cell's geometric parameters
    beampipe: str {"left", "right", "both", "none"}
        Specify if beam pipe is on one or both ends or at no end at all
    plot: bool
        If True, the cavity geometry is plotted for viewing

    Returns
    -------

    """

    if plot:
        plt.rcParams["figure.figsize"] = (12, 2)

    if end_cell_left is None:
        end_cell_left = mid_cell

    if end_cell_right is None:
        if end_cell_left is None:
            end_cell_right = mid_cell
        else:
            end_cell_right = end_cell_left

    A_m, B_m, a_m, b_m, Ri_m, L_m, Req, lft = np.array(mid_cell[:8])*unit*scale
    A_el, B_el, a_el, b_el, Ri_el, L_el, Req, lft_el = np.array(end_cell_left[:8])*unit*scale
    A_er, B_er, a_er, b_er, Ri_er, L_er, Req, lft_er = np.array(end_cell_right[:8])*unit*scale

    step = 0.001

    if beampipe.lower() == 'both':
        L_bp_l = 4 * L_m
        L_bp_r = 4 * L_m
    elif beampipe.lower() == 'none':
        L_bp_l = 0.000  # 4 * L_m  #
        L_bp_r = 0.000  # 4 * L_m  #
    elif beampipe.lower() == 'left':
        L_bp_l = 4 * L_m
        L_bp_r = 0.000
    elif beampipe.lower() == 'right':
        L_bp_l = 0.000
        L_bp_r = 4 * L_m
    else:
        L_bp_l = 0.000  # 4 * L_m  #
        L_bp_r = 0.000  # 4 * L_m  #

    # calculate shift
    shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er + (n_cell-2)*lft + lft_el + lft_er) / 2

    # calculate angles outside loop
    # CALCULATE x1_el, y1_el, x2_el, y2_el

    df = tangent_coords(A_el, B_el, a_el, b_el, Ri_el, L_el, Req, L_bp_l)
    x1el, y1el, x2el, y2el = df[0]

    # CALCULATE x1, y1, x2, y2
    df = tangent_coords(A_m, B_m, a_m, b_m, Ri_m, L_m, Req, L_bp_l)
    x1, y1, x2, y2 = df[0]

    # CALCULATE x1_er, y1_er, x2_er, y2_er
    df = tangent_coords(A_er, B_er, a_er, b_er, Ri_er, L_er, Req, L_bp_r)
    x1er, y1er, x2er, y2er = df[0]

    with open(file_path, 'w') as fil:
        fil.write("   2.0000000e-03   0.0000000e+00   0.0000001e+00   0.0000000e+00\n")
        fil.write("   1.25000000e-02   0.0000000e+00   0.0000001e+00   0.0000000e+00\n")  # a point inside the structure
        fil.write("  -3.1415927e+00  -2.7182818e+00   0.0000001e+00   0.0000000e+00\n")  # a point outside the structure

        # SHIFT POINT TO START POINT
        start_point = [-shift, 0]
        fil.write(f"  {start_point[1]:.16E}  {start_point[0]:.16E}   3.0000001e+00   0.0000000e+00\n")

        pts = lineTo(start_point, [-shift, Ri_el], step, plot)
        pt = [-shift, Ri_el]
        # for pp in pts:
        #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000001e+00   1.0000000e+00\n")

        # ADD BEAM PIPE LENGTH
        if L_bp_l != 0:
            pts = lineTo(pt, [L_bp_l - shift, Ri_el], step, plot)
            pt = [L_bp_l - shift, Ri_el]
            for pp in pts:
                fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000001e+00   1.0000000e+00\n")
            fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000001e+00   1.0000000e+00\n")

        for n in range(1, n_cell + 1):
            if n == 1:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el], plot)
                pt = [-shift + x1el, y1el]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2el, y2el], step, plot)
                pt = [-shift + x2el, y2el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_el + L_bp_l - shift, Req - B_el, A_el, B_el, step, pt, [L_bp_l + L_el - shift, Req], plot)
                pt = [L_bp_l + L_el - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # flat top
                pts = lineTo(pt, [L_bp_l + L_el + lft_el - shift, Req], step, plot)
                pt = [L_bp_l + L_el + lft_el - shift, Req]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                if n_cell == 1:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l + lft_el - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                                [L_el + lft_el + L_er - x2er + 2 * L_bp_l - shift, Req], plot)
                    pt = [L_el + lft_el + L_er - x2er + 2 * L_bp_l - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}0   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    pts = lineTo(pt, [L_el + lft_el + L_er - x1er + 2 * L_bp_l - shift, y1er], step, plot)
                    pt = [L_el + lft_el + L_er - x1er + 2 * L_bp_l - shift, y1er]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + lft_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                [L_bp_l + L_el + lft_el + L_er - shift, y1er], plot)

                    pt = [L_bp_l + lft_el + L_el + L_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")

                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_er + lft_el)
                else:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l + lft_el - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                                [L_el + lft_el + L_m - x2 + 2 * L_bp_l - shift, Req], plot)
                    pt = [L_el + lft_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    pts = lineTo(pt, [L_el + lft_el + L_m - x1 + 2 * L_bp_l - shift, y1], step, plot)
                    pt = [L_el + lft_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                    for pp in pts:
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + lft_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                [L_bp_l + L_el + lft_el + L_m - shift, y1], plot)
                    pt = [L_bp_l + L_el + lft_el + L_m - shift, Ri_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_m + lft_el)
                    # ic(shift)

            elif n > 1 and n != n_cell:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req], plot)
                pt = [L_bp_l + L_m - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # flat top
                pts = lineTo(pt, [L_bp_l + L_m + lft - shift, Req], step, plot)
                pt = [L_bp_l + L_el + lft - shift, Req]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l + lft - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                            [L_m + L_m + lft - x2 + 2 * L_bp_l - shift, Req], plot)
                pt = [L_m + L_m + lft - x2 + 2 * L_bp_l - shift, y2]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                pts = lineTo(pt, [L_m + L_m + lft - x1 + 2 * L_bp_l - shift, y1], step, plot)
                pt = [L_m + L_m + lft - x1 + 2 * L_bp_l - shift, y1]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_m + lft + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                            [L_bp_l + L_m + L_m + lft - shift, y1], plot)
                pt = [L_bp_l + L_m + L_m + lft - shift, Ri_m]
                ic(pt)
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # calculate new shift
                shift = shift - 2*L_m - (lft_el + lft)
            else:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                pts = lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req], plot)
                pt = [L_bp_l + L_m - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # flat top
                pts = lineTo(pt, [L_bp_l + L_m + lft_er - shift, Req], step, plot)
                pt = [L_bp_l + L_m + lft_er - shift, Req, Req]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + lft_er + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                            [L_m + L_er + lft_er - x2er + L_bp_l + L_bp_r - shift, Req], plot)
                pt = [L_m + L_er + lft_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                pts = lineTo(pt, [L_m + L_er + lft_er - x1er + L_bp_l + L_bp_r - shift, y1er], step, plot)
                pt = [L_m + L_er + lft_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_er + lft_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                            [L_bp_l + L_m + L_er + lft_er - shift, y1er], plot)
                pt = [L_bp_l + L_m + L_er + lft_er - shift, Ri_er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   {n}   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   {n}   1.0000000e+00\n")

        # BEAM PIPE
        # reset shift

        shift = (L_bp_r + L_bp_l + L_el + lft_el + (n_cell - 1) * 2 * L_m + (n_cell - 2)*lft + L_er + lft_er) / 2
        pts = lineTo(pt, [L_bp_r + L_bp_l + 2 * (n_cell-1) * L_m + (n_cell-2)*lft + lft_el + lft_er + L_el + L_er - shift, Ri_er], step, plot)

        if L_bp_r != 0:
            pt = [2 * (n_cell-1) * L_m + L_el + L_er + L_bp_l + L_bp_r + (n_cell-2)*lft + lft_el + lft_er - shift, Ri_er]
            for pp in pts:
                fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000001e+00   1.0000000e+00\n")
            fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000001e+00   1.0000000e+00\n")

        # END PATH
        pts = lineTo(pt, [2 * (n_cell-1) * L_m + L_el + L_er + (n_cell-2)*lft + lft_el + lft_er + L_bp_l + L_bp_r - shift, 0], step, plot)  # to add beam pipe to right
        pt = [2 * (n_cell-1) * L_m + L_el + L_er + (n_cell-2)*lft + lft_el + lft_er + L_bp_l + L_bp_r - shift, 0]
        # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
        # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
        # for pp in pts:
        #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   0.0000000e+00   0.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   0.0000001e+00   1.0000000e+00\n")

        # CLOSE PATH
        pts = lineTo(pt, start_point, step, plot)
        # for pp in pts:
        #     fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   0.0000000e+00   0.0000000e+00\n")
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   0.0000001e+00   0.0000000e+00\n")

    print(' plot state', plot)
    if plot:
        print(' in here plotting')
        plt.show()


def write_geometry_ngsolve(file_path, n_cell, mid_cell, end_cell_left=None, end_cell_right=None, beampipe='none',
                           plot=False):
    """
    Write cavity geometry to be used by Multipac for multipacting analysis
    Parameters
    ----------
    file_path: str
        File path to write geometry to
    n_cell: int
        Number of cavity cells
    mid_cell: list, ndarray
        Array of cavity middle cells' geometric parameters
    end_cell_left: list, ndarray
        Array of cavity left end cell's geometric parameters
    end_cell_right: list, ndarray
        Array of cavity left end cell's geometric parameters
    beampipe: str {"left", "right", "both", "none"}
        Specify if beam pipe is on one or both ends or at no end at all
    plot: bool
        If True, the cavity geometry is plotted for viewing

    Returns
    -------

    """
    if plot:
        plt.rcParams["figure.figsize"] = (12, 2)

    if end_cell_left is None:
        end_cell_left = mid_cell

    if end_cell_right is None:
        if end_cell_left is None:
            end_cell_right = mid_cell
        else:
            end_cell_right = end_cell_left

    # # TESLA end cell 2
    A_m, B_m, a_m, b_m, Ri_m, L_m, Req = np.array(mid_cell[:7])*1e-3
    A_el, B_el, a_el, b_el, Ri_el, L_el, Req = np.array(end_cell_left[:7])*1e-3
    A_er, B_er, a_er, b_er, Ri_er, L_er, Req = np.array(end_cell_right[:7])*1e-3

    step = 0.0005

    if beampipe.lower() == 'both':
        L_bp_l = 4 * L_m
        L_bp_r = 4 * L_m
    elif beampipe.lower() == 'none':
        L_bp_l = 0.0001  # 4 * L_m  #
        L_bp_r = 0.0001  # 4 * L_m  #
    elif beampipe.lower() == 'left':
        L_bp_l = 4 * L_m
        L_bp_r = 0.0001
    elif beampipe.lower() == 'right':
        L_bp_l = 0.0001
        L_bp_r = 4 * L_m
    else:
        L_bp_l = 0.0001  # 4 * L_m  #
        L_bp_r = 0.0001  # 4 * L_m  #

    # calculate shift
    shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er) / 2
    # shift = 0
    # shift = L_m  # for end cell

    df = tangent_coords(A_el, B_el, a_el, b_el, Ri_el, L_el, Req, L_bp_l)
    x1el, y1el, x2el, y2el = df[0]

    # CALCULATE x1, y1, x2, y2
    df = tangent_coords(A_m, B_m, a_m, b_m, Ri_m, L_m, Req, L_bp_l)
    x1, y1, x2, y2 = df[0]

    # CALCULATE x1_er, y1_er, x2_er, y2_er
    df = tangent_coords(A_er, B_er, a_er, b_er, Ri_er, L_er, Req, L_bp_r)
    x1er, y1er, x2er, y2er = df[0]
    default_folder = "D:\Dropbox\multipacting\MPGUI21"

    if file_path is None:
        file_path = default_folder

    with open(fr'{file_path}\geodata.n', 'w') as fil:
        fil.write("   2.0000000e-03   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")
        fil.write("   1.25000000e-02   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")  # a point inside the structure
        fil.write("  -3.1415927e+00  -2.7182818e+00   0.0000000e+00   0.0000000e+00\n")  # a point outside the structure

        # SHIFT POINT TO START POINT
        start_point = [-shift, 0]
        fil.write(f"  {start_point[1]:.16E}  {start_point[0]:.16E}   3.0000000e+00   0.0000000e+00\n")

        lineTo(start_point, [-shift, Ri_el], step)
        pt = [-shift, Ri_el]
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        # ADD BEAM PIPE LENGTH
        if L_bp_l != 0:
            lineTo(pt, [L_bp_l - shift, Ri_el], step)
            pt = [L_bp_l - shift, Ri_el]
            print(pt)
            fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        for n in range(1, n_cell + 1):
            if n == 1:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el])
                pt = [-shift + x1el, y1el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2el, y2el], step)
                pt = [-shift + x2el, y2el]
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_el + L_bp_l - shift, Req - B_el, A_el, B_el, step, pt, [L_bp_l + L_el - shift, Req])
                pt = [L_bp_l + L_el - shift, Req]
                for pp in pts:
                    fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                if n_cell == 1:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                                [L_el + L_er - x2er + 2 * L_bp_l - shift, Req])
                    pt = [L_el + L_er - x2er + 2 * L_bp_l - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_el + L_er - x1er + 2 * L_bp_l - shift, y1er], step)
                    pt = [L_el + L_er - x1er + 2 * L_bp_l - shift, y1er]
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                [L_bp_l + L_el + L_er - shift, y1er])

                    pt = [L_bp_l + L_el + L_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_er)
                    # ic(shift)
                else:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                                [L_el + L_m - x2 + 2 * L_bp_l - shift, Req])
                    pt = [L_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_el + L_m - x1 + 2 * L_bp_l - shift, y1], step)
                    pt = [L_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                [L_bp_l + L_el + L_m - shift, y1])
                    pt = [L_bp_l + L_el + L_m - shift, Ri_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                    fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_m)
                    # ic(shift)

            elif n > 1 and n != n_cell:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2, y2], step)
                pt = [-shift + x2, y2]
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req])
                pt = [L_bp_l + L_m - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, [pt[0], Req - B_m],
                            [L_m + L_m - x2 + 2 * L_bp_l - shift, Req])
                pt = [L_m + L_m - x2 + 2 * L_bp_l - shift, y2]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                lineTo(pt, [L_m + L_m - x1 + 2 * L_bp_l - shift, y1], step)
                pt = [L_m + L_m - x1 + 2 * L_bp_l - shift, y1]
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                            [L_bp_l + L_m + L_m - shift, y1])
                pt = [L_bp_l + L_m + L_m - shift, Ri_m]
                ic(pt)
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # calculate new shift
                shift = shift - 2*L_m
            else:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2, y2], step)
                pt = [-shift + x2, y2]
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req])
                pt = [L_bp_l + L_m - shift, Req]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l - shift, Req - B_er, A_er, B_er, step, [pt[0], Req - B_er],
                            [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, Req])
                pt = [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                lineTo(pt, [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step)
                pt = [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                            [L_bp_l + L_m + L_er - shift, y1er])
                pt = [L_bp_l + L_m + L_er - shift, Ri_er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.16E}  {pp[0]:.16E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   1.0000000e+00   1.0000000e+00\n")

        # BEAM PIPE
        # reset shift
        shift = (L_bp_r + L_bp_l + (n_cell - 1) * 2 * L_m + L_el + L_er) / 2

        if L_bp_r != 0:  # if there's a problem, check here.
            lineTo(pt, [L_bp_r + L_bp_l + 2 * (n_cell-1) * L_m + L_el + L_er - shift, Ri_er], step)
            pt = [2 * (n_cell-1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, Ri_er]
            fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   3.0000000e+00   0.0000000e+00\n")

        # END PATH
        lineTo(pt, [2 * (n_cell-1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0], step)  # to add beam pipe to right
        pt = [2 * (n_cell-1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0]
        # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
        # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
        fil.write(f"  {pt[1]:.16E}  {pt[0]:.16E}   0.0000000e+00   0.0000000e+00\n")

        # CLOSE PATH
        lineTo(pt, start_point, step)
        fil.write(f"  {start_point[1]:.16E}  {start_point[0]:.16E}   0.0000000e+00   0.0000000e+00\n")
    plt.gca().set_aspect('equal', 'box')
    plt.show()


def combine_params_output(folder, N):
    for i in range(N):
        if i == 0:
            df = pd.read_csv(f'{folder}/table_{i}.csv', sep='\t', engine='python')
        else:
            try:
                df = pd.concat([df, pd.read_csv(f'{folder}/table_{i}.csv', sep='\t', engine='python')])
            except:
                pass

    df.to_csv(fr'{folder}\table.csv', index=False, sep='\t', float_format='%.32f')
    df.to_excel(fr"{folder}\table.xlsx", index=False)


def create_multicell_random_variables(n_cell, mid_cell, end_cell_l, end_cell_r):
    mid_cell_pair = np.hstack([mid_cell, mid_cell])
    mid_muilticell = np.hstack([mid_cell_pair for _ in range(n_cell - 2)]).reshape(len(mid_cell), n_cell-2, 2)
    return mid_muilticell


def uq_from_excel(filepath):
    Ttab_val_f = pd.read_excel(fr'{filepath}\table.xlsx').to_numpy()
    # nodes, weights_  = cn_leg_05_2(28)
    nodes, weights_ = cn_leg_03_1(28)
    # nodes, weights_, _ = quad_stroud3(5, 2)
    ic(np.sum(weights_), nodes.shape)
    ic(np.max(Ttab_val_f, axis=0))
    ic(np.argmax(Ttab_val_f, axis=0))

    eigen_obj_list = ["freq [MHz]", "R/Q [Ohm]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "G [Ohm]", "kcc [%]", "ff [%]"]
    result_dict_eigen = {}
    for o in eigen_obj_list:
        result_dict_eigen[o] = {'expe': [], 'stdDev': []}

    ic(Ttab_val_f, weights_)
    v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(Ttab_val_f, weights_)

    # append results to dict
    for i, o in enumerate(eigen_obj_list):
        result_dict_eigen[o]['expe'].append(v_expe_fobj[i])
        result_dict_eigen[o]['stdDev'].append(v_stdDev_fobj[i])

        # pdf = normal_dist(np.sort(np.array(Ttab_val_f).T[i]), v_expe_fobj[i], v_stdDev_fobj[i])
        # plt.plot(np.sort(np.array(Ttab_val_f).T[i]), pdf)

    # plt.show()
    print(result_dict_eigen)
    with open(fr"{filepath}\uq.json", 'w') as file:
        file.write(json.dumps(result_dict_eigen, indent=4, separators=(',', ': ')))


def results_from_qois(file_folder, eigen_obj_list):
    dir_list = os.listdir(file_folder)
    qois_result_dict = {}
    Ttab_val_f = []
    keys = []
    for dirs in dir_list:
        if os.path.isdir(fr'{file_folder}\{dirs}'):
            if 'C3794_n2' in dirs:
                keys.append(int(dirs.split('C3794_n2_Q')[-1]))
                with open(fr'{file_folder}\{dirs}\monopole\qois.json') as json_file:
                    qois_result_dict.update(json.load(json_file))

                qois_result = get_qoi_value(qois_result_dict, eigen_obj_list)
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

    data_table = pd.DataFrame(Ttab_val_f, columns=list(eigen_obj_list))
    data_table.index = keys
    # ic(data_table)
    data_table = data_table.sort_index()
    # ic(data_table)

    data_table.to_csv(fr'{file_folder}\table.csv', index=False, sep='\t')
    data_table.to_excel(fr'{file_folder}\table.xlsx', index=False)


def get_qoi_value(d, obj):
    """
    Gets the quantities of interest from simulation results
    Parameters
    ----------
    d: dict
        Dictionary containing several figures of merits from eigenmode solver
    obj: list
        List of objective functions
    n_cells: int
        Number of cells
    norm_length: float
        Normalisation length for :math: `E_\mathrm{acc}`

    Returns
    -------

    """
    # Req = d['CAVITY RADIUS'][n_cells - 1] * 10  # convert to mm
    # Freq = d['FREQUENCY'][n_cells - 1]
    # E_stored = d['STORED ENERGY'][n_cells - 1]
    # # Rsh = d['SHUNT IMPEDANCE'][n_cells-1]  # MOhm
    # Q = d['QUALITY FACTOR'][n_cells - 1]
    # Epk = d['MAXIMUM ELEC. FIELD'][n_cells - 1]  # MV/m
    # Hpk = d['MAXIMUM MAG. FIELD'][n_cells - 1]  # A/m
    # # Vacc = dict['ACCELERATION'][0]
    # # Eavg = d['AVERAGE E.FIELD ON AXIS'][n_cells-1]  # MV/m
    # Rsh_Q = d['EFFECTIVE IMPEDANCE'][n_cells - 1]  # Ohm
    #
    # Vacc = np.sqrt(
    #     2 * Rsh_Q * E_stored * 2 * np.pi * Freq * 1e6) * 1e-6
    # # factor of 2, remember circuit and accelerator definition
    # # Eacc = Vacc / (374 * 1e-3)  # factor of 2, remember circuit and accelerator definition
    # Eacc = Vacc / (norm_length * 1e-3)  # for 1 cell factor of 2, remember circuit and accelerator definition
    # Epk_Eacc = Epk / Eacc
    # Bpk_Eacc = (Hpk * 4 * np.pi * 1e-7) * 1e3 / Eacc
    #
    # d = {
    #     "Req": Req,
    #     "freq": Freq,
    #     "Q": Q,
    #     "E": E_stored,
    #     "R/Q": 2 * Rsh_Q,
    #     "Epk/Eacc": Epk_Eacc,
    #     "Bpk/Eacc": Bpk_Eacc
    # }

    objective = []

    # append objective functions
    for o in obj:
        if o in d.keys():
            objective.append(d[o])

    return objective


def serialise(state_dict, widget=None, visited_widgets=None, marker=''):
    """
    Serialise PyQt5 GUI widgets
    Parameters
    ----------
    state_dict: dict
        Python dictionary of state of GUI widgets
    widget: QWidget
        PyQt5 GUI widget to be serialised
    visited_widgets: list
        List of QWidgets already serialised
    marker: str
        Unique marker for QWidget to be serialised

    Returns
    -------

    """
    if visited_widgets is None:
        visited_widgets = set()

    if widget in visited_widgets:
        return 0

    visited_widgets.add(widget)

    if widget.objectName() != "" and 'qt_' not in widget.objectName():
        if isinstance(widget, QPushButton):
            if widget.isCheckable():
                state_dict[f'{marker}_' + widget.objectName()] = widget.isChecked()

        if isinstance(widget, QComboBox):
            state_dict[f'{marker}_' + widget.objectName()] = widget.currentText()

        if isinstance(widget, QLineEdit):
            state_dict[f'{marker}_' + widget.objectName()] = widget.text()

        if isinstance(widget, QSpinBox):
            state_dict[f'{marker}_' + widget.objectName()] = widget.value()

        if isinstance(widget, QDoubleSpinBox):
            state_dict[f'{marker}_' + widget.objectName()] = widget.value()

        if isinstance(widget, QCheckBox):
            state_dict[f'{marker}_' + widget.objectName()] = widget.checkState()

    for child_widget in widget.findChildren(QWidget):
        serialise(state_dict, child_widget, visited_widgets, marker)


def deserialise(state_dict, widget, visited_widgets=None, marker=''):
    """
    Deserialise PyQt5 GUI widgets
    Parameters
    ----------
    state_dict: dict
        Python dictionary of state of GUI widgets
    widget: QWidget
        PyQt5 GUI widget to be deserialised
    visited_widgets: list
        List of QWidgets already deserialised
    marker: str
        Unique marker for QWidget to be deserialised

    Returns
    -------

    """
    try:
        if visited_widgets is None:
            visited_widgets = set()

        if widget in visited_widgets:
            return

        visited_widgets.add(widget)

        if widget.objectName() != "" and 'qt_' not in widget.objectName():
            if isinstance(widget, QPushButton):
                if widget.isCheckable():
                    if not widget.isChecked() == state_dict[f'{marker}_' + widget.objectName()]:
                        widget.toggle()

            if isinstance(widget, QComboBox):
                widget.setCurrentText(state_dict[f'{marker}_' + widget.objectName()])

                # if isinstance(widget, QCheckableComboBox):
                try:
                    selection = state_dict[f'{marker}_' + widget.objectName()].split(', ')
                    # try to select copied selection
                    for txt in selection:
                        if txt in widget.texts:
                            item = widget.model().item(widget.texts.index(txt))
                            item.setCheckState(2)
                except AttributeError:
                    pass

            if isinstance(widget, QLineEdit):
                widget.setText(state_dict[f'{marker}_' + widget.objectName()])

            if isinstance(widget, QSpinBox):
                widget.setValue(int(state_dict[f'{marker}_' + widget.objectName()]))

            if isinstance(widget, QDoubleSpinBox):
                widget.setValue(float(state_dict[f'{marker}_' + widget.objectName()]))

            if isinstance(widget, QCheckBox):
                widget.setCheckState(state_dict[f'{marker}_' + widget.objectName()])

        for child_widget in widget.findChildren(QWidget):
            deserialise(state_dict, child_widget, visited_widgets, marker)

    except (TypeError, KeyError) as e:
        print(f"Could not deserialize {widget.objectName()}: ", e)


if __name__ == '__main__':
    # nodes, weights, bpoly = quad_stroud3(5, 1)
    # ic(weights)
    # ic(2. * nodes.T - 1)
    # ic(nodes.shape)
    # #
    # ic()
    # nodes, weights  = cn_leg_03_1(5)
    # ic(weights)
    # ic(nodes.T)
    # ic(np.shape(nodes))
    #
    # ic()
    # nodes, weights  = cn_leg_03_xiu(5)
    # ic(weights)
    # ic(nodes.T)
    # ic(nodes.shape)
    #
    # ic()
    # nodes, weights = cn_gauss(6, 2)
    # ic(weights)
    # ic(nodes.T)
    # ic(nodes.shape)
    # nodes, weights  = cn_leg_05_2(5)
    # ic(weights)
    # ic(nodes)

    # folder = r'D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\NGSolveMEVP\C3794_n2'
    # # results_from_qois(fr'{folder}',
    # #                   ["freq [MHz]", "R/Q [Ohm]", "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "G [Ohm]", "kcc [%]", "ff [%]"])
    # combine_params_output(folder, 30)
    # uq_from_excel(fr'{folder}')

    file_path = fr'D:\Dropbox\CavityDesignHub\utils\geodata.n'
    l_end_cell = np.array([73.52, 131.75, 106.25, 118.7, 150, 187, 350])
    # #
    midcell = np.array([[[60, 60], [56, 54], [76, 56], [76, 56]],  # <- A
                         [[60, 43], [73, 65], [65, 45], [65, 45]],  # <- B
                         [[34, 50], [54, 50], [37, 40], [37, 40]],  # <- a
                         [[40, 36], [56, 57], [54, 56], [54, 56]],   # <- b
                         [[150, 151], [170, 165], [150, 155], [150, 155]],   # <- Ri
                         [[187, 187], [184, 176], [178, 170], [178, 170]],  # <- L
                         [[369.6321578127116, 345], [340, 350], [350, 360], [350, 360]]])

    r_end_cell = np.array([70, 120, 100, 110, 130, 176, 340])
    # l_end_cell = np.array([73.58423524, 132.0824365, 106.4198031, 118.3902662, 149.584055, 186.7246878, 369.6661589])
    # # #
    # midcell = np.array([[[73.14364045, 73.05957235]],  # <- A
    #                      [[132.0125027, 131.8810112]],  # <- B
    #                      [[106.1869176, 106.4751203]],  # <- a
    #                      [[118.3555127, 118.8357715]],   # <- b
    #                      [[149.7665359, 149.9423194]],   # <- Ri
    #                      [[186.7982958, 186.5009973]],  # <- L
    #                      [[369.9589997, 369.7163648]]])
    # #
    # r_end_cell = np.array([73.52142396, 131.329964, 106.1202883, 119.1231179, 150.6952404, 186.9420167, 370.0226815])

    writeCavityForMultipac_multicell(file_path, 5, midcell, l_end_cell, r_end_cell, beampipe='both',
                                     plot=True, unit=1e-3, scale=1)
    #
    # nodes_ = pd.read_csv(fr'C:\Users\sosoho\DakotaProjects\Cavity\C3794_lhs\sim_result_table.dat', sep='\s+').iloc[:,
    #          2:-2]
    # nodes_ = nodes_.to_numpy().T
    # ic(nodes_)
    # midcell = np.atleast_2d(np.array([62.22222222222222, 66.12612612612612, 30.22022022022022, 23.113113113113116,
    #                      71.98698698698699, 93.5, 171.1929])).T*1e-3
    # # ic(midcell)
    # x = create_multicell_random_variables(3, midcell, midcell, midcell)
    # ic(x)

    # moved_val_indx = 7
    # proc_nodes_len = 25
    # nodes = nodes[:, 14]
    # ic(nodes, nodes.shape)
    # for i2 in range(2 * 2 - 1):
    #     if i2 % 2 == 0:
    #         # print("num", proc_nodes_len - moved_val_indx, "->", proc_nodes_len - moved_val_indx + 7)
    #         nodes = np.insert(nodes, proc_nodes_len - moved_val_indx + 7,
    #                                nodes[proc_nodes_len - moved_val_indx])
    #         # update index
    #         moved_val_indx += 7
    #     else:
    #         # print("num", proc_nodes_len - moved_val_indx, "->", proc_nodes_len - moved_val_indx + 6)
    #         nodes = np.insert(nodes, proc_nodes_len - moved_val_indx + 6,
    #                                nodes[proc_nodes_len - moved_val_indx])
    #         moved_val_indx += 5
    # ic(nodes, nodes.shape)

    # n_cells = 5
    # cav_var_list = ['A', 'B', 'a', 'b', 'Ri', 'L', 'Req']
    # midcell_var_dict = dict()
    # for i1 in range(len(midcell)):
    #     for i2 in range(n_cells-1):
    #         # if i2 == 0:
    #         #     midcell_var_dict[f'{cav_var_list[i1]}_el'] = [i1, i2]
    #         for i3 in range(2):
    #             midcell_var_dict[f'{cav_var_list[i1]}_{i2}_m{i3}'] = [i1, i2, i3]
    #         # if i2 == n_cells-1:
    #         #     midcell_var_dict[f'{cav_var_list[i1]}_er'] = [i1, i2]
    #
    # ic(midcell_var_dict.keys())
    # ic(len(midcell_var_dict.keys()))
