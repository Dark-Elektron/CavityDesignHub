import ast
import json
import os
from PyQt5.QtWidgets import QFileDialog
from icecream import ic
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from scipy.optimize import fsolve
import numpy as np
from PyQt5.QtCore import *
from utils.file_reader import FileReader

fr = FileReader()

AN_DURATION = 250
global animation


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
        A, B, a, b, Ri, L, Req, _ = cell
    else:
        A, B, a, b, Ri, L, Req = cell
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
    data = ([0 + L_bp, Ri + b, L + L_bp, Req - B],
            [a, b, A, B])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])

    df = fsolve(ellipse_tangent,
                np.array([a + L_bp, Ri + 0.85 * b, L - A + L_bp, Req - 0.85 * B]),
                args=data, fprime=jac, xtol=1.49012e-12, full_output=True)

    error_msg = df[-2]
    if error_msg == 4:
        df = fsolve(ellipse_tangent, np.array([a + L_bp, Ri + 1.15 * b, L - A + L_bp, Req - 1.15 * B]),
                    args=data, fprime=jac, xtol=1.49012e-12, full_output=True)
        ic(df)

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
    df1_dy1 = - A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)**2)
    df1_dx2 = - A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p)**2 * (y1 - k))
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
    df4_dy1 = -b ** 2 * (x1 - x2) * (x1 - h) * ((y1 - y2) + (y1 - k)) / (a ** 2 * ((y1 - y2) * (y1 - k))**2)
    df4_dx2 = b ** 2 * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k))
    df4_dy2 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2)**2 * (y1 - k))

    J = [[df1_dx1, df1_dy1, df1_dx2, df1_dy2],
         [df2_dx1, df2_dy1, df2_dx2, df2_dy2],
         [df3_dx1, df3_dy1, df3_dx2, df3_dy2],
         [df4_dx1, df4_dy1, df4_dx2, df4_dy2]]

    return J


def perform_geometry_checks(par_mid, par_end):
    """
    Checks geometry to filter out low loss cavity geometries

    Parameters
    ----------
    par_mid: list, array like
        Mid cell parameters of cavity
    par_end: list, array like
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


def write_cst_paramters(key, ic_, oc, projectDir, cell_type):
    """
    Writes cavity geometric data that can be imported into CST Studio

    Parameters
    ----------
    key: str, int
        Cavity marker
    ic_: list, array like
        Inner cavity cell geometric variables
    oc: list, array like
        Outer cavity cell geometric variables
    projectDir: str
        Project directory
    cell_type: str
        Single cell or multicell

    Returns
    -------

    """
    ic_ = update_alpha(ic_)
    oc = update_alpha(oc)
    if cell_type is None:
        # print("Writing parameters to file")
        path = fr'{projectDir}/SimulationData/SLANS/{key}/{key}.txt'

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
        path = fr'{projectDir}/SimulationData/SLANS/{key}/{key}.txt'
        path_mc = fr'{projectDir}/SimulationData/SLANS/{key}/{key}_Multicell.txt'

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


def quad_stroud5(rdim, degree):
    # data for Stroud-5 quadrature in [0,1]**rdim
    # nodes and weights
    # o, nodestr, weights = quadpy.cn.stroud_cn_5_2(rdim)
    nodes = 0.5 * nodestr + 0.5
    weights = weights / (2 ** rdim)
    dummy, nnodes = np.size(nodes)

    # evaluation of Legendre polynomials
    bpoly = np.zeros((degree + 1, rdim, nnodes))
    for l in range(rdim):
        for j in range(nnodes):
            bpoly[1, l, j] = 1.
            bpoly[2, l, j] = nodestr(l, j)
            for i in range(1, degree):
                bpoly[i + 1, l, j] = ((2 * i - 1) * nodestr(l, j) * bpoly[i, l, j] - (i - 1) * bpoly[i - 1, l, j]) / i

    # standardisation of Legendre polynomials
    for i in range(1, degree + 1):
        bpoly[i, :, :] = bpoly[i, :, :] * np.sqrt(2 * i - 1)

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
    global animation
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
            animation = QPropertyAnimation(widget, b"maximumWidth")
        else:
            animation = QPropertyAnimation(widget, b"minimumWidth")

        animation.setDuration(AN_DURATION)
        animation.setStartValue(width)
        animation.setEndValue(widthCollapsed)
        animation.setEasingCurve(QEasingCurve.InOutQuart)
        animation.start()


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

    path = path.replace(r"/", "\\")
    return path


def load_shape_space(filename):
    """
    Loads shape space file to Python dictionary object

    Appears to have no usage

    Parameters
    ----------
    filename: str
        Absolute path of shape space file

    Returns
    -------
    Python dictionary object of shape space

    """
    directory = filename

    # check if extension is included
    if directory.split('.')[-1] != 'json':
        directory = f'{dir}.json'

    df = fr.json_reader(directory)

    return df.to_dict()


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


def get_geometric_parameters(frame_control, code):
    shape_space = {}
    if frame_control.ui.cb_Shape_Entry_Mode.currentIndex() == 0:
        try:
            # get selected keys
            frame_control.selected_keys = frame_control.ui.cb_Shape_Space_Keys.currentText()
            print("selected keys: ", frame_control.ui.cb_Shape_Space_Keys.currentText())
            # print("Selected keys: ", frame_control.selected_keys, type(frame_control.selected_keys[0]))

            # check keys of shape space if results already exist
            to_all = None
            for key, val in frame_control.loaded_shape_space.items():
                # process for only keys selected in combobox
                if frame_control.ui.cb_Shape_Space_Keys.currentText() == "" \
                        or frame_control.ui.cb_Shape_Space_Keys.currentText() == "All":
                    pass
                else:
                    if isinstance(frame_control.selected_keys, str):
                        if key != frame_control.selected_keys:
                            continue
                    else:
                        if key not in frame_control.selected_keys:
                            continue

                if not to_all:
                    ans = frame_control.prompt(code, key)
                    if ans == 'Yes':
                        shape_space[key] = val

                    if ans == 'No':
                        continue

                    if ans == 'YesToAll':
                        shape_space[key] = val
                        to_all = 'YesToAll'

                    if ans == 'NoToAll':
                        to_all = 'NoToAll'

                    if ans == "Does not exist":
                        shape_space[key] = val
                        to_all = None
                else:
                    if to_all == 'YesToAll':
                        shape_space[key] = val
                    else:
                        path = f'{frame_control.main_control.projectDir}/SimulationData/{code}/{key}'
                        if os.path.exists(path):
                            continue
                        else:
                            shape_space[key] = val

            return shape_space
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

            inner_cell_space = [A_i_space, B_i_space, a_i_space, b_i_space, Ri_i_space, L_i_space, Req_i_space, alpha_i_space]
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

            outer_cell_L_space = [A_ol_space, B_ol_space, a_ol_space, b_ol_space, Ri_ol_space, L_ol_space, Req_ol_space, alpha_ol_space]
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

            outer_cell_R_space = [A_or_space, B_or_space, a_or_space, b_or_space, Ri_or_space, L_or_space, Req_or_space, alpha_or_space]
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
                                            shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_L, 'BP': 'both', 'FREQ': None}
                                        elif frame_control.ui.cb_LBP.checkState() == 2 and frame_control.ui.cb_RBP.checkState() == 0:
                                            shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_L, 'BP': 'left', 'FREQ': None}
                                        elif frame_control.ui.cb_LBP.checkState() == 0 and frame_control.ui.cb_RBP.checkState() == 2:
                                            shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_L, 'BP': 'right', 'FREQ': None}
                                        else:
                                            shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_L, 'BP': 'none', 'FREQ': None}

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
                                                                    inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, 0]
                                                                    outer_cell_L = [A_ol, B_ol, a_ol, b_ol, Ri_ol, L_ol, Req_i, 0]
                                                                    outer_cell_R = outer_cell_L
                                                                    if frame_control.ui.cb_LBP.checkState() == 2 and frame_control.ui.cb_RBP.checkState() == 0:
                                                                        shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'left', 'FREQ': None}
                                                                    elif frame_control.ui.cb_LBP.checkState() == 0 and frame_control.ui.cb_RBP.checkState() == 2:
                                                                        shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'right', 'FREQ': None}
                                                                    elif frame_control.ui.cb_LBP.checkState() == 2 and frame_control.ui.cb_RBP.checkState() == 2:
                                                                        shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'both', 'FREQ': None}
                                                                    else:
                                                                        shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'none', 'FREQ': None}

                                                                    count += 1
                                                                else:
                                                                    for A_or in outer_cell_R_space[0]:
                                                                        for B_or in outer_cell_R_space[1]:
                                                                            for a_or in outer_cell_R_space[2]:
                                                                                for b_or in outer_cell_R_space[3]:
                                                                                    for Ri_or in outer_cell_R_space[4]:
                                                                                        for L_or in outer_cell_R_space[5]:
                                                                                            # for Req_or in outer_cell_R_space[6]:
                                                                                            inner_cell = [A_i, B_i, a_i, b_i, Ri_i, L_i, Req_i, 0]
                                                                                            outer_cell_L = [A_ol, B_ol, a_ol, b_ol, Ri_ol, L_ol, Req_i, 0]
                                                                                            outer_cell_R = [A_or, B_or, a_or, b_or, Ri_or, L_or, Req_i, 0]
                                                                                            if frame_control.ui.cb_LBP.checkState() == 2 and frame_control.ui.cb_RBP.checkState() == 0:
                                                                                                shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'left', 'FREQ': None}
                                                                                            elif frame_control.ui.cb_LBP.checkState() == 0 and frame_control.ui.cb_RBP.checkState() == 2:
                                                                                                shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'right', 'FREQ': None}
                                                                                            elif frame_control.ui.cb_LBP.checkState() == 2 and frame_control.ui.cb_RBP.checkState() == 2:
                                                                                                shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'both', 'FREQ': None}
                                                                                            else:
                                                                                                shape_space[count] = {'IC': inner_cell, 'OC': outer_cell_L, 'OC_R': outer_cell_R, 'BP': 'none', 'FREQ': None}

                                                                                            count += 1
        return shape_space

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


def open_file(frame_control, le, cb, start_folder=''):
    # clear combobox
    frame_control.ui.cb_Shape_Space_Keys.clear()
    frame_control.ui.cb_Shape_Space_Keys.addItem('All')
    # self.selected_keys.clear()

    filename, _ = QFileDialog.getOpenFileName(None, "Open File", start_folder, "Json Files (*.json)")
    try:
        le.setText(filename)
        with open(filename, 'r') as file:
            dd = json.load(file)

        # populate checkboxes with key
        for col in dd.keys():
            cb.addItem(fr'{col}')

        frame_control.loaded_shape_space = dd

    except Exception as e:
        print('Failed to open file:: ', e)


def uq(key, shape, qois, n_cells, n_modules, n_modes, f_shift, bc, parentDir, projectDir):
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
            norm_length = 2*n_cells*shape['IC'][5]

            qois_result = get_qoi_value(params, slans_obj_list, n_cells, norm_length)
            print_(qois_result)
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
        print_(fr"There was a problem running UQ analysis for {key}")


def write_cavity_for_custom_eig_solver(file_path, n_cell, mid_cell, end_cell_left=None, end_cell_right=None,
                                       beampipe='none', plot=False):
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
    if len(mid_cell) == 7:
        A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = mid_cell
    else:
        A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, _ = mid_cell

    if len(mid_cell) == 7:
        A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = end_cell_left
    else:
        A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, _ = end_cell_left

    if len(mid_cell) == 7:
        A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = end_cell_right
    else:
        A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, _ = end_cell_right

    step = 2  # step in boundary points in mm

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

    # calculate angles outside loop
    # CALCULATE x1_el, y1_el, x2_el, y2_el
    df = tangent_coords(A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, L_bp_l)
    x1el, y1el, x2el, y2el = df[0]

    # CALCULATE x1, y1, x2, y2
    df = tangent_coords(A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, L_bp_l)
    x1, y1, x2, y2 = df[0]

    # CALCULATE x1_er, y1_er, x2_er, y2_er
    df = tangent_coords(A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, L_bp_r)
    x1er, y1er, x2er, y2er = df[0]

    with open(file_path, 'w') as fil:
        fil.write("   2.0000000e-03   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")
        fil.write("   1.25000000e-02   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")  # a point inside the structure
        fil.write("  -3.1415927e+00  -2.7182818e+00   0.0000000e+00   0.0000000e+00\n")  # a point outside the structure

        # SHIFT POINT TO START POINT
        start_point = [-shift, 0]
        fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   3.0000000e+00   0.0000000e+00\n")

        lineTo(start_point, [-shift, Ri_el], step, plot)
        pt = [-shift, Ri_el]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        # ADD BEAM PIPE LENGTH
        lineTo(pt, [L_bp_l - shift, Ri_el], step, plot)
        pt = [L_bp_l - shift, Ri_el]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        for n in range(1, n_cell + 1):
            if n == 1:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el], plot)
                pt = [-shift + x1el, y1el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2el, y2el], step, plot)
                pt = [-shift + x2el, y2el]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_el + L_bp_l - shift, Req_el - B_el, A_el, B_el, step, pt, [L_bp_l + L_el - shift, Req_el], plot)
                pt = [L_bp_l + L_el - shift, Req_el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                if n_cell == 1:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l - shift, Req_er - B_er, A_er, B_er, step, [pt[0], Req_er - B_er],
                                [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, Req_er], plot)
                    pt = [L_el + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_el + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step, plot)
                    pt = [L_el + L_er - x1er + + L_bp_l + L_bp_r  - shift, y1er]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                [L_bp_l + L_el + L_er - shift, y1er], plot)

                    pt = [L_bp_l + L_el + L_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")

                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_er)
                else:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, [pt[0], Req_m - B_m],
                                [L_el + L_m - x2 + 2 * L_bp_l - shift, Req_m], plot)
                    pt = [L_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_el + L_m - x1 + 2 * L_bp_l - shift, y1], step, plot)
                    pt = [L_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                [L_bp_l + L_el + L_m - shift, y1], plot)
                    pt = [L_bp_l + L_el + L_m - shift, Ri_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_m)

            elif n > 1 and n != n_cell:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req_m], plot)
                pt = [L_bp_l + L_m - shift, Req_m]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, [pt[0], Req_m - B_m],
                            [L_m + L_m - x2 + 2 * L_bp_l - shift, Req_m], plot)
                pt = [L_m + L_m - x2 + 2 * L_bp_l - shift, y2]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                lineTo(pt, [L_m + L_m - x1 + 2 * L_bp_l - shift, y1], step, plot)
                pt = [L_m + L_m - x1 + 2 * L_bp_l - shift, y1]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                            [L_bp_l + L_m + L_m - shift, y1], plot)
                pt = [L_bp_l + L_m + L_m - shift, Ri_m]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # calculate new shift
                shift = shift - 2*L_m
            else:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1], plot)
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2, y2], step, plot)
                pt = [-shift + x2, y2]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req_m], plot)
                pt = [L_bp_l + L_m - shift, Req_m]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l - shift, Req_er - B_er, A_er, B_er, step, [pt[0], Req_er - B_er],
                            [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, Req_er], plot)
                pt = [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                lineTo(pt, [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step, plot)
                pt = [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                            [L_bp_l + L_m + L_er - shift, y1er], plot)
                pt = [L_bp_l + L_m + L_er - shift, Ri_er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        # BEAM PIPE
        # reset shift
        shift = (L_bp_r + L_bp_l + (n_cell - 1) * 2 * L_m + L_el + L_er) / 2
        lineTo(pt, [L_bp_r + L_bp_l + 2 * (n_cell-1) * L_m + L_el + L_er - shift, Ri_er], step, plot)
        pt = [2 * (n_cell-1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, Ri_er]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   3.0000000e+00   0.0000000e+00\n")

        # END PATH
        lineTo(pt, [2 * (n_cell-1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0], step, plot)  # to add beam pipe to right
        pt = [2 * (n_cell-1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0]
        # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
        # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

        # CLOSE PATH
        lineTo(pt, start_point, step, plot)
        fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

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
        ll = np.linspace(start, stop, int((stop - start) / abs(step) + 1))
        if stop not in ll:
            ll = np.append(ll, stop)

        return ll
    else:
        ll = np.linspace(stop, start, int((start - stop) / abs(step) + 1))
        if start not in ll:
            ll = np.append(ll, start)
        return ll


def lineTo(prevPt, nextPt, step, plot=False):
    if prevPt[0] == nextPt[0]:
        # vertical line
        # chwxk id nextPt is greater
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
        plt.plot(px, py)


def arcTo2(x_center, y_center, a, b, step, start_angle, end_angle, plot=False):
    u = x_center  # x-position of the center
    v = y_center  # y-position of the center
    a = a  # radius on the x-axis
    b = b  # radius on the y-axis
    sa = (start_angle / 360) * 2 * np.pi  # convert angle to radians
    ea = (end_angle / 360) * 2 * np.pi  # convert angle to radians

    if ea < sa:
        # end point of curve
        x_end, y_end = u + a * np.cos(sa), v + b * np.sin(sa)

        t = np.arange(ea, sa, np.pi / 70)
        # t = np.linspace(ea, sa, 100)
        # check if end angle is included, include if not
        if sa not in t:
            t = np.append(t, sa)
        t = t[::-1]
    else:
        # end point of curve
        x_end, y_end = u + a * np.cos(ea), v + b * np.sin(ea)

        t = np.arange(sa, ea, np.pi / 70)
        # t = np.linspace(ea, sa, 100)
        if ea not in t:
            t = np.append(t, ea)

    # print("t0 ", [(u + a * np.cos(t))[0], (v + b * np.sin(t))[0]])
    # ic([u + a * np.cos(t), v + b * np.sin(t)])
    # ic()
    if plot:
        plt.plot(u + a * np.cos(t), v + b * np.sin(t))

    return [x_end, y_end]


def arcTo(x_center, y_center, a, b, step, start, end, plot=False):
    u = x_center  # x-position of the center
    v = y_center  # y-position of the center
    a = a  # radius on the x-axis
    b = b  # radius on the y-axis

    t = np.arange(0, 2 * np.pi, np.pi / 70)

    x = u + a * np.cos(t)
    y = v + b * np.sin(t)
    pts = np.column_stack((x, y))
    inidx = np.all(np.logical_and(np.array(start) < pts, pts < np.array(end)), axis=1)
    inbox = pts[inidx]
    inbox = inbox[inbox[:, 0].argsort()]

    if plot:
        plt.plot(inbox[:, 0], inbox[:, 1])

    return inbox


if __name__ == '__main__':
    nodes, weights, bpoly = quad_stroud3(5, 2)
    ic(weights)
    ic(nodes)

    # scheme = quadpy.cn.stroud_cn_5_2(5)
    # ic(scheme.points)
