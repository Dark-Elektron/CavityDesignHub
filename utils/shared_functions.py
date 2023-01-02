import ast
import json
import os
from PyQt5.QtWidgets import QFileDialog
from icecream import ic
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
    A, B, a, b, Ri, L, Req, _ = cell
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

            # print_(shape_space)
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


def open_file(frame_control, le, cb):
    # clear combobox
    frame_control.ui.cb_Shape_Space_Keys.clear()
    frame_control.ui.cb_Shape_Space_Keys.addItem('All')
    # self.selected_keys.clear()

    filename, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Json Files (*.json)")
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


if __name__ == '__main__':
    nodes, weights, bpoly = quad_stroud3(5, 2)
    ic(weights)
    ic(nodes)

    # scheme = quadpy.cn.stroud_cn_5_2(5)
    # ic(scheme.points)
