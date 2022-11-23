import matplotlib.pyplot as plt
import numpy as np
from icecream import ic
from scipy.optimize import fsolve


def writeCavity(key, mid_cell, lend_cell, rend_cell, beampipe):

    # # C4123
    A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = np.array(mid_cell)*1e-3
    A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = np.array(lend_cell)*1e-3
    A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = np.array(rend_cell)*1e-3

    n_cell = 1
    step = 2  # step in boundary points in mm
    L_bp_l, L_bp_r = beampipe

    # calculate shift
    shift = (L_bp_r + L_bp_l + (n_cell - 1) * 2 * L_m + L_el + L_er) / 2

    # calculate angles outside loop
    # CALCULATE x1_el, y1_el, x2_el, y2_el
    data = ([0 + L_bp_l, Ri_el + b_el, L_el + L_bp_l, Req_el - B_el],
            [a_el, b_el, A_el, B_el])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])

    x1el, y1el, x2el, y2el = fsolve(f, np.array(
        [a_el + L_bp_l, Ri_el + 0.85 * b_el, L_el - A_el + L_bp_l, Req_el - 0.85 * B_el]),
                                    args=data,
                                    xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

    ### CALCULATE x1, y1, x2, y2
    data = ([0 + L_bp_l, Ri_m + b_m, L_m + L_bp_l, Req_m - B_m],
            [a_m, b_m, A_m, B_m])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    x1, y1, x2, y2 = fsolve(f, np.array([a_m + L_bp_l, Ri_m + 0.85 * b_m, L_m - A_m + L_bp_l, Req_m - 0.85 * B_m]),
                            args=data, xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

    ### CALCULATE x1_er, y1_er, x2_er, y2_er
    data = ([0 + L_bp_l, Ri_er + b_er, L_er + L_bp_l, Req_er - B_er],
            [a_er, b_er, A_er, B_er])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    x1er, y1er, x2er, y2er = fsolve(f, np.array(
        [a_er + L_bp_l, Ri_er + 0.85 * b_er, L_er - A_er + L_bp_l, Req_er - 0.85 * B_er]),
                                    args=data,
                                    xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

    # with open(r'D:\Dropbox\multipacting\MPGUI21\geodata.n', 'w') as fil:
    with open(fr'D:\Dropbox\CEMCodesHub\C800MHz\PostprocessingData\Data\{key}_geom.txt', 'w') as fil:
        # SHIFT POINT TO START POINT
        start_point = [-shift, 0]
        fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}\n")

        lineTo(start_point, [-shift, Ri_el], step)
        pt = [-shift, Ri_el]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}\n")

        # ADD BEAM PIPE LENGTH
        lineTo(pt, [L_bp_l - shift, Ri_el], step)
        pt = [L_bp_l - shift, Ri_el]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}\n")

        # DRAW ARC:
        pts = arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el])
        pt = [-shift + x1el, y1el]
        for pp in pts:
            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}\n")
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}\n")

        # DRAW LINE CONNECTING ARCS
        lineTo(pt, [-shift + x2el, y2el], step)
        pt = [-shift + x2el, y2el]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}\n")

        # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
        pts = arcTo(L_el + L_bp_l - shift, Req_el - B_el, A_el, B_el, step, pt, [L_bp_l + L_el - shift, Req_el])
        pt = [L_bp_l + L_el - shift, Req_el]
        for pp in pts:
            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}\n")
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}\n")

        # EQUATOR ARC TO NEXT POINT
        # half of bounding box is required, start is the lower coordinate of the bounding box and end is the upper
        pts = arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, [pt[0], pt[1] - B_er],
                    [+ 2 * L_er - x2er + 2 * L_bp_l - shift, Req_er])
        pt = [+ 2 * L_er - x2er + 2 * L_bp_l - shift, y2er]
        for pp in pts:
            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}\n")
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}\n")

        # STRAIGHT LINE TO NEXT POINT
        lineTo(pt, [+ 2 * L_er - x1er + 2 * L_bp_l - shift, y1er], step)
        pt = [+ 2 * L_er - x1er + 2 * L_bp_l - shift, y1er]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}\n")

        # ARC
        # half of bounding box is required, start is the lower coordinate of the bounding box and end is the upper
        pts = arcTo(2 * L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                    [L_bp_l + L_el + L_er - shift, y1er])
        pt = [L_bp_l + L_el + L_er - shift, Ri_er]
        for pp in pts:
            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}\n")
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}\n")

        # BEAM PIPE
        lineTo(pt, [2 * n_cell * L_er + L_bp_l + L_bp_r - shift, Ri_er], step)
        pt = [2 * n_cell * L_er + L_bp_l + L_bp_r - shift, Ri_er]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}\n")

        # END PATH
        lineTo(pt, [2 * n_cell * L_er + L_bp_l + L_bp_r - shift, 0], step)  # to add beam pipe to right
        pt = [2 * n_cell * L_er + L_bp_l + L_bp_r - shift, 0]
        # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
        # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}\n")

        # CLOSE PATH
        lineTo(pt, start_point, step)
        fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}\n")


def f(z, *data):
    coord, dim = data
    h, k, p, q = coord
    a, b, A, B = dim
    x1, y1, x2, y2 = z

    f1 = (x1 - h) ** 2 / a ** 2 + (y1 - k) ** 2 / b ** 2 - 1
    f2 = (x2 - p) ** 2 / A ** 2 + (y2 - q) ** 2 / B ** 2 - 1
    f3 = A ** 2 * b ** 2 * (x1 - h) * (y2 - q) / (a ** 2 * B ** 2 * (x2 - p) * (y1 - k)) - 1
    f4 = -b ** 2 * (x1 - x2) * (x1 - h) / (a ** 2 * (y1 - y2) * (y1 - k)) - 1

    return f1, f2, f3, f4


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


def lineTo(prevPt, nextPt, step):
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


def arcTo2(x_center, y_center, a, b, step, start_angle, end_angle):
    u = x_center  # x-position of the center
    v = y_center  # y-position of the center
    a = a  # radius on the x-axis
    b = b  # radius on the y-axis
    sa = (start_angle / 360) * 2 * np.pi  # convert angle to radians
    ea = (end_angle / 360) * 2 * np.pi  # convert angle to radians

    if ea < sa:
        # end point of curve
        x_end, y_end = u + a * np.cos(sa), v + b * np.sin(sa)

        t = np.arange(ea, sa, np.pi / 100)
        # t = np.linspace(ea, sa, 100)
        # check if end angle is included, include if not
        if sa not in t:
            t = np.append(t, sa)
        t = t[::-1]
    else:
        # end point of curve
        x_end, y_end = u + a * np.cos(ea), v + b * np.sin(ea)

        t = np.arange(sa, ea, np.pi / 100)
        # t = np.linspace(ea, sa, 100)
        if ea not in t:
            t = np.append(t, ea)

    # print("t0 ", [(u + a * np.cos(t))[0], (v + b * np.sin(t))[0]])
    # ic([u + a * np.cos(t), v + b * np.sin(t)])
    # ic()

    return [x_end, y_end]


def arcTo(x_center, y_center, a, b, step, start, end):
    u = x_center  # x-position of the center
    v = y_center  # y-position of the center
    a = a  # radius on the x-axis
    b = b  # radius on the y-axis

    t = np.arange(0, 2 * np.pi, np.pi / 100)

    x = u + a * np.cos(t)
    y = v + b * np.sin(t)
    pts = np.column_stack((x, y))
    inidx = np.all(np.logical_and(np.array(start) < pts, pts < np.array(end)), axis=1)
    inbox = pts[inidx]
    inbox = inbox[inbox[:, 0].argsort()]

    return inbox


mid40866 = [65.0, 30, 25, 20, 60, 93.5, 160.014]
mid3794 = [36.76, 65.875, 53.125, 59.35, 75, 93.5, 185.7276]
midG6 = [51.0510633611297, 50.8492563995682, 37.94024533955164, 27.22177461866947, 73.83844219184923, 93.5, 173.2763]
g1_c10_ch_q3 = [54.879053291574124, 53.90857921086279, 33.15859388359441, 15.488719027059194, 81.21661594222239, 93.5, 172.559]
key = ['C40866', 'C3794_800MHz', "G6_C170_M", 'g1_c10_ch_q3']
mids = [mid40866, mid3794, midG6]
for i, mid in enumerate(mids):
    writeCavity(key[i], mid, mid, mid, [0.001, 0.001])