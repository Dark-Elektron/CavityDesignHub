import matplotlib.pyplot as plt
import numpy as np
from icecream import ic
from scipy.optimize import fsolve


def drawCavity():

    # # TESLA end cell 1
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 42, 42, 12, 19, 35, 57.6524, 103.353
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 40.34, 40.34, 10, 13.5, 39, 55.716, 103.353
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 40.34, 40.34, 10, 13.5, 39, 55.716, 103.353

    midC3795 = np.array([62.22222222222222, 66.12612612612612, 30.22022022022022, 23.113113113113116,
                         71.98698698698699, 93.5, 171.1929])*1e-3
    endC3795 = np.array([62.58258258258258, 57.53753753753754, 17.207207207207208, 12.002002002002001,
                         80.38038038038039, 93.31191678718535, 171.1929]) * 1e-3

    midFCCUROS5 = np.array([67.72, 57.45, 21.75, 35.95, 60, 93.5, 166.591])*1e-3
    endFCCUROS5 = np.array([66.5, 51, 17, 23, 78, 85.77, 166.591]) * 1e-3

    midTESLA = np.array([68.12, 68.12, 19.46, 30.81, 56.76, 93.5, 167.62]) * 1e-3
    endTESLA_l = np.array([65.42, 65.42, 16.22, 21.89, 63.25, 90.36, 167.62]) * 1e-3
    endTESLA_r = np.array([68.12, 68.12, 14.60, 20.76, 63.25, 92.14, 167.62]) * 1e-3

    # # TESLA end cell 2
    A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = midTESLA
    A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = endTESLA_r
    A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = endTESLA_l

    n_cell = 5
    step = 2  # step in boundary points in mm
    L_bp_l = 4 * L_m  # 0.0001  #
    L_bp_r = 4 * L_m  # 0.0001  #

    # calculate shift
    shift = (L_bp_r + L_bp_l + L_el + (n_cell - 1) * 2 * L_m + L_er) / 2
    # shift = 0
    # shift = L_m  # for end cell

    # calculate angles outside loop
    # CALCULATE x1_el, y1_el, x2_el, y2_el
    data = ([0 + L_bp_l, Ri_el + b_el, L_el + L_bp_l, Req_el - B_el],
            [a_el, b_el, A_el, B_el])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])

    x1el, y1el, x2el, y2el = fsolve(f, np.array(
        [a_el + L_bp_l, Ri_el + 0.85 * b_el, L_el - A_el + L_bp_l, Req_el - 0.85 * B_el]),
                                    args=data,
                                    xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

    # CALCULATE x1, y1, x2, y2
    data = ([0 + L_bp_l, Ri_m + b_m, L_m + L_bp_l, Req_m - B_m],
            [a_m, b_m, A_m, B_m])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    x1, y1, x2, y2 = fsolve(f, np.array([a_m + L_bp_l, Ri_m + 0.85 * b_m, L_m - A_m + L_bp_l, Req_m - 0.85 * B_m]),
                            args=data, xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

    # CALCULATE x1_er, y1_er, x2_er, y2_er
    data = ([0 + L_bp_r, Ri_er + b_er, L_er + L_bp_r, Req_er - B_er],
            [a_er, b_er, A_er, B_er])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    x1er, y1er, x2er, y2er = fsolve(f, np.array(
        [a_er + L_bp_r, Ri_er + 0.85 * b_er, L_er - A_er + L_bp_r, Req_er - 0.85 * B_er]),
                                    args=data,
                                    xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

    with open(r'D:\Dropbox\multipacting\MPGUI21\geodata.n', 'w') as fil:
        fil.write("   2.0000000e-03   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")
        fil.write("   1.25000000e-02   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")  # a point inside the structure
        fil.write("  -3.1415927e+00  -2.7182818e+00   0.0000000e+00   0.0000000e+00\n")  # a point outside the structure

        # SHIFT POINT TO START POINT
        start_point = [-shift, 0]
        fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   3.0000000e+00   0.0000000e+00\n")

        lineTo(start_point, [-shift, Ri_el], step)
        pt = [-shift, Ri_el]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        # ADD BEAM PIPE LENGTH
        lineTo(pt, [L_bp_l - shift, Ri_el], step)
        pt = [L_bp_l - shift, Ri_el]
        print(pt)
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        for n in range(1, n_cell + 1):
            if n == 1:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el])
                pt = [-shift + x1el, y1el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2el, y2el], step)
                pt = [-shift + x2el, y2el]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_el + L_bp_l - shift, Req_el - B_el, A_el, B_el, step, pt, [L_bp_l + L_el - shift, Req_el])
                pt = [L_bp_l + L_el - shift, Req_el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                if n_cell == 1:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l - shift, Req_er - B_er, A_er, B_er, step, [pt[0], Req_er - B_er],
                                [L_el + L_er - x2er + 2 * L_bp_l - shift, Req_er])
                    pt = [L_el + L_er - x2er + 2 * L_bp_l - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_el + L_er - x1er + 2 * L_bp_l - shift, y1er], step)
                    pt = [L_el + L_er - x1er + 2 * L_bp_l - shift, y1er]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                [L_bp_l + L_el + L_er - shift, y1er])

                    pt = [L_bp_l + L_el + L_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")

                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_er)
                    ic(shift)
                else:
                    print("if else")
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, [pt[0], Req_m - B_m],
                                [L_el + L_m - x2 + 2 * L_bp_l - shift, Req_m])
                    pt = [L_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_el + L_m - x1 + 2 * L_bp_l - shift, y1], step)
                    pt = [L_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                [L_bp_l + L_el + L_m - shift, y1])
                    pt = [L_bp_l + L_el + L_m - shift, Ri_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_m)
                    ic(shift)

            elif n > 1 and n != n_cell:
                print("elif")
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2, y2], step)
                pt = [-shift + x2, y2]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req_m])
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
                            [L_m + L_m - x2 + 2 * L_bp_l - shift, Req_m])
                pt = [L_m + L_m - x2 + 2 * L_bp_l - shift, y2]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                lineTo(pt, [L_m + L_m - x1 + 2 * L_bp_l - shift, y1], step)
                pt = [L_m + L_m - x1 + 2 * L_bp_l - shift, y1]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                            [L_bp_l + L_m + L_m - shift, y1])
                pt = [L_bp_l + L_m + L_m - shift, Ri_m]
                ic(pt)
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # calculate new shift
                shift = shift - 2*L_m
            else:
                print("else")
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2, y2], step)
                pt = [-shift + x2, y2]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req_m])
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
                            [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, Req_er])
                pt = [L_m + L_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                lineTo(pt, [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er], step)
                pt = [L_m + L_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                            [L_bp_l + L_m + L_er - shift, y1er])
                pt = [L_bp_l + L_m + L_er - shift, Ri_er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        # BEAM PIPE
        # reset shift
        print("pt before", pt)
        shift = (L_bp_r + L_bp_l + (n_cell - 1) * 2 * L_m + L_el + L_er) / 2
        lineTo(pt, [L_bp_r + L_bp_l + 2 * (n_cell-1) * L_m + L_el + L_er - shift, Ri_er], step)
        pt = [2 * (n_cell-1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, Ri_er]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   3.0000000e+00   0.0000000e+00\n")
        print("pt after", pt)

        # END PATH
        lineTo(pt, [2 * (n_cell-1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0], step)  # to add beam pipe to right
        pt = [2 * (n_cell-1) * L_m + L_el + L_er + L_bp_l + L_bp_r - shift, 0]
        # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
        # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

        # CLOSE PATH
        lineTo(pt, start_point, step)
        fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

    plt.show()


def drawCavity_flat_top():
    plt.rcParams["figure.figsize"] = (12, 3)
    midJlab = np.array([64.45, 54.58, 19.1, 25.92, 65, 83.554, 163.975, 90, 20]) * 1e-3  # [A, B, a, b, Ri, L, Req, alpha, l]
    endJlab = np.array([64.45, 54.58, 19.1, 25.92, 65, 83.554, 163.975, 90, 11.19]) * 1e-3
    endJlab_r = np.array([64.45, 54.58, 19.1, 25.92, 65, 83.554, 163.975, 90, 11.19]) * 1e-3

    # # TESLA end cell 2
    A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m, alpha, lft = midJlab
    A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el, alpha_el, lft_el = endJlab
    A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er, alpha_er, lft_er = endJlab

    n_cell = 5
    step = 2  # step in boundary points in mm
    L_bp_l = 4 * L_m  # 0.0001  #
    L_bp_r = 4 * L_m  # 0.0001  #

    # calculate shift
    shift = (L_bp_r + L_bp_l + L_el + lft_el + (n_cell - 1) * 2 * L_m + (n_cell - 2)*lft + L_er + lft_er) / 2
    # shift = 0
    # shift = L_m  # for end cell

    # calculate angles outside loop
    # CALCULATE x1_el, y1_el, x2_el, y2_el
    data = ([0 + L_bp_l, Ri_el + b_el, L_el + L_bp_l, Req_el - B_el],
            [a_el, b_el, A_el, B_el])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])

    x1el, y1el, x2el, y2el = fsolve(f, np.array(
        [a_el + L_bp_l, Ri_el + 0.85 * b_el, L_el - A_el + L_bp_l, Req_el - 0.85 * B_el]),
                                    args=data,
                                    xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

    # CALCULATE x1, y1, x2, y2
    data = ([0 + L_bp_l, Ri_m + b_m, L_m + L_bp_l, Req_m - B_m],
            [a_m, b_m, A_m, B_m])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    x1, y1, x2, y2 = fsolve(f, np.array([a_m + L_bp_l, Ri_m + 0.85 * b_m, L_m - A_m + L_bp_l, Req_m - 0.85 * B_m]),
                            args=data, xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

    # CALCULATE x1_er, y1_er, x2_er, y2_er
    data = ([0 + L_bp_r, Ri_er + b_er, L_er + L_bp_r, Req_er - B_er],
            [a_er, b_er, A_er, B_er])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    x1er, y1er, x2er, y2er = fsolve(f, np.array(
        [a_er + L_bp_r, Ri_er + 0.85 * b_er, L_er - A_er + L_bp_r, Req_er - 0.85 * B_er]),
                                    args=data,
                                    xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

    ic(alpha, alpha_el, alpha_er)
    with open(r'D:\Dropbox\multipacting\MPGUI21\geodata.n', 'w') as fil:
        fil.write("   2.0000000e-03   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")
        fil.write("   1.25000000e-02   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")  # a point inside the structure
        fil.write("  -3.1415927e+00  -2.7182818e+00   0.0000000e+00   0.0000000e+00\n")  # a point outside the structure

        # SHIFT POINT TO START POINT
        start_point = [-shift, 0]
        fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   3.0000000e+00   0.0000000e+00\n")

        lineTo(start_point, [-shift, Ri_el], step)
        pt = [-shift, Ri_el]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        # ADD BEAM PIPE LENGTH
        lineTo(pt, [L_bp_l - shift, Ri_el], step)
        pt = [L_bp_l - shift, Ri_el]
        print(pt)
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        for n in range(1, n_cell + 1):
            if n == 1:
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_el + b_el, a_el, b_el, step, pt, [-shift + x1el, y1el])
                pt = [-shift + x1el, y1el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2el, y2el], step)
                pt = [-shift + x2el, y2el]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_el + L_bp_l - shift, Req_el - B_el, A_el, B_el, step, pt, [L_bp_l + L_el - shift, Req_el])
                pt = [L_bp_l + L_el - shift, Req_el]
                for pp in pts:
                    fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # flat top
                lineTo(pt, [L_bp_l + L_el + lft_el - shift, Req_el], step)
                pt = [L_bp_l + L_el + lft_el - shift, Req_el]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                if n_cell == 1:
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l + lft_el - shift, Req_er - B_er, A_er, B_er, step, [pt[0], Req_er - B_er],
                                [L_el + lft_el + L_er - x2er + 2 * L_bp_l - shift, Req_er])
                    pt = [L_el + lft_el + L_er - x2er + 2 * L_bp_l - shift, y2er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_el + L_er - x1er + 2 * L_bp_l - shift, y1er], step)
                    pt = [L_el + L_er - x1er + 2 * L_bp_l - shift, y1er]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                                [L_bp_l + L_el + L_er - shift, y1er])

                    pt = [L_bp_l + L_el + L_er - shift, Ri_er]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")

                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_er + lft_el)
                    ic(shift)
                else:
                    print("if else")
                    # EQUATOR ARC TO NEXT POINT
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + L_bp_l + lft_el - shift, Req_m - B_m, A_m, B_m, step, [pt[0], Req_m - B_m],
                                [L_el + L_m + lft_el - x2 + 2 * L_bp_l - shift, Req_m])
                    pt = [L_el + lft_el + L_m - x2 + 2 * L_bp_l - shift, y2]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # STRAIGHT LINE TO NEXT POINT
                    lineTo(pt, [L_el + lft_el + L_m - x1 + 2 * L_bp_l - shift, y1], step)
                    pt = [L_el + lft_el + L_m - x1 + 2 * L_bp_l - shift, y1]
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # ARC
                    # half of bounding box is required,
                    # start is the lower coordinate of the bounding box and end is the upper
                    pts = arcTo(L_el + lft_el + L_m + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                                [L_bp_l + L_el + lft_el + L_m - shift, y1])
                    pt = [L_bp_l + L_el + lft_el + L_m - shift, Ri_m]
                    for pp in pts:
                        if (np.around(pp, 12) != np.around(pt, 12)).all():
                            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                        else:
                            print("Found one")
                    fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                    # calculate new shift
                    shift = shift - (L_el + L_m + lft_el)
                    ic(shift)

            elif n > 1 and n != n_cell:
                print("elif")
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2, y2], step)
                pt = [-shift + x2, y2]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req_m])
                pt = [L_bp_l + L_m - shift, Req_m]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # flat top
                lineTo(pt, [L_bp_l + L_m + lft - shift, Req_m], step)
                pt = [L_bp_l + L_el + lft - shift, Req_el]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_bp_l + lft - shift, Req_m - B_m, A_m, B_m, step, [pt[0], Req_m - B_m],
                            [L_m + L_m + lft - x2 + 2 * L_bp_l - shift, Req_m])
                pt = [L_m + L_m + lft - x2 + 2 * L_bp_l - shift, y2]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                lineTo(pt, [L_m + L_m + lft - x1 + 2 * L_bp_l - shift, y1], step)
                pt = [L_m + L_m + lft - x1 + 2 * L_bp_l - shift, y1]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_m + lft + L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, [pt[0], Ri_m],
                            [L_bp_l + L_m + L_m + lft - shift, y1])
                pt = [L_bp_l + L_m + L_m + lft - shift, Ri_m]
                ic(pt)
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # calculate new shift
                shift = shift - 2*L_m - lft
            else:
                print("else")
                # DRAW ARC:
                pts = arcTo(L_bp_l - shift, Ri_m + b_m, a_m, b_m, step, pt, [-shift + x1, y1])
                pt = [-shift + x1, y1]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW LINE CONNECTING ARCS
                lineTo(pt, [-shift + x2, y2], step)
                pt = [-shift + x2, y2]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # DRAW ARC, FIRST EQUATOR ARC TO NEXT POINT
                pts = arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, pt, [L_bp_l + L_m - shift, Req_m])
                pt = [L_bp_l + L_m - shift, Req_m]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # flat top
                lineTo(pt, [L_bp_l + L_m + lft_er - shift, Req_m], step)
                pt = [L_bp_l + L_el + lft_er - shift, Req_el]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # EQUATOR ARC TO NEXT POINT
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + lft_er + L_bp_l - shift, Req_er - B_er, A_er, B_er, step, [pt[0], Req_er - B_er],
                            [L_m + L_er + lft_er - x2er + L_bp_l + L_bp_r - shift, Req_er])
                pt = [L_m + L_er + lft_er - x2er + L_bp_l + L_bp_r - shift, y2er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # STRAIGHT LINE TO NEXT POINT
                lineTo(pt, [L_m + L_er + lft_er - x1er + L_bp_l + L_bp_r - shift, y1er], step)
                pt = [L_m + L_er + lft_er - x1er + L_bp_l + L_bp_r - shift, y1er]
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

                # ARC
                # half of bounding box is required,
                # start is the lower coordinate of the bounding box and end is the upper
                pts = arcTo(L_m + L_er + lft_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                            [L_bp_l + L_m + L_er + lft_er - shift, y1er])
                pt = [L_bp_l + L_m + L_er + lft_er - shift, Ri_er]
                for pp in pts:
                    if (np.around(pp, 12) != np.around(pt, 12)).all():
                        fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
                    else:
                        print("Found one")
                fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        # BEAM PIPE
        # reset shift
        print("pt before", pt)
        shift = (L_bp_r + L_bp_l + L_el + lft_el + (n_cell - 1) * 2 * L_m + (n_cell - 2)*lft + L_er + lft_er) / 2
        lineTo(pt, [L_bp_r + L_bp_l + 2 * (n_cell-1) * L_m + (n_cell-2)*lft + lft_el + lft_er + L_el + L_er - shift, Ri_er], step)

        pt = [2 * (n_cell-1) * L_m + L_el + L_er + L_bp_l + L_bp_r + (n_cell-2)*lft + lft_el + lft_er - shift, Ri_er]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   3.0000000e+00   0.0000000e+00\n")
        print("pt after", pt)

        # END PATH
        lineTo(pt, [2 * (n_cell-1) * L_m + L_el + L_er + (n_cell-2)*lft + lft_el + lft_er + L_bp_l + L_bp_r - shift, 0], step)  # to add beam pipe to right
        pt = [2 * (n_cell-1) * L_m + L_el + L_er + (n_cell-2)*lft + lft_el + lft_er + L_bp_l + L_bp_r - shift, 0]
        # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
        # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

        # CLOSE PATH
        lineTo(pt, start_point, step)
        fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

    plt.show()


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

    plt.plot(px, py)


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

    plt.plot(u + a * np.cos(t), v + b * np.sin(t))

    return [x_end, y_end]


def arcTo(x_center, y_center, a, b, step, start, end):
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

    plt.plot(inbox[:, 0], inbox[:, 1])

    return inbox

if __name__ == '__main__':

    # drawCavity()
    drawCavity_flat_top()
    # drawCapacitor()
