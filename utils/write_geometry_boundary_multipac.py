import matplotlib.pyplot as plt
import numpy as np
from icecream import ic
from scipy.optimize import fsolve


def drawCavity():
    #
    #
    # # cornell erl
    # #### DEFINE VARIABLES
    # A_m =   43.99*1e-3
    # B_m =   35.06*1e-3
    # a_m =   12.53*1e-3
    # b_m =   20.95*1e-3
    # Ri_m =  35*1e-3
    # L_m =   57.6524*1e-3
    # Req_m = 101.205*1e-3
    #
    # ### left end cell
    # A_el =      43.99*1e-3
    # B_el =      35.06*1e-3
    # a_el =      12.53*1e-3
    # b_el =      20.95*1e-3
    # Ri_el =     35*1e-3
    # L_el =      57.6524*1e-3
    # Req_el =    101.205*1e-3
    #
    # ### right end cell
    # A_er =      43.99*1e-3
    # B_er =      35.06*1e-3
    # a_er =      12.53*1e-3
    # b_er =      20.95*1e-3
    # Ri_er =     35*1e-3
    # L_er =      57.6524*1e-3
    # Req_er =    101.205*1e-3

    # # FCCUROS Z working point 1 cell
    # #### DEFINE VARIABLES
    # A_m =   70*1e-3
    # B_m =   70*1e-3
    # a_m =   25*1e-3
    # b_m =   25*1e-3
    # Ri_m =  156*1e-3
    # L_m =   120*1e-3
    # Req_m = 350.574*1e-3
    #
    # ### left end cell
    # A_el =   70*1e-3
    # B_el =   70*1e-3
    # a_el =   25*1e-3
    # b_el =   25*1e-3
    # Ri_el =  156*1e-3
    # L_el =   120*1e-3
    # Req_el = 350.574*1e-3
    #
    # ### right end cell
    # A_er =      70*1e-3
    # B_er =      70*1e-3
    # a_er =      25*1e-3
    # b_er =      25*1e-3
    # Ri_er =     156*1e-3
    # L_er =      120*1e-3
    # Req_er =    350.574*1e-3

    # # C1092 Z working point 1 cell
    # #### DEFINE VARIABLES
    # A_m =   73.35*1e-3
    # B_m =   55.39*1e-3
    # a_m =   100*1e-3
    # b_m =   71.68*1e-3
    # Ri_m =  150*1e-3
    # L_m =   187*1e-3
    # Req_m = 349.0943*1e-3
    #
    # ### left end cell
    # A_el =   73.35*1e-3
    # B_el =   55.39*1e-3
    # a_el =   100*1e-3
    # b_el =   71.68*1e-3
    # Ri_el =  150*1e-3
    # L_el =   187*1e-3
    # Req_el = 349.0943*1e-3
    #
    # ### right end cell
    # A_er =      73.35*1e-3
    # B_er =      55.39*1e-3
    # a_er =      100*1e-3
    # b_er =      71.68*1e-3
    # Ri_er =     150*1e-3
    # L_er =      187*1e-3
    # Req_er =    349.0943*1e-3

    # # FCCUROS1.0 1 cell
    # ### DEFINE VARIABLES
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 70*1e-3, 70*1e-3, 25*1e-3, 25*1e-3, 156*1e-3, 120*1e-3, 350.574*1e-3
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 70*1e-3, 70*1e-3, 25*1e-3, 25*1e-3, 156*1e-3, 120*1e-3, 350.574*1e-3
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 70*1e-3, 70*1e-3, 25*1e-3, 25*1e-3, 156*1e-3, 120*1e-3, 350.574*1e-3

    # # FCCUROS1.1 (C1092V2)
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 79.19*1e-3, 70*1e-3, 90*1e-3, 60*1e-3, 150*1e-3, 187*1e-3, 349.0943*1e-3
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 79.19*1e-3, 70*1e-3, 90*1e-3, 60*1e-3, 150*1e-3, 187*1e-3, 349.0943*1e-3
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 79.19*1e-3, 70*1e-3, 90*1e-3, 60*1e-3, 150*1e-3, 187*1e-3, 349.0943*1e-3

    # # C1092V2
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 80.55*1e-3, 69.2*1e-3, 90.41*1e-3, 72.34*1e-3, 150*1e-3, 187*1e-3, 349.0943*1e-3
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 80.55*1e-3, 69.2*1e-3, 90.41*1e-3, 72.34*1e-3, 150*1e-3, 187*1e-3, 349.0943*1e-3
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 80.55*1e-3, 69.2*1e-3, 90.41*1e-3, 72.34*1e-3, 150*1e-3, 187*1e-3, 349.0943*1e-3

    # # LHC
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 104*1e-3, 104*1e-3, 25*1e-3, 25*1e-3, 150*1e-3, 160*1e-3, 344.398*1e-3
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 104*1e-3, 104*1e-3, 25*1e-3, 25*1e-3, 150*1e-3, 160*1e-3, 344.398*1e-3
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 104*1e-3, 104*1e-3, 25*1e-3, 25*1e-3, 150*1e-3, 160*1e-3, 344.398*1e-3

    # # C3794
    A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 73.52 * 1e-3, 131.75 * 1e-3, 106.25 * 1e-3, 118.7 * 1e-3, 150 * 1e-3, 187 * 1e-3, 369.6321578127116 * 1e-3
    A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 73.52 * 1e-3, 131.75 * 1e-3, 106.25 * 1e-3, 118.7 * 1e-3, 150 * 1e-3, 187 * 1e-3, 369.6321578127116 * 1e-3
    A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 73.52 * 1e-3, 131.75 * 1e-3, 106.25 * 1e-3, 118.7 * 1e-3, 150 * 1e-3, 187 * 1e-3, 369.6321578127116 * 1e-3

    # # # C770
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 66.66 * 1e-3, 128.77 * 1e-3, 99.47 * 1e-3, 116.2 * 1e-3, 150 * 1e-3, 187 * 1e-3, 375.2096572945935 * 1e-3
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 66.66 * 1e-3, 128.77 * 1e-3, 99.47 * 1e-3, 116.2 * 1e-3, 150 * 1e-3, 187 * 1e-3, 375.2096572945935 * 1e-3
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 66.66 * 1e-3, 128.77 * 1e-3, 99.47 * 1e-3, 116.2 * 1e-3, 150 * 1e-3, 187 * 1e-3, 375.2096572945935 * 1e-3

    # # # C2183
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 66.96 * 1e-3, 86.8 * 1e-3, 103.28 * 1e-3, 102.36 * 1e-3, 150 * 1e-3, 187 * 1e-3, 365.6418482259177 * 1e-3
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 66.96 * 1e-3, 86.8 * 1e-3, 103.28 * 1e-3, 102.36 * 1e-3, 150 * 1e-3, 187 * 1e-3, 365.6418482259177 * 1e-3
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 66.96 * 1e-3, 86.8 * 1e-3, 103.28 * 1e-3, 102.36 * 1e-3, 150 * 1e-3, 187 * 1e-3, 365.6418482259177 * 1e-3

    # # # C650
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 61.36 * 1e-3, 69.88 * 1e-3, 101.36 * 1e-3, 111.55 * 1e-3, 150 * 1e-3, 187 * 1e-3, 367.3227022198889 * 1e-3
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 61.36 * 1e-3, 69.88 * 1e-3, 101.36 * 1e-3, 111.55 * 1e-3, 150 * 1e-3, 187 * 1e-3, 367.3227022198889 * 1e-3
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 61.36 * 1e-3, 69.88 * 1e-3, 101.36 * 1e-3, 111.55 * 1e-3, 150 * 1e-3, 187 * 1e-3, 367.3227022198889 * 1e-3

    # # # C3345
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 77.58 * 1e-3, 136.64 * 1e-3, 104.57 * 1e-3, 22.93 * 1e-3, 150 * 1e-3, 187 * 1e-3, 352.9841846643084 * 1e-3
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 77.58 * 1e-3, 136.64 * 1e-3, 104.57 * 1e-3, 22.93 * 1e-3, 150 * 1e-3, 187 * 1e-3, 352.9841846643084 * 1e-3
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 77.58 * 1e-3, 136.64 * 1e-3, 104.57 * 1e-3, 22.93 * 1e-3, 150 * 1e-3, 187 * 1e-3, 352.9841846643084 * 1e-3

    # # # C4618
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 60.93 * 1e-3, 112.08 * 1e-3, 119.74 * 1e-3, 28.92 * 1e-3, 150 * 1e-3, 187 * 1e-3, 357.9968262845081 * 1e-3
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 60.93 * 1e-3, 112.08 * 1e-3, 119.74 * 1e-3, 28.92 * 1e-3, 150 * 1e-3, 187 * 1e-3, 357.9968262845081 * 1e-3
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 60.93 * 1e-3, 112.08 * 1e-3, 119.74 * 1e-3, 28.92 * 1e-3, 150 * 1e-3, 187 * 1e-3, 357.9968262845081 * 1e-3
    #
    # # # C4250
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 112.19 * 1e-3, 100.21 * 1e-3, 69.53 * 1e-3, 118.7 * 1e-3, 150 * 1e-3, 187 * 1e-3, 345.8553672307123 * 1e-3
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 112.19 * 1e-3, 100.21 * 1e-3, 69.53 * 1e-3, 118.7 * 1e-3, 150 * 1e-3, 187 * 1e-3, 345.8553672307123 * 1e-3
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 112.19 * 1e-3, 100.21 * 1e-3, 69.53 * 1e-3, 118.7 * 1e-3, 150 * 1e-3, 187 * 1e-3, 345.8553672307123 * 1e-3
    #
    # # # C4123
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 71.43 * 1e-3, 142.78 * 1e-3, 81.48 * 1e-3, 119.27 * 1e-3, 150 * 1e-3, 187 * 1e-3, 375.2660345511048 * 1e-3
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 71.43 * 1e-3, 142.78 * 1e-3, 81.48 * 1e-3, 119.27 * 1e-3, 150 * 1e-3, 187 * 1e-3, 375.2660345511048 * 1e-3
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 71.43 * 1e-3, 142.78 * 1e-3, 81.48 * 1e-3, 119.27 * 1e-3, 150 * 1e-3, 187 * 1e-3, 375.2660345511048 * 1e-3
    #
    # # # TESLA end cell 1
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = 42 * 1e-3, 42 * 1e-3, 12 * 1e-3, 19 * 1e-3, 35 * 1e-3, 57.6524 * 1e-3, 103.353 * 1e-3
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = 40.34 * 1e-3, 40.34 * 1e-3, 10 * 1e-3, 13.5 * 1e-3, 39 * 1e-3, 55.716 * 1e-3, 103.353 * 1e-3
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = 40.34 * 1e-3, 40.34 * 1e-3, 10 * 1e-3, 13.5 * 1e-3, 39 * 1e-3, 55.716 * 1e-3, 103.353 * 1e-3
    #
    # mid40866 = np.array([65.0, 30, 25, 20, 60, 93.5, 160.014])*1e-3
    # mid3794 = np.array([36.76, 65.875, 53.125, 59.35, 75, 93.5, 185.7276])*1e-3
    # midG6 = np.array([51.0510633611297, 50.8492563995682, 37.94024533955164, 27.22177461866947, 73.83844219184923, 93.5, 173.2763])*1e-3
    # midG6 = np.array([54.879053291574124, 53.90857921086279, 33.15859388359441, 15.488719027059194, 81.21661594222239, 93.5, 172.559])*1e-3
    # midG6 = np.array([34.51099669, 54.3707791, 55.51561063, 10.05453999, 76.23570072, 93.5, 176.902])*1e-3
    # midG12 = np.array([62.22222222222222, 66.12612612612612, 30.22022022022022, 23.113113113113116, 71.98698698698699, 93.5, 171.1929])*1e-3
    # endG12 = np.array([62.58258258258258, 57.53753753753754, 17.207207207207208, 12.002002002002001, 80.38038038038039, 93.31191678718535, 171.1929])*1e-3
    #
    # A_m, B_m, a_m, b_m, Ri_m, L_m, Req_m = midG12
    # A_el, B_el, a_el, b_el, Ri_el, L_el, Req_el = midG12
    # A_er, B_er, a_er, b_er, Ri_er, L_er, Req_er = midG12
    # print(A_m)

    n_cell = 3
    step = 2  # step in boundary points in mm
    L_bp_r = 4 * L_m  # 0.0001  #
    L_bp_l = 4 * L_m  # 0.0001  #

    # calculate shift
    shift = (L_bp_r + L_bp_l + (n_cell - 1) * 2 * L_m + L_el + L_er) / 2

    # calculate angles outside loop
    ### CALCULATE x1_el, y1_el, x2_el, y2_el
    data = ([0 + L_bp_l, Ri_el + b_el, L_el + L_bp_l, Req_el - B_el],
            [a_el, b_el, A_el, B_el])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])

    x1el, y1el, x2el, y2el = fsolve(f, np.array(
        [a_el + L_bp_l, Ri_el + 0.5 * b_el, L_el - A_el + L_bp_l, Req_el - 0.5 * B_el]),
                                    args=data,
                                    xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

    ### CALCULATE x1, y1, x2, y2
    data = ([0 + L_bp_l, Ri_m + b_m, L_m + L_bp_l, Req_m - B_m],
            [a_m, b_m, A_m, B_m])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    x1, y1, x2, y2 = fsolve(f, np.array([a_m + L_bp_l, Ri_m + 0.5 * b_m, L_m - A_m + L_bp_l, Req_m - 0.5 * B_m]),
                            args=data, xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

    ### CALCULATE x1_er, y1_er, x2_er, y2_er
    data = ([0 + L_bp_l, Ri_er + b_er, L_er + L_bp_l, Req_er - B_er],
            [a_er, b_er, A_er, B_er])  # data = ([h, k, p, q], [a_m, b_m, A_m, B_m])
    x1er, y1er, x2er, y2er = fsolve(f, np.array(
        [a_er + L_bp_l, Ri_er + 0.5 * b_er, L_er - A_er + L_bp_l, Req_er - 0.5 * B_er]),
                                    args=data,
                                    xtol=1.49012e-12)  # [a_m, b_m-0.3*b_m, L_m-A_m, Req_m-0.7*B_m] initial guess

    with open(r'D:\Dropbox\multipacting\MPGUI21\geodata.n', 'w') as fil:
    # with open(r'D:\Dropbox\CavityDesignHub\C1092V\PostprocessingData\Data\TESLA_End_cell1.txt', 'w') as fil:
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
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

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

        # EQUATOR ARC TO NEXT POINT
        # half of bounding box is required, start is the lower coordinate of the bounding box and end is the upper
        pts = arcTo(L_m + L_bp_l - shift, Req_m - B_m, A_m, B_m, step, [pt[0], pt[1] - B_er],
                    [+ 2 * L_er - x2er + 2 * L_bp_l - shift, Req_er])
        pt = [+ 2 * L_er - x2er + 2 * L_bp_l - shift, y2er]
        for pp in pts:
            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        # STRAIGHT LINE TO NEXT POINT
        lineTo(pt, [+ 2 * L_er - x1er + 2 * L_bp_l - shift, y1er], step)
        pt = [+ 2 * L_er - x1er + 2 * L_bp_l - shift, y1er]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        # ARC
        # half of bounding box is required, start is the lower coordinate of the bounding box and end is the upper
        pts = arcTo(2 * L_er + L_bp_l - shift, Ri_er + b_er, a_er, b_er, step, [pt[0], Ri_er],
                    [L_bp_l + L_el + L_er - shift, y1er])
        pt = [L_bp_l + L_el + L_er - shift, Ri_er]
        for pp in pts:
            fil.write(f"  {pp[1]:.7E}  {pp[0]:.7E}   1.0000000e+00   1.0000000e+00\n")
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        # BEAM PIPE
        lineTo(pt, [2 * n_cell * L_er + L_bp_l + L_bp_r - shift, Ri_er], step)
        pt = [2 * n_cell * L_er + L_bp_l + L_bp_r - shift, Ri_er]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   3.0000000e+00   0.0000000e+00\n")

        # END PATH
        lineTo(pt, [2 * n_cell * L_er + L_bp_l + L_bp_r - shift, 0], step)  # to add beam pipe to right
        pt = [2 * n_cell * L_er + L_bp_l + L_bp_r - shift, 0]
        # lineTo(pt, [2 * n_cell * L_er + L_bp_l - shift, 0], step)
        # pt = [2 * n_cell * L_er + L_bp_l - shift, 0]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

        # CLOSE PATH
        lineTo(pt, start_point, step)
        fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

    plt.show()


def drawCapacitor():
    # l = 18*1e-3
    # d = 12*1e-3
    l = 12 * 1e-3
    d = 18 * 1e-3
    step = 1  # step in boundary points in mm

    shift = l

    with open(r'D:\Dropbox\multipacting\MPGUI21\geodata.n', 'w') as fil:
        fil.write("   2.0000000e-03   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")
        fil.write("   2.5000000e-03   0.0000000e+00   0.0000000e+00   0.0000000e+00\n")  # a point inside the structure
        fil.write("  -3.1415927e+00  -2.7182818e+00   0.0000000e+00   0.0000000e+00\n")  # a point outside the structure
        # SHIFT POINT TO START POINT
        start_point = [-shift, 0]
        fil.write(f"  {start_point[1]:.7E}  {start_point[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        # Add intermediate point
        lineTo(start_point, [-shift, d / 2], step)
        pt = [-shift, d / 2]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

        lineTo(pt, [-shift, d], step)
        pt = [-shift, d]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        # ADD BEAM PIPE LENGTH
        lineTo(pt, [l, d], step)
        pt = [l, d]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

        # Add intermediate point
        lineTo(pt, [l, d / 2], step)
        pt = [l, d / 2]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   1.0000000e+00   1.0000000e+00\n")

        # ADD BEAM PIPE LENGTH
        lineTo(pt, [l, 0], step)
        pt = [l, 0]
        fil.write(f"  {pt[1]:.7E}  {pt[0]:.7E}   0.0000000e+00   0.0000000e+00\n")

        # ADD BEAM PIPE LENGTH
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

    plt.plot(u + a * np.cos(t), v + b * np.sin(t))

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

    plt.plot(inbox[:, 0], inbox[:, 1])

    return inbox


if __name__ == '__main__':
    drawCavity()
    # drawCapacitor()
