# Function program mp_cavity_coupler(s)
# --------------------------------------------------------------------------
# Carries out the multipacing analysis in a cavity and a coupler geometry.
#
# --------------------------------------------------------------------------
# CALLS TO : plot_mixed_fields.m, calculate_counters.m, plot_triplot.m,
#            plot_triplot_coupler.m, calculate_distance.m, plot_distance.m,
#            calculate_trajectory.m, plot_trajectory.m, error_message.m
# 09/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
# --------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
    
def mp_cavity_coupler(s = None): 
    if len(varargin) == 0:
        s = 1
    
    scipy.io.loadmat('fieldparam')
    gtype = fieldparam(1)
    # plot the electromagnetic fields
    if gtype == 2:
        plot_mixed_fields(0,'c')
    else:
        plot_FEM_fields(0,'c')
    
    # calculate the counter functions
    calculate_counters(0)
    scipy.io.loadmat('Ccounter')
    scipy.io.loadmat('Acounter')
    scipy.io.loadmat('counter_flevels')
    plt.figure(5)
    if gtype == 1:
        plot_triplot(2)
    else:
        if gtype == 2:
            plot_triplot_coupler(1)
            R = fieldparam(5) + i * fieldparam(6)
            Z = fieldparam(8)
            Pow = volt2pow(flevel,R,Z)
    
    cl = 0
    if np.amax(C) == 0:
        cl = error_message(np.array(['Counter function is dentically zero and multipacting ','analysis is completed.']))
    else:
        val,ind = np.amax(A)
        subplot(3,1,1)
        ax = axis
        hold('on')
        if gtype == 1:
            plt.plot(np.array([flevel(ind) / 1000.0,flevel(ind) / 1000.0]),np.array([ax(3),ax(4)]),'-r')
        else:
            if gtype == 2:
                plt.plot(np.array([Pow(ind) / 1000.0,Pow(ind) / 1000.0]),np.array([ax(3),ax(4)]),'-r')
        hold('off')
        subplot(3,1,2)
        ax = axis
        hold('on')
        if gtype == 1:
            plt.plot(np.array([flevel(ind) / 1000.0,flevel(ind) / 1000.0]),np.array([ax(3),ax(4)]),'-r')
        else:
            if gtype == 2:
                plt.plot(np.array([Pow(ind) / 1000.0,Pow(ind) / 1000.0]),np.array([ax(3),ax(4)]),'-r')
        hold('off')
        subplot(3,1,3)
        ax = axis
        hold('on')
        if gtype == 1:
            plt.plot(np.array([flevel(ind) / 1000.0,flevel(ind) / 1000.0]),np.array([ax(3),ax(4)]),'-r')
        else:
            if gtype == 2:
                plt.plot(np.array([Pow(ind) / 1000.0,Pow(ind) / 1000.0]),np.array([ax(3),ax(4)]),'-r')
        hold('off')
        # calculate the distance map
        calculate_distance(0)
        # plot the distance map
        plot_distance(0,'c',0)
        scipy.io.loadmat('Ddistance')
        val,yi0 = np.amin(np.abs(D))
        if s > 0:
            cl = error_message('                                  ')
        # calculate an electron trajectory
        calculate_trajectory(0,'c')
        # plot the trajectory
        plot_trajectory('c',0)
    
    # ---------------------------------------------------------------------