# Function program mp_window(s)
# -------------------------------------------------------------------------
# Carries out the multipacing analysis in a window geometry.
#
# -------------------------------------------------------------------------
# CALLS TO : plot_mixed_fields.m, error_message.m, calculate_counters.m,
#            plot_triplot_win.m, calculate_distance_win.m, plot_distance.m,
#            calculate_trajectory.m, plot_trajectory.m
# 14/03/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
# -------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
    
def mp_window(s = None): 
    if len(varargin) < 1:
        s = 0
    
    # plot the electromagnetic fields
    plot_mixed_fields(0,'l')
    # calculate the counter functions (for both sides)
    calculate_counters_win(0)
    if s > 0:
        error_message('Triplot for the warm side.')
    
    plot_triplot_win(1,'l')
    scipy.io.loadmat('Acounterl')
    scipy.io.loadmat('Ccounterl')
    if np.amax(C) == 0:
        error_message('Counter function on the warm side is identically zero.')
    else:
        scipy.io.loadmat('counter_flevels')
        val,ind = np.amax(A)
        scipy.io.loadmat('fieldparam')
        Pow = volt2pow(flevel,fieldparam(5) + i * fieldparam(6),fieldparam(8))
        Window = plt.figure(5)
        set(Window,'name','MULTIPAC - Output Window I')
        plt.figure(5)
        subplot(3,1,1)
        ax = axis
        hold('on')
        plt.plot(np.array([Pow(ind) / 1000.0,Pow(ind) / 1000.0]),np.array([ax(3),ax(4)]),'-r')
        hold('off')
        subplot(3,1,2)
        ax = axis
        hold('on')
        plt.plot(np.array([Pow(ind) / 1000.0,Pow(ind) / 1000.0]),np.array([ax(3),ax(4)]),'-r')
        hold('off')
        subplot(3,1,3)
        ax = axis
        hold('on')
        plt.plot(np.array([Pow(ind) / 1000.0,Pow(ind) / 1000.0]),np.array([ax(3),ax(4)]),'-r')
        hold('off')
        # calculate the distance map
        calculate_distance_win(0,'l')
        # plot the distance map
        plot_distance(0,'l',0)
        scipy.io.loadmat('Ddistancel')
        val,yi0 = np.amin(np.abs(D))
        if s > 0:
            error_message('                                  ')
        # calculate an electron trajectory
        calculate_trajectory(0,'l')
        # plot the trajectory
        plot_trajectory('l',0)
    
    if s > 0:
        error_message('                                    ')
    
    error_message('To see the cold side, press any key.')
    pause
    plot_triplot_win(1,'r')
    if s > 0:
        error_message('Triplot for the cold side.')
    
    scipy.io.loadmat('Acounterr')
    scipy.io.loadmat('Ccounterr')
    if np.amax(C) == 0:
        error_message('Counter function on the cold side is identically zero.')
    else:
        scipy.io.loadmat('counter_flevels')
        val,ind = np.amax(A)
        scipy.io.loadmat('fieldparam')
        Pow = volt2pow(flevel,fieldparam(5) + i * fieldparam(6),fieldparam(8))
        Window = plt.figure(5)
        set(Window,'name','MULTIPAC - Output Window I')
        plt.figure(5)
        subplot(3,1,1)
        ax = axis
        hold('on')
        plt.plot(np.array([Pow(ind) / 1000.0,Pow(ind) / 1000.0]),np.array([ax(3),ax(4)]),'-r')
        hold('off')
        subplot(3,1,2)
        ax = axis
        hold('on')
        plt.plot(np.array([Pow(ind) / 1000.0,Pow(ind) / 1000.0]),np.array([ax(3),ax(4)]),'-r')
        hold('off')
        subplot(3,1,3)
        ax = axis
        hold('on')
        plt.plot(np.array([Pow(ind) / 1000.0,Pow(ind) / 1000.0]),np.array([ax(3),ax(4)]),'-r')
        hold('off')
        # calculate the distance map
        calculate_distance_win(0,'r')
        # plot the distance map
        plot_distance(0,'r',0)
        scipy.io.loadmat('Ddistancer')
        val,yi0 = np.amin(np.abs(D))
        if s > 0:
            error_message('                                  ')
        # calculate an electron trajectory
        calculate_trajectory(0,'r')
        # plot the trajectory
        plot_trajectory('r',0)
    
    # ---------------------------------------------------------------------