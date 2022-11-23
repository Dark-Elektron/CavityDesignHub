# Function program plot_counter(side)
# --------------------------------------------------------------------
# Plots the counter function.
# INPUT  side : 'l'eft or 'r'ight, for a window, 'c' otherwise
#
# --------------------------------------------------------------------
# CALLS TO : check_inputs.m, check_outputs.m, check_inputs_win.m,
#            check_outputs_win.m, load_output_data.m, error_message.m,
#            load_output_coupler.m, load_output_window.m, plotcn.m
#            check_geotype.m
# 10/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
# 13/12/00 :                 - check_geotype.m added
# ---------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
    
def plot_counter(side = None): 
    if len(varargin) < 1:
        if exist('side'):
            scipy.io.loadmat('side')
        else:
            side = 'c'
    
    save('side','side')
    ok3 = check_geotype(side)
    if ok3 == 0:
        return
    
    if side == np.logical_or('l',side) == 'r':
        ok1 = check_inputs_win
        ok2 = check_outputs_win
    else:
        ok1 = check_inputs
        ok2 = check_outputs
    
    cl = 0
    if ok1 * ok2 > 0:
        scipy.io.loadmat('fieldparam')
        gtype = fieldparam(1)
        if gtype == 1:
            load_output_data
            if n == 0:
                cl = error_message('Unable to plot the counters. No initial points.')
                ok = 0
            if ok > 0:
                cl = error_message('Plotting the counter function.')
                Window = plt.figure(5)
                clf
                set(Window,'name','MULTIPAC - Output Window I')
                subplot(1,1,1)
                plt.plot(flevel / 1000.0,C)
                grid
                plt.ylabel(np.array(['c_{',num2str(N),'} ']))
                plt.xlabel('Peak Electric Field  [kV/m]')
                plt.title(np.array(['MultiPac 2.1             Counter Function               ',date]))
                ax = axis
                plt.axis(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0,0,ax(4)]))
                htyp = uicontrol('Style','Popup','String','lin|rlin|rlin%','Units','Normalized','Position',np.array([0.78,0.85,0.12,0.05]),'Callback','plotcn','Tag','Plot type')
                #      hdim = uicontrol('Style','Popup','String','power|field|voltage',...
#	        'Units','Normalized','Position',[0.78 0.78 0.12 0.05],...
#	        'Callback','plotcn','Tag','Dimension');
        else:
            if gtype == 2:
                load_output_coupler
                if n == 0:
                    cl = error_message('Unable to plot the counters. No initial points.')
                    ok = 0
                if ok > 0:
                    cl = error_message('Plotting the counter function.')
                    Window = plt.figure(5)
                    clf
                    set(Window,'name','MULTIPAC - Output Window I')
                    subplot(1,1,1)
                    plt.plot(Pow / 1000.0,C)
                    grid
                    plt.ylabel(np.array(['c_{',num2str(N),'} ']))
                    plt.xlabel('RF Power  [kW]')
                    plt.title(np.array(['MultiPac 2.1            Counter Function                ',date]))
                    ax = axis
                    plt.axis(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0,0,ax(4)]))
                    htyp = uicontrol('Style','Popup','String','lin|rlin|rlin%','Units','Normalized','Position',np.array([0.78,0.85,0.12,0.05]),'Callback','plotcn','Tag','Plot type')
                    hdim = uicontrol('Style','Popup','String','power|voltage','Units','Normalized','Position',np.array([0.78,0.78,0.12,0.05]),'Callback','plotcn','Tag','Dimension')
            else:
                if gtype == 3:
                    load_output_window
                    if ok > 0:
                        if side == 'l':
                            C = Cl
                            n = nl
                        else:
                            if side == 'r':
                                C = Cr
                                n = nr
                        if n == 0:
                            cl = error_message('Unable to plot the counters. No initial points.')
                            ok = 0
                            return
                        if side == 'l':
                            cl = error_message('Plotting the counter function. Warm side.')
                        else:
                            if side == 'r':
                                cl = error_message('Plotting the counter function. Cold side.')
                            else:
                                cl = error_message(np.array(['To plot the counter function in window, ','choose Warm Side or Cold Side.']))
                                return
                        Window = plt.figure(5)
                        clf
                        set(Window,'name','MULTIPAC - Output Window I')
                        subplot(1,1,1)
                        plt.plot(Pow / 1000.0,C)
                        grid
                        plt.ylabel(np.array(['c_{',num2str(N),'} ']))
                        plt.xlabel('RF Power  [kW]')
                        plt.title(np.array(['MultiPac 2.1            Counter Function                ',date]))
                        ax = axis
                        plt.axis(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0,0,ax(4)]))
                        htyp = uicontrol('Style','Popup','String','lin|rlin|rlin%','Units','Normalized','Position',np.array([0.78,0.85,0.12,0.05]),'Callback','plotcn','Tag','Plot type')
                        hdim = uicontrol('Style','Popup','String','power|voltage','Units','Normalized','Position',np.array([0.78,0.78,0.12,0.05]),'Callback','plotcn','Tag','Dimension')
    
    if cl == 1:
        clear_window
    
    # --------------------------------------------------------------------