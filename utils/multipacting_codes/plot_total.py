# Function program plot_total(side)
# --------------------------------------------------------------------
# Plots the total counter function.
#
# INPUT  : 'l'eft or 'r'ight, for a window
#
# --------------------------------------------------------------------
# CALLS TO : check_inputs.m, check_outputs.m, check_inputs_win.m,
#            check_outputs_win.m, load_output_data.m, error_message.m,
#            load_output_coupler.m, load_output_window.m, plotto.m,
#            check_geotype.m
# 10/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
# 25/09/00 :                - plots the total counter
# 13/12/00 :                - check_geotype.m added
# ---------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
    
def plot_total(side = None): 
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
                cl = error_message('Unable to plot the total counter. No initial points.')
                ok = 0
            if ok > 0:
                cl = error_message('Plotting the total counter function.')
                Window = plt.figure(5)
                clf
                set(Window,'name','MULTIPAC - Output Window I')
                subplot(1,1,1)
                plt.plot(flevel / 1000.0,At)
                grid
                plt.ylabel(np.array(['t_{',num2str(N),'} ']))
                plt.xlabel('Peak Electric Field  [kV/m]')
                plt.title(np.array(['MultiPac 2.1           Total Counter Function         ',date]))
                ax = axis
                plt.axis(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0,ax(3),ax(4)]))
                htyp = uicontrol('Style','Popup','String','lin|rlin|log|rlog','Units','Normalized','Position',np.array([0.78,0.85,0.12,0.05]),'Callback','plotto','Tag','Plot type')
                #      hdim = uicontrol('Style','Popup','String','field|voltage|power',...
#	       'Units','Normalized','Position',[0.78 0.78 0.12 0.05],...
#	       'Callback','plotto','Tag','Dimension');
        else:
            if gtype == 2:
                load_output_coupler
                if n == 0:
                    cl = error_message('Unable to plot the total counter. No initial points.')
                    ok = 0
                if ok > 0:
                    cl = error_message('Plotting the total counter function.')
                    Window = plt.figure(5)
                    clf
                    set(Window,'name','MULTIPAC - Output Window I')
                    subplot(1,1,1)
                    plt.plot(Pow / 1000.0,At)
                    grid
                    plt.ylabel(np.array(['t_{',num2str(N),'} ']))
                    plt.xlabel('RF power  [kW]')
                    plt.title(np.array(['MultiPac 2.1           Total Counter Function         ',date]))
                    ax = axis
                    plt.axis(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0,ax(3),ax(4)]))
                    htyp = uicontrol('Style','Popup','String','lin|rlin|log|rlog','Units','Normalized','Position',np.array([0.78,0.85,0.12,0.05]),'Callback','plotto','Tag','Plot type')
                    hdim = uicontrol('Style','Popup','String','power|voltage','Units','Normalized','Position',np.array([0.78,0.78,0.12,0.05]),'Callback','plotto','Tag','Dimension')
            else:
                if gtype == 3:
                    load_output_window
                    if ok > 0:
                        if side == 'l':
                            A = Al
                            n = nl
                        else:
                            if side == 'r':
                                A = Ar
                                n = nr
                        if n == 0:
                            cl = error_message(np.array(['Unable to plot the total counter. No initial ','points.']))
                            ok = 0
                            return
                        if side == 'l':
                            cl = error_message('Plotting the total counter function. Warm side.')
                        else:
                            if side == 'r':
                                cl = error_message('Plotting the total counter function. Cold side.')
                            else:
                                cl = error_message(np.array(['To plot the total counter function in ','window, choose warm or cold side.']))
                                return
                        Window = plt.figure(5)
                        clf
                        set(Window,'name','MULTIPAC - Output Window I')
                        subplot(1,1,1)
                        plt.plot(Pow / 1000.0,At)
                        grid
                        plt.ylabel(np.array(['t_{',num2str(N),'} ']))
                        plt.xlabel('RF power  [kW]')
                        plt.title(np.array(['MultiPac 2.1           Total Counter Function         ',date]))
                        ax = axis
                        plt.axis(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0,ax(3),ax(4)]))
                        htyp = uicontrol('Style','Popup','String','lin|rlin|log|rlog','Units','Normalized','Position',np.array([0.78,0.85,0.12,0.05]),'Callback','plotto','Tag','Plot type')
                        hdim = uicontrol('Style','Popup','String','power|voltage','Units','Normalized','Position',np.array([0.78,0.78,0.12,0.05]),'Callback','plotto','Tag','Dimension')
    
    if cl == 1:
        clear_window
    
    # --------------------------------------------------------------------