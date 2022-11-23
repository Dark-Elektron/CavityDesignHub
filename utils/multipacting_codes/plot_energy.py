# Function program plot_energy(side)
# ---------------------------------------------------------------------------
# Plots the final impact energy.
# INPUT  : 'l'eft or 'r'ight, for a window
# ---------------------------------------------------------------------------
# CALLS TO : clear_window.m, check_inputs.m, check_outputs.m, error_message.m,
#            check_inputs_win.m, check_outputs_win.m, load_output_data.m,
#            load_output_coupler.m, load_output_window.m, plotef.m,
#            check_geotype.m
# 10/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
# 13/12/00 :                - check_geotype.m added
# ---------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
    
def plot_energy(side = None): 
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
                cl = error_message('Unable to plot the impact energy. No initial points.')
                ok = 0
            if ok > 0:
                cl = error_message('Plotting the final impact energy.')
                Window = plt.figure(5)
                clf
                set(Window,'name','MULTIPAC - Output Window I')
                subplot(1,1,1)
                plt.plot(flevel / 1000.0,Efq)
                grid
                hold('on')
                #      plot([min(flevel)/1e3,max(flevel)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
                e0 = interp1(secy1(np.arange(1,e1+1),2),secy1(np.arange(1,e1+1),1),1)
                plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([e0,e0]),'-r')
                plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([secy1(e2,1),secy1(e2,1)]),'-r')
                plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([secy1(e3,1),secy1(e3,1)]),'--r')
                if exist('secy2'):
                    scipy.io.loadmat('secy2')
                    e21 = np.amin(find(secy2(:,2) >= 1))
                    e20 = interp1(secy2(np.arange(1,e21+1),2),secy2(np.arange(1,e21+1),1),1)
                    e22 = np.amax(find(secy2(:,2) >= 1))
                    val,e23 = np.amax(secy2(:,2))
                    plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([e20,e20]),'-b')
                    plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([secy2(e22,1),secy2(e22,1)]),'-b')
                    plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([secy2(e23,1),secy2(e23,1)]),'--b')
                if exist('secy3'):
                    scipy.io.loadmat('secy3')
                    e31 = np.amin(find(secy3(:,2) >= 1))
                    e30 = interp1(secy3(np.arange(1,e31+1),2),secy3(np.arange(1,e31+1),1),1)
                    e32 = np.amax(find(secy3(:,2) >= 1))
                    val,e33 = np.amax(secy3(:,2))
                    plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([e30,e30]),'-g')
                    plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([secy3(e32,1),secy3(e32,1)]),'-g')
                    plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([secy3(e33,1),secy3(e33,1)]),'--g')
                if exist('secy4'):
                    scipy.io.loadmat('secy4')
                    e41 = np.amin(find(secy4(:,2) >= 1))
                    e40 = interp1(secy4(np.arange(1,e41+1),2),secy1(np.arange(1,e41+1),1),1)
                    e42 = np.amax(find(secy4(:,2) >= 1))
                    val,e43 = np.amax(secy4(:,2))
                    plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([e40,e40]),'-k')
                    plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([secy4(e42,1),secy4(e42,1)]),'-k')
                    plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([secy4(e43,1),secy4(e43,1)]),'--k')
                hold('off')
                plt.ylabel(np.array(['Ef_{',num2str(N),'} ']))
                plt.xlabel('Peak Electric Field  [kV/m]')
                plt.title(np.array(['MultiPac 2.1           Final Impact Energy in ev        ',date]))
                htyp = uicontrol('Style','Popup','String','lin|sqrt|loq','Units','Normalized','Position',np.array([0.78,0.85,0.12,0.05]),'Callback','plotef','Tag','Plot type')
                #     hdim = uicontrol('Style','Popup','String','field|voltage|power',...
#	       'Units','Normalized','Position',[0.78 0.78 0.12 0.05],...
#	       'Callback','plotef','Tag','Dimension');
        else:
            if gtype == 2:
                load_output_coupler
                if n == 0:
                    cl = error_message('Unable to plot the impact energy. No initial points.')
                    ok = 0
                if ok > 0:
                    cl = error_message('Plotting the final impact energy.')
                    Window = plt.figure(5)
                    clf
                    set(Window,'name','MULTIPAC - Output Window I')
                    subplot(1,1,1)
                    plt.plot(Pow / 1000.0,Efq)
                    grid
                    hold('on')
                    #      plot([min(Pow)/1e3,max(Pow)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
                    e0 = interp1(secy1(np.arange(1,e1+1),2),secy1(np.arange(1,e1+1),1),1)
                    plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([e0,e0]),'-r')
                    plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy1(e2,1),secy1(e2,1)]),'-r')
                    plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy1(e3,1),secy1(e3,1)]),'--r')
                    if exist('secy2'):
                        scipy.io.loadmat('secy2')
                        e21 = np.amin(find(secy2(:,2) >= 1))
                        e20 = interp1(secy2(np.arange(1,e21+1),2),secy2(np.arange(1,e21+1),1),1)
                        e22 = np.amax(find(secy2(:,2) >= 1))
                        val,e23 = np.amax(secy2(:,2))
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([e20,e20]),'-k')
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy2(e22,1),secy2(e22,1)]),'-k')
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy2(e23,1),secy2(e23,1)]),'--k')
                    if exist('secy3'):
                        scipy.io.loadmat('secy3')
                        e31 = np.amin(find(secy3(:,2) >= 1))
                        e30 = interp1(secy3(np.arange(1,e31+1),2),secy3(np.arange(1,e31+1),1),1)
                        e32 = np.amax(find(secy3(:,2) >= 1))
                        val,e33 = np.amax(secy3(:,2))
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([e30,e30]),'-g')
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy3(e32,1),secy3(e32,1)]),'-g')
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy3(e33,1),secy3(e33,1)]),'--g')
                    if exist('secy4'):
                        scipy.io.loadmat('secy4')
                        e41 = np.amin(find(secy4(:,2) >= 1))
                        e40 = interp1(secy4(np.arange(1,e41+1),2),secy1(np.arange(1,e41+1),1),1)
                        e42 = np.amax(find(secy4(:,2) >= 1))
                        val,e43 = np.amax(secy4(:,2))
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([e40,e40]),'-m')
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy4(e42,1),secy4(e42,1)]),'-m')
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy4(e43,1),secy4(e43,1)]),'--m')
                    hold('off')
                    plt.ylabel(np.array(['Ef_{',num2str(N),'} ']))
                    plt.xlabel('RF power  [kW]')
                    plt.title(np.array(['MultiPac 2.1          Final Impact Energy in eV         ',date]))
                    htyp = uicontrol('Style','Popup','String','lin|sqrt|loq','Units','Normalized','Position',np.array([0.78,0.85,0.12,0.05]),'Callback','plotef','Tag','Plot type')
                    hdim = uicontrol('Style','Popup','String','power|voltage','Units','Normalized','Position',np.array([0.78,0.78,0.12,0.05]),'Callback','plotef','Tag','Dimension')
            else:
                if gtype == 3:
                    load_output_window
                    if ok > 0:
                        if side == 'l':
                            Efq = Eql
                            n = nl
                        else:
                            if side == 'r':
                                Efq = Eqr
                                n = nr
                        if n == 0:
                            cl = error_message(np.array(['Unable to plot the impact energy. No ','initial points.']))
                            ok = 0
                            return
                        if side == 'l':
                            cl = error_message('Plotting the final impact energy. Warm side.')
                        else:
                            if side == 'r':
                                cl = error_message('Plotting the final impact energy. Cold side.')
                            else:
                                cl = error_message(np.array(['To plot the impact energy in window, ','choose warm or cold side.']))
                                return
                        Window = plt.figure(5)
                        clf
                        set(Window,'name','MULTIPAC - Output Window I')
                        subplot(1,1,1)
                        plt.plot(Pow / 1000.0,Efq)
                        grid
                        hold('on')
                        #      plot([min(Pow)/1e3,max(Pow)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
                        e0 = interp1(secy1(np.arange(1,e1+1),2),secy1(np.arange(1,e1+1),1),1)
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([e0,e0]),'-r')
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy1(e2,1),secy1(e2,1)]),'-r')
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy1(e3,1),secy1(e3,1)]),'--r')
                        if exist('secy2'):
                            scipy.io.loadmat('secy2')
                            e21 = np.amin(find(secy2(:,2) >= 1))
                            e20 = interp1(secy2(np.arange(1,e21+1),2),secy2(np.arange(1,e21+1),1),1)
                            e22 = np.amax(find(secy2(:,2) >= 1))
                            val,e23 = np.amax(secy2(:,2))
                            plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([e20,e20]),'-k')
                            plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy2(e22,1),secy2(e22,1)]),'-k')
                            plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy2(e23,1),secy2(e23,1)]),'--k')
                        if exist('secy3'):
                            scipy.io.loadmat('secy3')
                            e31 = np.amin(find(secy3(:,2) >= 1))
                            e30 = interp1(secy3(np.arange(1,e31+1),2),secy3(np.arange(1,e31+1),1),1)
                            e32 = np.amax(find(secy3(:,2) >= 1))
                            val,e33 = np.amax(secy3(:,2))
                            plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([e30,e30]),'-g')
                            plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy3(e32,1),secy3(e32,1)]),'-g')
                            plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy3(e33,1),secy3(e33,1)]),'--g')
                        if exist('secy4'):
                            scipy.io.loadmat('secy4')
                            e41 = np.amin(find(secy4(:,2) >= 1))
                            e40 = interp1(secy4(np.arange(1,e41+1),2),secy1(np.arange(1,e41+1),1),1)
                            e42 = np.amax(find(secy4(:,2) >= 1))
                            val,e43 = np.amax(secy4(:,2))
                            plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([e40,e40]),'-m')
                            plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy4(e42,1),secy4(e42,1)]),'-m')
                            plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy4(e43,1),secy4(e43,1)]),'--m')
                        hold('off')
                        ax = axis
                        plt.axis(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0,0,ax(4)]))
                        plt.ylabel(np.array(['Ef_{',num2str(N),'} ']))
                        plt.xlabel('RF power  [kW]')
                        plt.title(np.array(['MultiPac 2.1           Final Impact Energy in eV        ',date]))
                        htyp = uicontrol('Style','Popup','String','lin|sqrt|loq','Units','Normalized','Position',np.array([0.78,0.85,0.12,0.05]),'Callback','plotef','Tag','Plot type')
                        hdim = uicontrol('Style','Popup','String','power|voltage','Units','Normalized','Position',np.array([0.78,0.78,0.12,0.05]),'Callback','plotef','Tag','Dimension')
    
    if cl == 1:
        clear_window
    
    # --------------------------------------------------------------------