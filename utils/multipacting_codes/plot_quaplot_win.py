# Function program plot_quaplot_win(s,side,ss)
# --------------------------------------------------------------------------
# Plots the counter function, the impact energy, the total counter function
# and the enhanced counter function for a window geometry.
# INPUT  s : 3 - voltage
#            2 - field level
#            1 - power
#        side : 'l'eft or 'r'ight for a window
# -------------------------------------------------------------------------
# CALLS TO : load_output_window.m, error_message.m
# 09/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
# 25/09/00 :                - plot for the total counter added
# -------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
    
def plot_quaplot_win(s = None,side = None,ss = None): 
    if len(varargin) < 1:
        s = 1
    
    if len(varargin) < 2:
        scipy.io.loadmat('side')
    else:
        if len(varargin) == 2:
            save('side','side')
    
    if len(varargin) < 3:
        ss = 1
    
    load_output_window
    if ok > 0:
        if side == 'l':
            C = Cl
            A = Al
            Eq = Eql
            n = nl
            At = Atl
        else:
            if side == 'r':
                C = Cr
                A = Ar
                Eq = Eqr
                n = nr
                At = Atr
        if n == 0:
            error_message('Unable to plot the quaplot. No initial points.')
            return
        if ok1 * ok2 == 0:
            error_message('Counter functions or impact energy missing.')
        else:
            if ss > np.logical_and(0,side) == 'l':
                cl = error_message(np.array(['Plotting the quaplot (counter functions ','and impact energy). Warm side.']))
            else:
                if ss > np.logical_and(0,side) == 'r':
                    cl = error_message(np.array(['Plotting the quaplot (counter functions ','and impact energy). Cold side.']))
            Window = plt.figure(5)
            clf
            set(Window,'name','MULTIPAC - Output Window I')
            plt.figure(5)
            if s == 3:
                subplot(4,1,1)
                plt.plot(U / 1000.0,C / n)
                grid
                plt.ylabel(np.array(['c_{',num2str(N),'} / c_0 ']))
                plt.xlabel('Voltage  [kV]')
                if side == 'l':
                    plt.title(np.array(['MultiPac 2.1       Counter Function   Warm Side       ',date]))
                else:
                    if side == 'r':
                        plt.title(np.array(['MultiPac 2.1       Counter Function   Cold Side       ',date]))
                ax = axis
                plt.axis(np.array([np.amin(U) / 1000.0,np.amax(U) / 1000.0,0,np.amax(np.array([0.1,ax(4)]))]))
                subplot(4,1,2)
                plt.plot(U / 1000.0,Eq)
                grid
                hold('on')
                #      plot([min(flevel)/1e3,max(flevel)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
                e0 = interp1(secy1(np.arange(1,e1+1),2),secy1(np.arange(1,e1+1),1),1)
                plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([e0,e0]),'-r')
                plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([secy1(e2,1),secy1(e2,1)]),'-r')
                plt.plot(np.array([np.amin(flevel) / 1000.0,np.amax(flevel) / 1000.0]),np.array([secy1(e3,1),secy1(e3,1)]),'--r')
                hold('off')
                plt.ylabel(np.array(['Ef_{',num2str(N),'} ']))
                plt.xlabel('Voltage [kV]')
                plt.title('Final Impact Energy in eV')
                ax = axis
                plt.axis(np.array([np.amin(U) / 1000.0,np.amax(U) / 1000.0,0,ax(4)]))
                subplot(4,1,3)
                semilogy(U / 1000.0,(A + 1) / n)
                grid
                plt.xlabel('Voltage   [kV]')
                hold('on')
                plt.plot(np.array([np.amin(U) / 1000.0,np.amax(U) / 1000.0]),np.array([1,1]),'-r')
                hold('off')
                ax = axis
                plt.axis(np.array([np.amin(U) / 1000.0,np.amax(U) / 1000.0,np.amin((A + 1) / n),ax(4)]))
                plt.ylabel(np.array(['e_{',num2str(N),'} / c_0 ']))
                plt.xlabel('Voltage  [kV]')
                plt.title('Enhanced Counter Function')
                subplot(4,1,4)
                semilogy(U / 1000.0,(At + 1) / n)
                grid
                plt.xlabel('Voltage   [kV]')
                hold('on')
                plt.plot(np.array([np.amin(U) / 1000.0,np.amax(U) / 1000.0]),np.array([1,1]),'-r')
                hold('off')
                ax = axis
                plt.axis(np.array([np.amin(U) / 1000.0,np.amax(U) / 1000.0,np.amin((At + 1) / n),ax(4)]))
                plt.ylabel(np.array(['t_{',num2str(N),'} / c_0 ']))
                plt.xlabel('Voltage  [kV]')
                plt.title('Total Counter Function')
            else:
                if s == 2:
                    subplot(4,1,1)
                    plt.plot(Efl / 1000.0,C / n)
                    grid
                    plt.ylabel(np.array(['c_{',num2str(N),'} / c_0 ']))
                    plt.xlabel('Peak electric field  [kV/m]')
                    if side == 'l':
                        plt.title(np.array(['MultiPac 2.1       Counter Function   Warm Side       ',date]))
                    else:
                        if side == 'r':
                            plt.title(np.array(['MultiPac 2.1       Counter Function   Cold Side       ',date]))
                    ax = axis
                    plt.axis(np.array([np.amin(Efl) / 1000.0,np.amax(Efl) / 1000.0,0,np.amax(np.array([0.1,ax(4)]))]))
                    subplot(4,1,2)
                    plt.plot(Efl / 1000.0,Eq)
                    grid
                    hold('on')
                    #      plot([min(Efl)/1e3,max(Efl)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
                    e0 = interp1(secy1(np.arange(1,e1+1),2),secy1(np.arange(1,e1+1),1),1)
                    plt.plot(np.array([np.amin(Efl) / 1000.0,np.amax(Efl) / 1000.0]),np.array([e0,e0]),'-r')
                    plt.plot(np.array([np.amin(Efl) / 1000.0,np.amax(Efl) / 1000.0]),np.array([secy1(e2,1),secy1(e2,1)]),'-r')
                    plt.plot(np.array([np.amin(Efl) / 1000.0,np.amax(Efl) / 1000.0]),np.array([secy1(e3,1),secy1(e3,1)]),'--r')
                    hold('off')
                    plt.ylabel(np.array(['Ef_{',num2str(N),'} ']))
                    plt.xlabel('Peak electric field  [kV/m]')
                    plt.title('Final Impact Energy in eV')
                    ax = axis
                    plt.axis(np.array([np.amin(Efl) / 1000.0,np.amax(Efl) / 1000.0,0,ax(4)]))
                    subplot(4,1,3)
                    semilogy(Efl / 1000.0,(A + 1) / n)
                    grid
                    plt.xlabel('Voltage   [kV]')
                    hold('on')
                    plt.plot(np.array([np.amin(Efl) / 1000.0,np.amax(Efl) / 1000.0]),np.array([1,1]),'-r')
                    hold('off')
                    ax = axis
                    plt.axis(np.array([np.amin(Efl) / 1000.0,np.amax(Efl) / 1000.0,np.amin((A + 1) / n),ax(4)]))
                    plt.ylabel(np.array(['e_{',num2str(N),'} / c_0 ']))
                    plt.xlabel('Peak electric field   [kV/m]')
                    plt.title('Enhanced counter function')
                    subplot(4,1,4)
                    semilogy(Efl / 1000.0,(At + 1) / n)
                    grid
                    plt.xlabel('Voltage   [kV]')
                    hold('on')
                    plt.plot(np.array([np.amin(Efl) / 1000.0,np.amax(Efl) / 1000.0]),np.array([1,1]),'-r')
                    hold('off')
                    ax = axis
                    plt.axis(np.array([np.amin(Efl) / 1000.0,np.amax(Efl) / 1000.0,np.amin((At + 1) / n),ax(4)]))
                    plt.ylabel(np.array(['t_{',num2str(N),'} / c_0 ']))
                    plt.xlabel('Peak electric field   [kV/m]')
                    plt.title('Total counter function')
                else:
                    if s == 1:
                        subplot(4,1,1)
                        plt.plot(Pow / 1000.0,C / n)
                        grid
                        plt.ylabel(np.array(['c_{',num2str(N),'} / c_0 ']))
                        plt.xlabel('Power  [kW]')
                        if side == 'l':
                            plt.title(np.array(['MultiPac 2.1       Counter Function   Warm Side       ',date]))
                        else:
                            if side == 'r':
                                plt.title(np.array(['MultiPac 2.1       Counter Function   Cold Side       ',date]))
                        ax = axis
                        plt.axis(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0,0,np.amax(np.array([0.1,ax(4)]))]))
                        subplot(4,1,2)
                        plt.plot(Pow / 1000.0,Eq)
                        grid
                        hold('on')
                        #      plot([min(Pow)/1e3,max(Pow)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
                        e0 = interp1(secy1(np.arange(1,e1+1),2),secy1(np.arange(1,e1+1),1),1)
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([e0,e0]),'-r')
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy1(e2,1),secy1(e2,1)]),'-r')
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([secy1(e3,1),secy1(e3,1)]),'--r')
                        hold('off')
                        plt.ylabel(np.array(['Ef_{',num2str(N),'} ']))
                        plt.xlabel('Power [kW]')
                        plt.title('Final Impact Energy in eV')
                        ax = axis
                        plt.axis(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0,0,ax(4)]))
                        subplot(4,1,3)
                        semilogy(Pow / 1000.0,(A + 1) / n)
                        grid
                        plt.xlabel('Voltage   [kV]')
                        hold('on')
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([1,1]),'-r')
                        hold('off')
                        ax = axis
                        plt.axis(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0,np.amin((A + 1) / n),ax(4)]))
                        plt.ylabel(np.array(['e_{',num2str(N),'} / c_0 ']))
                        plt.xlabel('Power  [kW]')
                        plt.title('Enhanced counter function')
                        subplot(4,1,4)
                        semilogy(Pow / 1000.0,(At + 1) / n)
                        grid
                        plt.xlabel('Voltage   [kV]')
                        hold('on')
                        plt.plot(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0]),np.array([1,1]),'-r')
                        hold('off')
                        ax = axis
                        plt.axis(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0,np.amin((At + 1) / n),ax(4)]))
                        plt.ylabel(np.array(['t_{',num2str(N),'} / c_0 ']))
                        plt.xlabel('Power  [kW]')
                        plt.title('Total counter function')
    
    # --------------------------------------------------------------------