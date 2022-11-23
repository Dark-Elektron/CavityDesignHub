# Function program plot_triplot(ss)
# --------------------------------------------------------------------
# Plots the counter function, the impact energy and the enhanced
# counter function.
#
# --------------------------------------------------------------------
# CALLS TO : load_output_data.m, error_message.m, clear_window.m
# 06/03/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
# ---------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
    
def plot_triplot(ss = None): 
    load_output_data
    cl = 0
    if ok > 0:
        if n == 0:
            error_message('Unable to plot the counters. No initial points.')
            return
        if ok1 * ok2 == 0:
            cl = error_message('Counter functions or impact energy missing.')
        else:
            if ss > 0:
                cl = error_message(np.array(['Plotting the triplot (counter, enhanced ','counter and impact energy).']))
            Window = plt.figure(5)
            clf
            set(Window,'name','MULTIPAC - Output Window I')
            plt.figure(5)
            subplot(3,1,1)
            plt.plot(Efl / 1000000.0,C / n)
            grid
            plt.ylabel(np.array(['c_{',num2str(N),'} / c_0 ']))
            plt.xlabel('Peak electric field  [MV/m]')
            plt.title(np.array(['MultiPac 2.1             Counter function             ',date]))
            ax = axis
            plt.axis(np.array([np.amin(Efl) / 1000000.0,np.amax(Efl) / 1000000.0,0,np.amax(np.array([0.1,ax(4)]))]))
            subplot(3,1,2)
            semilogy(Efl / 1000000.0,Efq)
            grid
            hold('on')
            #    plot([min(Efl)/1e3,max(Efl)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
            e0 = interp1(secy1(np.arange(1,e1+1),2),secy1(np.arange(1,e1+1),1),1)
            plt.plot(np.array([np.amin(Efl) / 1000000.0,np.amax(Efl) / 1000000.0]),np.array([e0,e0]),'-r')
            plt.plot(np.array([np.amin(Efl) / 1000000.0,np.amax(Efl) / 1000000.0]),np.array([secy1(e2,1),secy1(e2,1)]),'-r')
            plt.plot(np.array([np.amin(Efl) / 1000000.0,np.amax(Efl) / 1000000.0]),np.array([secy1(e3,1),secy1(e3,1)]),'--r')
            hold('off')
            plt.ylabel(np.array(['Ef_{',num2str(N),'} ']))
            plt.xlabel('Peak electric field  [MV/m]')
            plt.title('Final Impact Energy in eV')
            ax = axis
            plt.axis(np.array([np.amin(Efl) / 1000000.0,np.amax(Efl) / 1000000.0,0,ax(4)]))
            subplot(3,1,3)
            semilogy(Efl / 1000000.0,(A + 1) / n)
            grid
            plt.xlabel('Voltage   [MV]')
            hold('on')
            plt.plot(np.array([np.amin(Efl) / 1000000.0,np.amax(Efl) / 1000000.0]),np.array([1,1]),'-r')
            hold('off')
            ax = axis
            plt.axis(np.array([np.amin(Efl) / 1000000.0,np.amax(Efl) / 1000000.0,np.amin((A + 1) / n),ax(4)]))
            plt.ylabel(np.array(['e_{',num2str(N),'} / c_0 ']))
            plt.xlabel('Peak electric field   [MV/m]')
            plt.title('Enhanced counter function')
    
    if cl == 1:
        clear_window
    
    # --------------------------------------------------------------------