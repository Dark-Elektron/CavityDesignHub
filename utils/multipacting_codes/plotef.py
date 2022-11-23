# Function program plotef.m
# --------------------------------------------------------------------------
# Plots the final impact energy.
#
# --------------------------------------------------------------------------
# CALLS TO : load_output_data.m, load_output_coupler.m, load_output_window.m
# 06/03/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
# --------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
    
def plotef(side = None): 
    if len(varargin) < 1:
        scipy.io.loadmat('side')
    
    scipy.io.loadmat('fieldparam')
    gtype = fieldparam(1)
    if gtype == 1:
        load_output_data
    else:
        if gtype == 2:
            load_output_coupler
        else:
            if gtype == 3:
                load_output_window
                if side == 'l':
                    Efq = Eql
                    n = nl
                else:
                    if side == 'r':
                        Efq = Eqr
                        n = nr
    
    htyp = findobj('Tag','Plot type')
    if gtype == 1:
        dval = 3
    else:
        hdim = findobj('Tag','Dimension')
        dval = get(hdim,'value')
    
    tval = get(htyp,'value')
    if tval == 1:
        if dval == 2:
            plt.plot(U / 1000.0,Efq)
            grid
            plt.xlabel('Peak voltage   [kV]')
            ax = axis
            plt.axis(np.array([np.amin(U) / 1000.0,np.amax(U) / 1000.0,0,ax(4)]))
        else:
            if dval == 3:
                plt.plot(Efl / 1000.0,Efq)
                grid
                plt.xlabel('Peak Electric Field   [kV/m]')
                ax = axis
                plt.axis(np.array([np.amin(Efl) / 1000.0,np.amax(Efl) / 1000.0,0,ax(4)]))
            else:
                if dval == 1:
                    plt.plot(Pow / 1000.0,Efq)
                    grid
                    plt.xlabel('Forward Power   [kW]')
                    ax = axis
                    plt.axis(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0,0,ax(4)]))
        plt.ylabel(np.array(['Ef_{',num2str(N),'} ']))
        plt.title(np.array(['MultiPac 2.0            Final Impact Energy in ev          ',date]))
    else:
        if tval == 2:
            if dval == 2:
                plt.plot(U / 1000.0,np.sqrt(Efq))
                grid
                plt.xlabel('Peak voltage   [kV]')
                ax = axis
                plt.axis(np.array([np.amin(U) / 1000.0,np.amax(U) / 1000.0,0,ax(4)]))
            else:
                if dval == 3:
                    plt.plot(Efl / 1000.0,np.sqrt(Efq))
                    grid
                    plt.xlabel('Peak Electric Field   [kV/m]')
                    ax = axis
                    plt.axis(np.array([np.amin(Efl) / 1000.0,np.amax(Efl) / 1000.0,0,ax(4)]))
                else:
                    if dval == 1:
                        plt.plot(Pow / 1000.0,np.sqrt(Efq))
                        grid
                        plt.xlabel('Forward Power   [kW]')
                        ax = axis
                        plt.axis(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0,0,ax(4)]))
            plt.ylabel(np.array(['sqrt of Ef_{',num2str(N),'} ']))
            plt.title(np.array(['MultiPac 2.0        Sqrt of Final Impact Energy in ev^{1/2}       ',date]))
        else:
            if tval == 3:
                if dval == 2:
                    semilogy(U / 1000.0,Efq + 1)
                    grid
                    plt.xlabel('Peak voltage   [kV]')
                    ax = axis
                    plt.axis(np.array([np.amin(U) / 1000.0,np.amax(U) / 1000.0,0,ax(4)]))
                else:
                    if dval == 3:
                        semilogy(Efl / 1000.0,Efq + 1)
                        grid
                        plt.xlabel('Peak Electric Field   [kV/m]')
                        ax = axis
                        plt.axis(np.array([np.amin(Efl) / 1000.0,np.amax(Efl) / 1000.0,0,ax(4)]))
                    else:
                        if dval == 1:
                            semilogy(Pow / 1000.0,Efq + 1)
                            grid
                            plt.xlabel('Forward Power   [kW]')
                            ax = axis
                            plt.axis(np.array([np.amin(Pow) / 1000.0,np.amax(Pow) / 1000.0,0,ax(4)]))
                plt.ylabel(np.array(['Ef_{',num2str(N),'} ']))
                plt.title(np.array(['MultiPac 2.0            Final Impact Energy [ev] in log_{10}-scale         ',date]))
    
    if tval == np.logical_or(1,tval) == 3:
        #  sec1 = secy1(e1,1);
        sec1 = interp1(secy1(np.arange(1,e1+1),2),secy1(np.arange(1,e1+1),1),1)
        sec2 = secy1(e2,1)
        sec3 = secy1(e3,1)
    else:
        if tval == 2:
            #  sec1 = sqrt(secy1(e1,1));
            sec1 = np.sqrt(interp1(secy1(np.arange(1,e1+1),2),secy1(np.arange(1,e1+1),1),1))
            sec2 = np.sqrt(secy1(e2,1))
            sec3 = np.sqrt(secy1(e3,1))
    
    if dval == 2:
        fl1 = np.amin(U) / 1000.0
        fl2 = np.amax(U) / 1000.0
    else:
        if dval == 3:
            fl1 = np.amin(Efl) / 1000.0
            fl2 = np.amax(Efl) / 1000.0
        else:
            if dval == 1:
                fl1 = np.amin(Pow) / 1000.0
                fl2 = np.amax(Pow) / 1000.0
    
    # plot the secondaty yield tresholds and the maxium yield
    hold('on')
    plt.plot(np.array([fl1,fl2]),np.array([sec1,sec1]),'-r')
    plt.plot(np.array([fl1,fl2]),np.array([sec2,sec2]),'-r')
    plt.plot(np.array([fl1,fl2]),np.array([sec3,sec3]),'--r')
    hold('off')
    # plot the other secondary yield curves
    if exist('secy2'):
        scipy.io.loadmat('secy2')
        e1 = np.amin(find(secy2(:,2) >= 1))
        e2 = np.amax(find(secy2(:,2) >= 1))
        val,e3 = np.amax(secy2(:,2))
        if tval == np.logical_or(1,tval) == 3:
            sec1 = interp1(secy2(np.arange(1,e1+1),2),secy2(np.arange(1,e1+1),1),1)
            sec2 = secy2(e2,1)
            sec3 = secy2(e3,1)
        else:
            if tval == 2:
                sec1 = np.sqrt(interp1(secy2(np.arange(1,e1+1),2),secy2(np.arange(1,e1+1),1),1))
                sec2 = np.sqrt(secy2(e2,1))
                sec3 = np.sqrt(secy2(e3,1))
        # plot the secondaty yield tresholds and the maxium yield
        hold('on')
        plt.plot(np.array([fl1,fl2]),np.array([sec1,sec1]),'-k')
        plt.plot(np.array([fl1,fl2]),np.array([sec2,sec2]),'-k')
        plt.plot(np.array([fl1,fl2]),np.array([sec3,sec3]),'--k')
        hold('off')
    
    if exist('secy3'):
        scipy.io.loadmat('secy3')
        e1 = np.amin(find(secy3(:,2) >= 1))
        e2 = np.amax(find(secy3(:,2) >= 1))
        val,e3 = np.amax(secy3(:,2))
        if tval == np.logical_or(1,tval) == 3:
            sec1 = interp1(secy3(np.arange(1,e1+1),2),secy2(np.arange(1,e1+1),1),1)
            sec2 = secy3(e2,1)
            sec3 = secy3(e3,1)
        else:
            if tval == 2:
                sec1 = np.sqrt(interp1(secy3(np.arange(1,e1+1),2),secy2(np.arange(1,e1+1),1),1))
                sec2 = np.sqrt(secy3(e2,1))
                sec3 = np.sqrt(secy3(e3,1))
        # plot the secondaty yield tresholds and the maxium yield
        hold('on')
        plt.plot(np.array([fl1,fl2]),np.array([sec1,sec1]),'-g')
        plt.plot(np.array([fl1,fl2]),np.array([sec2,sec2]),'-g')
        plt.plot(np.array([fl1,fl2]),np.array([sec3,sec3]),'--g')
        hold('off')
    
    if exist('secy4'):
        scipy.io.loadmat('secy4')
        e1 = np.amin(find(secy4(:,2) >= 1))
        e2 = np.amax(find(secy4(:,2) >= 1))
        val,e3 = np.amax(secy4(:,2))
        if tval == np.logical_or(1,tval) == 3:
            sec1 = interp1(secy4(np.arange(1,e1+1),2),secy2(np.arange(1,e1+1),1),1)
            sec2 = secy4(e2,1)
            sec3 = secy4(e3,1)
        else:
            if tval == 2:
                sec1 = np.sqrt(interp1(secy4(np.arange(1,e1+1),2),secy2(np.arange(1,e1+1),1),1))
                sec2 = np.sqrt(secy4(e2,1))
                sec3 = np.sqrt(secy4(e3,1))
        # plot the secondaty yield tresholds and the maxium yield
        hold('on')
        plt.plot(np.array([fl1,fl2]),np.array([sec1,sec1]),'-m')
        plt.plot(np.array([fl1,fl2]),np.array([sec2,sec2]),'-m')
        plt.plot(np.array([fl1,fl2]),np.array([sec3,sec3]),'--m')
        hold('off')
    
    # --------------------------------------------------------------------