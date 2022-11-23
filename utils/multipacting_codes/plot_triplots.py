# Function program plot_triplots(side,s)
# -----------------------------------------------------------------------
# To plot triplots.
#
# -----------------------------------------------------------------------
# CALLS TO : plot_triplot.m, plot_triplot_coupler.m, plot_triplot_win.m
#            check_geotype.m
# 09/05/00 : Pasi Ylä-Oijala - Rolf Nevanlinna Institute
# 13/12/00 :                 - check_geotype.m added
# -----------------------------------------------------------------------

import numpy as np
    
def plot_triplots(side = None,s = None): 
    if len(varargin) < 2:
        s = 1
    
    ok3 = check_geotype(side)
    if ok3 == 0:
        return
    
    cl = 0
    ok = check_fieldparam
    if ok > 0:
        scipy.io.loadmat('fieldparam')
        gtype = fieldparam(1)
        if gtype == 1:
            plot_triplot(s)
        else:
            if gtype == 2:
                plot_triplot_coupler(1,s)
            else:
                if gtype == 3:
                    if len(varargin) < 1:
                        scipy.io.loadmat('side')
                    if side == 'c':
                        cl = error_message(np.array(['To plot the triplot in window, choose warm side ','or cold side']))
                        return
                    save('side','side')
                    plot_triplot_win(1,side,s)
    
    if cl == 1:
        clear_window
    
    # ---------------------------------------------------------------