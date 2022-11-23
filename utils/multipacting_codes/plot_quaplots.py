# Function program plot_quaplots(side,s)
# -----------------------------------------------------------------------
# To plot quaplots.
#
# -----------------------------------------------------------------------
# CALLS TO : plot_quaplot.m, plot_quaplot_coupler.m, plot_quaplot_win.m
#            check_geotype.m
# 25/09/00 : Pasi Ylä-Oijala - Rolf Nevanlinna Institute
# 13/12/00 :                 - check_geotype.m added
# -----------------------------------------------------------------------

import numpy as np
    
def plot_quaplots(side = None,s = None): 
    if len(varargin) < 2:
        s = 1
    
    ok3 = check_geotype(side)
    if ok3 == 0:
        return
    
    ok = check_fieldparam
    if ok > 0:
        scipy.io.loadmat('fieldparam')
        gtype = fieldparam(1)
        if gtype == 1:
            plot_quaplot(s)
        else:
            if gtype == 2:
                plot_quaplot_coupler(1,s)
            else:
                if gtype == 3:
                    if len(varargin) < 1:
                        scipy.io.loadmat('side')
                    if side == 'c':
                        cl = error_message(np.array(['To plot the quaplot in window, choose warm side ','or cold side']))
                        return
                    save('side','side')
                    plot_quaplot_win(1,side,s)
    
    # ---------------------------------------------------------------