# Function program run_field_solver.m
# ----------------------------------------------------------------------------
# Runs the field solver.
#
# ----------------------------------------------------------------------------
# CALLS TO : clear_window.m, test_geometry.m, error_message.m,
#            check_geometry.m, cavity_field.m, coupler_field.m, window_field.m
# 04/04/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
# ----------------------------------------------------------------------------

import matplotlib.pyplot as plt
    
def run_field_solver(s = None): 
    if len(varargin) < 1:
        s = 1
    
    if s == 1:
        clear_window
    
    okt = test_geometry(0)
    if okt == 0:
        return
    
    ok = check_geometry
    if ok > 0:
        plt.figure(1)
        error_message('---------- Field Solver ----------')
        error_message('                                  ')
        scipy.io.loadmat('fieldparam')
        gtype = fieldparam(1)
        if gtype == 1:
            cavity_field(s)
        else:
            if gtype == 2:
                coupler_field(s)
            else:
                if gtype == 3:
                    window_field(s)
        if s > 0:
            #    error_message('                           ');
            error_message('--------- The end ---------')
        error_message('                          ')
    
    # ---------------------------------------------------------------------