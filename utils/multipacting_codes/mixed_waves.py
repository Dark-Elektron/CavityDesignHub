# Function program mixed_waves(s)
# ---------------------------------------------------------------------
# Generates mixed waves.
#
# ---------------------------------------------------------------------
# CALLS TO : clear_window.m, check_fieldparam.m, error_message.m, mix.m
# 08/05/00 : Pasi Ylä-Oijala - Rolf Nevalinna Institute
# 17/10/00 :                 - coefficient files removed
# ---------------------------------------------------------------------

import numpy as np
    
def mixed_waves(s = None): 
    if len(varargin) < 1:
        s = 1
    
    cl = 0
    ok1 = check_fieldparam
    if ok1 > 0:
        scipy.io.loadmat('fieldparam')
        gtype = fieldparam(1)
        if gtype == 1:
            cl = error_message('Mixed waves are not defined for a cavity.')
            return
        else:
            R = fieldparam(5) + i * fieldparam(6)
            # mix the coefficient files
            alp1,alp2 = mix(R)
            print('Ref coeff         alpha1              alpha2')
            print(np.array(['  ',num2str(R),'   ',num2str(alp1),'    ',num2str(alp2)]))
            print('                            ')
            scipy.io.loadmat('fields1')
            Er1 = i * Er
            Ez1 = i * Ez
            Bp1 = 4e-07 * pi * H
            scipy.io.loadmat('fields2')
            Er2 = i * Er
            Ez2 = i * Ez
            Bp2 = 4e-07 * pi * H
            r = I2
            z = I1
            Er = alp1 * Er1 + alp2 * Er2
            Ez = alp1 * Ez1 + alp2 * Ez2
            Bp = alp1 * Bp1 + alp2 * Bp2
            save('mixfields','Er','Ez','Bp','rr','zz','r','z')
            if gtype == 2:
                if s > 1:
                    cl = error_message('                                     ')
                    cl = error_message('Generating mixed waves for a coupler.')
                else:
                    if s > 0:
                        cl = error_message('Generating mixed waves for a coupler.')
            else:
                if s > 0:
                    cl = error_message('Generating mixed waves for a window.')
            if s > 1:
                cl = error_message('Generation of mixed waves finished.')
                cl = error_message('                                   ')
    
    if cl == 1:
        clear_window
    
    # ----------------------------------------------------------