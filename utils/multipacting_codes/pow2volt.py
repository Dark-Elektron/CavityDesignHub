# Function program U = pow2volt(P,R,Z)
# --------------------------------------------------------------------
# Transfors rf power to peak voltage for a coaxial coupler.
# INPUT  P : power(s) [W]
#        R : reflection coefficient
#        Z : impedance of the line
# --------------------------------------------------------------------
# CALLS TO : error_message.m
# 10/05/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
# --------------------------------------------------------------------

import numpy as np
    
def pow2volt(P = None,R = None,Z = None): 
    U = np.zeros((P.shape,P.shape))
    if Z > 0:
        U = (1 + np.abs(R)) * np.sqrt(2 * Z * P)
    else:
        error_message('Impendace = 0. Voltage is not defined.')
    
    # --------------------------------------------------------------------
    return U