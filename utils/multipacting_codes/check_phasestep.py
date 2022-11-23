# Function program ok = check_phasestep(dphi)
# ----------------------------------------------------------------
# Checks the value of the phase step for creating initial phases.
#
# ----------------------------------------------------------------
# CALLS TO : error_message.m, clear_window.m
# 19/05/00 : PYO / RNI
# ----------------------------------------------------------------

    
def check_phasestep(dphi = None): 
    cl = 0
    ok = 1
    if dphi <= 0:
        cl = error_message('Error: Phase step must be > 0.')
        ok = 0
    else:
        if dphi < 1:
            cl = error_message('Warning: Phase step rather small ( < 1 ).')
        else:
            if dphi > 10:
                cl = error_message('Warning: Phase step rather large ( > 10 ).')
    
    if len(dphi) != 1:
        cl = error_message('Error: Phase step must be a positive real number.')
        ok = 0
    
    if cl == 1:
        clear_window
    
    # ----------------------------------------------------------------
    return ok