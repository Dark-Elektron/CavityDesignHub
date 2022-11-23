# Function program ok = check_spacestep(dx)
# ----------------------------------------------------------------
# Checks the value of the space step for creating initial points.
#
# ----------------------------------------------------------------
# CALLS TO : error_message.m, clear_window.m
# 19/05/00 : PYO / RNI
# ----------------------------------------------------------------

    
def check_spacestep(dx = None): 
    cl = 0
    ok = 1
    if dx <= 0:
        cl = error_message('Error: Space step must be > 0.')
        ok = 0
    else:
        if dx < 0.1:
            cl = error_message('Warning: Space step small ( < 0.1 ).')
        else:
            if dx > 10:
                cl = error_message('Warning: Space step large ( > 10 ).')
    
    if len(dx) != 1:
        cl = error_message('Error: Space step must be a positive real number.')
        ok = 0
    
    if cl == 1:
        clear_window
    
    # ----------------------------------------------------------------
    return ok