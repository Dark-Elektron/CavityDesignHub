# Function program [nodes,edges,patches] = make_patches(geodata,releps)
# ---------------------------------------------------------------------
# Converts geodata-formt into (nodes,edges,pacthes)-format.
#
# ---------------------------------------------------------------------
# CALLS TO : error_message.m
# 12/04/00 : Pasi Ylä-Oijala - RNI
# ---------------------------------------------------------------------

import numpy as np
    
def make_patches(geodata = None,releps = None): 
    n = len(geodata(:,1)) - 1
    gr = geodata(np.arange(4,n+1),1)
    gz = geodata(np.arange(4,n+1),2)
    gn = geodata(np.arange(4,n+1),3)
    n = len(gr)
    ind0 = find(gn == 0)
    
    if len(ind0) == 0:
        error_message('Can not find the artificial walls of a window?')
    
    # First the nodes of the geometry
    nodes = np.array([gz,gr])
    # Then the edges
    edges = np.array([np.transpose((np.arange(1,n+1))),np.transpose(np.array([np.arange(2,n+1),1])),gn])
    # And finally the pacthes
    ne = len(edges(:,1))
    patches = np.array([np.arange(1,ne+1),releps,1,1])
    # ------------------------------------------------------------------------
    return nodes,edges,patches