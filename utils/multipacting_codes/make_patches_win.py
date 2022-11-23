# Function program [nodes,edges,patches] = make_patches_win(geodata,nodes1,
# edges1,releps)
# -------------------------------------------------------------------------
# Converts geodata-formt into (nodes,edges,pacthes)-format. For a window.
#
# -------------------------------------------------------------------------
# CALLS TO : error_message.m, clear_window.m
# 10/05/00 : Pasi Ylä-Oijala - Rolf Nevanlinna Institute
# 15/06/00 :                 - bugs fixed
# -------------------------------------------------------------------------

import numpy as np
    
def make_patches_win(geodata = None,nodes1 = None,edges1 = None,releps = None): 
    cl = 0
    if releps == 1:
        cl = error_message('Relative pervittyvity of the window is 1?')
    
    n = len(geodata(:,1))
    gr = geodata(np.arange(4,n+1),1)
    gz = geodata(np.arange(4,n+1),2)
    gn = geodata(np.arange(4,n+1),3)
    n = len(gr)
    nodes = []
    edges = []
    patches = []
    nd = len(nodes1)
    tol = np.amin(np.sqrt((nodes1(np.arange(2,nd+1),1) - nodes1(np.arange(1,nd - 1+1),1)) ** 2 + (nodes1(np.arange(2,nd+1),2) - nodes1(np.arange(1,nd - 1+1),2)) ** 2)) / 10000.0
    # add new nodes
    for j in np.arange(1,n+1).reshape(-1):
        dis = np.amin(np.sqrt((nodes1(:,1) - gz(j)) ** 2 + (nodes1(:,2) - gr(j)) ** 2))
        if dis > tol:
            nodes = np.array([[nodes],[np.array([gz,gr])]])
            cl = error_message('Warning: Window has new nodes.')
    
    # find the edges
    n = len(gr)
    nod = np.array([gz,gr])
    edg = np.array([np.transpose((np.arange(1,n - 1+1))),np.transpose(np.array([np.arange(2,n - 1+1),1])),gn(np.arange(1,n - 1+1))])
    # nodes due to the edges in the grid
    z1 = nodes1(edges1(:,1),1)
    r1 = nodes1(edges1(:,1),2)
    z2 = nodes1(edges1(:,2),1)
    r2 = nodes1(edges1(:,2),2)
    zw1 = nod(edg(:,1),1)
    rw1 = nod(edg(:,1),2)
    zw2 = nod(edg(:,2),1)
    rw2 = nod(edg(:,2),2)
    # patch for the window
    new = np.amax(np.amax(edges1)) + 1
    for j in np.arange(1,len(edg)+1).reshape(-1):
        dis1,ind1 = np.amin(np.sqrt((z1 - zw1(j)) ** 2 + (r1 - rw1(j)) ** 2) + np.sqrt((z2 - zw2(j)) ** 2 + (r2 - rw2(j)) ** 2))
        dis2,ind2 = np.amin(np.sqrt((z1 - zw2(j)) ** 2 + (r1 - rw2(j)) ** 2) + np.sqrt((z2 - zw1(j)) ** 2 + (r2 - rw1(j)) ** 2))
        #  [dis1,ind1,dis2,ind2,j]
        dis,jj = np.amin(np.array([dis1,dis2]))
        tmp = np.array([ind1,ind2])
        ind = tmp(jj)
        if dis > tol:
            newedge = edg(j,:)
            newnode = nod(newedge(np.arange(1,2+1)),:)
            dis1,ind1 = np.amin(np.sqrt((nodes1(:,1) - newnode(1,1)) ** 2 + (nodes1(:,2) - newnode(1,2)) ** 2))
            dis2,ind2 = np.amin(np.sqrt((nodes1(:,1) - newnode(2,1)) ** 2 + (nodes1(:,2) - newnode(2,2)) ** 2))
            edges = np.array([[edges],[ind1,ind2,2]])
            patches = np.array([patches,new])
            new = new + 1
        else:
            patches = np.array([patches,ind])
    
    patches = np.array([patches,releps,1,1])
    if cl == 1:
        clear_window
    
    # ------------------------------------------------------------------------
    return nodes,edges,patches