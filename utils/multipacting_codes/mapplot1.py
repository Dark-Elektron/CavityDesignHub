# Function program iy0 = mapplot1(map,q,par,y00,bo,trajec)
# -------------------------------------------------------------------------
# Pplots D_N, Ea_N or Ef_N mapping.
# INPUT  map : (1,n) vector giving the place-time distance between the
#              place of initial emission and final place of electron after
#              N impacts, n is the number of initial electrons
#        q   : type of parameter map, q == 1, map = D_N, q == 2,
#              map = Ea_N and q == 3, map = Ef_N
#        par : (7,1) vector listing the MP parameters
#        y00 : initial sites
#        bo  : the boundary curve
# OUTPUT iy0 : indexes of the chosen initial electrons referring to matrix
#              y00
# -------------------------------------------------------------------------
# CALLS TO : none
# 01/02/99 : Marko Ukkola - Rolf Nevanlinna Institute (mapplot.m)
# 07/03/00 : Pasi Ylä-Oijala - modifications (mapplot1.m)
# -------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
    
def mapplot1(D_N = None,q = None,par = None,y00 = None,bo = None,trajec = None): 
    N = par(5)
    if q != 1:
        if q == 2:
            ti = np.array(['Ea_'])
        else:
            ti = np.array(['Ef_'])
        q = 1.6e-19
    else:
        ti = np.array([' Distance map   d_{',num2str(N),'}'])
    
    k = 0
    Prz = np.zeros((y00.shape,y00.shape))
    Prz[1,:] = y00(1,:)
    ky = 1
    yt = np.zeros((1,y00.shape[1-1]))
    yt[1] = y00(1,4)
    for i in np.arange(1,y00.shape[1-1]+1).reshape(-1):
        if y00(i,4) == Prz(1,4):
            k = k + 1
            Prz[k,:] = y00(i,:)
        if y00(i,4) < yt(1):
            yt[np.arange[2,[ky + 1]+1]] = yt(np.arange(1,ky+1))
            yt[1] = y00(i,4)
            ky = ky + 1
        else:
            if y00(i,4) > yt(ky):
                yt[ky + 1] = y00(i,4)
                ky = ky + 1
            else:
                if (np.logical_and((y00(i,4) > np.logical_and(yt(1),y00(i,4)) < yt(ky)),np.all(y00(i,4) != yt))):
                    pla = sum(yt(np.arange(1,ky+1)) < y00(i,4))
                    yt[np.arange[[pla + 2],[ky + 1]+1]] = yt(np.arange((pla + 1),(ky)+1))
                    yt[pla + 1] = y00(i,4)
                    ky = ky + 1
    
    yt = yt(np.arange(1,ky+1))
    Prz = Prz(np.arange(1,k+1),:)
    map = np.zeros((len(yt),Prz.shape[1-1]))
    for i in np.arange(1,y00.shape[1-1]+1).reshape(-1):
        iy = find(Prz(:,1) == np.logical_and(y00(i,1),Prz(:,2)) == y00(i,2))
        ix = find(yt == y00(i,4))
        map[ix,iy] = D_N(i)
    
    if np.all(np.all(map == - 2)):
        map = 100 * np.ones((len(yt),Prz.shape[1-1]))
    else:
        map[find[map == - 2]] = np.ones((map(find(map == - 2)).shape,map(find(map == - 2)).shape)) * np.amax(np.amax(map)) * 1.1
    
    subplot(2,1,2)
    plt.plot(bo(2,:),bo(1,:),'r',Prz(:,2),Prz(:,1),'or')
    vec = np.arange(1,Prz.shape[1-1]+3,3)
    if vec(len(vec)) != Prz.shape[1-1]:
        vec = np.array([vec,Prz.shape[1-1]])
    
    dip = 0.0004
    for i in vec.reshape(-1):
        text(Prz(i,2) + np.abs(dip * 0.5),Prz(i,1) + dip,num2str(i))
        dip = - dip
    
    plt.xlabel('z axis [m]')
    plt.ylabel('r axis [m]')
    plt.title('Initial points')
    map[1,1] = map(1,2) + 2 * eps
    subplot(2,1,1)
    pcolor(np.arange(1,Prz.shape[1-1]+1),yt * par(1) * 360,map)
    shading('flat')
    colormap(hot)
    matr = colormap
    iso = 1.1 * np.amax(np.amax(D_N))
    caxis(np.array([0,iso]))
    colorbar
    colormap(matr ** (1 / 3))
    plt.ylabel('Initial phase [deg]')
    plt.xlabel('Place referring to picture below')
    plt.title(np.array(['MultiPac 2.0          ',ti,'              ',date]))
    iy0 = []
    if trajec == 1:
        x,y = ginput
        x = x(len(x))
        iy0 = np.zeros((1,len(x)))
        for i in np.arange(1,len(x)+1).reshape(-1):
            tmp = np.abs(yt * par(1) * 360 - y(i))
            y[i] = np.amin(find(tmp == np.amin(tmp)))
            xd = np.array([np.amax(np.array([np.round(x(i) - 2),1])),np.amin(np.array([np.round(x(i) + 2),map.shape[2-1]]))])
            yd = np.array([np.amax(np.array([np.round(y(i) - 2),1])),np.amin(np.array([np.round(y(i) + 2),map.shape[1-1]]))])
            smap = map(np.arange(yd(1),yd(2)+1),np.arange(xd(1),xd(2)+1))
            mi = np.amin(np.amin(smap))
            J,I = find(smap == mi)
            indx = xd(1) + I(1) - 1
            indy = yd(1) + J(1) - 1
            tmp = find(D_N == mi)
            if map(indy,indx) == mi:
                iy0[i] = find(y00(:,1) == np.logical_and(Prz(indx,1),y00(:,2)) == np.logical_and(Prz(indx,2),y00(:,4)) == yt(indy))
    
    # -------------------------------------------------------------------------
    return iy0