# Function program plot_mesh(s)
# --------------------------------------------------------------------
# Plots the mesh.
#
# --------------------------------------------------------------------
# CALLS TO : error_message.m, clear_window.m
# xx/yy/00 : Seppo J‰rvemp‰‰
# 04/04/00 : Pasi Yl‰-Oijala - m-file
# --------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
    
def plot_mesh(s = None): 
    if len(varargin) < 1:
        s = 1
    
    cl = 0
    ok1 = exist('mesh.mat')
    ok2 = exist('fieldparam')
    if ok1 == 0:
        cl = error_message(np.array(['The mesh does not exists. Choose Mesh Generator ','in menu Run.']))
    else:
        if ok2 == 0:
            cl = error_message('Parameter file fieldparam does not exist.')
        else:
            scipy.io.loadmat('fieldparam')
            gtype = fieldparam(1)
            if gtype == 1:
                if s > 0:
                    cl = error_message('Plotting the mesh.')
                    #      error_message('                  ');
            else:
                if s > 0:
                    cl = error_message('Plotting the mesh. Blue area for streching.')
                    #      error_message('                  ');
            Window = plt.figure(3)
            clf
            set(Window,'name','MULTIPAC - Input Window II')
            # plots 2-d mesh in mesh.mat, which includes coord and etopol -arrays
            scipy.io.loadmat('mesh')
            xmin = np.amin(coord(:,1))
            xmax = np.amax(coord(:,1))
            ymin = np.amin(coord(:,2))
            ymax = np.amax(coord(:,2))
            plt.plot(np.array([xmin,xmax,xmax,xmin]),np.array([ymin,ymin,ymax,ymax]),'k.')
            hold('on')
            koko,pois = etopol.shape
            ala = 0.0
            X = np.zeros((4,koko))
            Y = np.zeros((4,koko))
            C = np.zeros((1,koko))
            for i in np.arange(1,koko+1).reshape(-1):
                n1 = etopol(i,1)
                x1 = coord(n1,1)
                y1 = coord(n1,2)
                n2 = etopol(i,2)
                x2 = coord(n2,1)
                y2 = coord(n2,2)
                n3 = etopol(i,3)
                x3 = coord(n3,1)
                y3 = coord(n3,2)
                X[:,i] = np.transpose(np.array([x1,x2,x3,x1]))
                Y[:,i] = np.transpose(np.array([y1,y2,y3,y1]))
                osa = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
                ala = ala + osa
                plt.plot(np.array([x1,x2,x3,x1]),np.array([y1,y2,y3,y1]),'k')
            I = find(alue == 0)
            fill(X(:,I),Y(:,I),'b')
            I = find(tyyppi == 1)
            fill(X(:,I),Y(:,I),'r')
            I = find(tyyppi == 2)
            fill(X(:,I),Y(:,I),'g')
            I = find(edges(:,3) > 0)
            for i in np.arange(1,len(I)+1).reshape(-1):
                i1 = edges(I(i),1)
                i2 = edges(I(i),2)
                plt.plot(np.array([coord(i1,1),coord(i2,1)]),np.array([coord(i1,2),coord(i2,2)]),'r')
            I = find(boundary(1,4) == 3)
            for i in np.arange(1,len(I)+1).reshape(-1):
                plt.plot(coord(I(i),1),coord(I(i),2),'b*')
            I = find(boundary(1,4) == 0)
            plt.plot(coord(I,1),coord(I,2),'w*')
            ala = ala / 2
            plt.title(np.array(['MultiPac 2.0                    Mesh                 ',date]))
            hold('off')
            plt.xlabel('z axis [m]')
            plt.ylabel('r axis [m]')
            plt.axis('image')
            colormap(jet)
    
    if cl == 1:
        clear_window
    
    # --------------------------------------------------------------------