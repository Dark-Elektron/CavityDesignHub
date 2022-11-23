# FUNCTION program X = POINTS3(x,y,z)
# -----------------------------------------------------------------------------
# Program to generate a field-point-matrix X in R^3 for given coordinate-
# vectors x,y and z. Each column of X is a point in R^3 as follows:
#              X(:,h) = [x(i),y(j),z(k)]', h=1,...,nx*ny*nz.
# The indexing system is as follows:
#              leading index k runs from 1 to nz,
#              for every k, the secondary index j runs from 1 to ny and
#              for every j, the index i runs from 1 to nx.
# The vectors x, y and z can be given in rectangular, spherical (x=r,y=theta,
# z=phi) or cylindrical (x=r,y=phi,z=z) coordinates.
# INPUT   x,y,z : (1,nx), (1,ny) and (1,nz)-vectors specifying the x-,y- and
#                  z-coordinates, [m] or [rad]
# OUTPUT    X   : (3,nx*ny*nz)-matrix giving the field points, [m] or [rad]
# ---------------------------------------------------------------------------
# CALLS TO : None
# 08/10/91 : Pasi Yla-Oijala: Rolf Nevanlinna Institute
# ---------------------------------------------------------------------------

import numpy as np
    
def points3(x = None,y = None,z = None): 
    nx = len(x)
    ny = len(y)
    nz = len(z)
    X = np.zeros((3,nx * ny * nz))
    if nz == 1:
        zx = np.ones((1,nx)) * z
        for j in np.arange(1,ny+1).reshape(-1):
            yx = np.ones((1,nx)) * y(j)
            X[:,np.arange[[j - 1] * nx + 1,j * nx+1]] = np.array([[x],[yx],[zx]])
    else:
        if ny == 1:
            yx = np.ones((1,nx)) * y
            for k in np.arange(1,nz+1).reshape(-1):
                zx = np.ones((1,nx)) * z(k)
                X[:,np.arange[[k - 1] * nx + 1,k * nx+1]] = np.array([[x],[yx],[zx]])
        else:
            if nx == 1:
                xx = np.ones((1,ny)) * x
                for k in np.arange(1,nz+1).reshape(-1):
                    zy = np.ones((1,ny)) * z(k)
                    X[:,np.arange[[k - 1] * ny + 1,k * ny+1]] = np.array([[xx],[y],[zy]])
            else:
                ind = 1
                for k in np.arange(1,nz+1).reshape(-1):
                    zx = np.ones((1,nx)) * z(k)
                    for j in np.arange(1,ny+1).reshape(-1):
                        yx = np.ones((1,nx)) * y(j)
                        X[:,np.arange[nx * [ind - 1] + 1,nx * ind+1]] = np.array([[x],[yx],[zx]])
                        ind = ind + 1
    
    # -----------------------------------------------------------------------
    return X