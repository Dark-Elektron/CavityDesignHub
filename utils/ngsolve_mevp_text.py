#!/usr/bin/env python
# coding: utf-8

# In[17]:


from ngsolve import *
from ngsolve.webgui import Draw
from netgen.occ import *
from netgen.geom2d import SplineGeometry
import pandas as pd
import matplotlib.pyplot as plt

# In[10]:


# define geometry

geo = SplineGeometry()
# get cavity geometry
cav_geom = pd.read_csv('D:\Dropbox\multipacting\MPGUI21\mid_TESLA\geodata.n',
                       header=None, skiprows=3, skipfooter=1, sep='\s+', engine='python')[[1, 0]]

# points have to be defined anticlockwise
cav_geom = cav_geom.iloc[::-1]
print(cav_geom)

pnts = list(cav_geom.itertuples(index=False, name=None))
pts = [geo.AppendPoint(*pnt) for pnt in pnts]

curves = []
for count in range(1, len(pts)):
    curves.append([["line", pts[count - 1], pts[count]], f"l{count}"])

    # add last point to first
    if count + 1 == len(pts):
        curves.append([["line", pts[count], pts[0]], f"l{0}"])

[geo.Append(c) for c, bc in curves]

# In[11]:


ngmesh = geo.GenerateMesh(maxh=0.1)

# In[12]:


mesh = Mesh(ngmesh)
Draw(mesh)

# In[13]:
# SetHeapSize(100*1000*1000)

fes = HCurl(mesh, order=3)
print("ndof =", fes.ndof)
u, v = fes.TnT()

a = BilinearForm(curl(u) * curl(v) * dx)
m = BilinearForm(u * v * dx)

apre = BilinearForm(curl(u) * curl(v) * dx + u * v * dx)
pre = Preconditioner(apre, "direct", inverse="sparsecholesky")

# In[14]:


with TaskManager():
    a.Assemble()
    m.Assemble()
    apre.Assemble()

    # build gradient matrix as sparse matrix (and corresponding scalar FESpace)
    gradmat, fesh1 = fes.CreateGradient()

    gradmattrans = gradmat.CreateTranspose()  # transpose sparse matrix
    math1 = gradmattrans @ m.mat @ gradmat  # multiply matrices
    math1[0, 0] += 1  # fix the 1-dim kernel
    invh1 = math1.Inverse(inverse="sparsecholesky")

    # build the Poisson projector with operator Algebra:
    proj = IdentityMatrix() - gradmat @ invh1 @ gradmattrans @ m.mat

    projpre = proj @ pre.mat

    evals, evecs = solvers.PINVIT(a.mat, m.mat, pre=projpre, num=12, maxit=20)

# In[15]:


print("Eigenvalues")
for lam in evals:
    print(lam)

# In[16]:


gfu = GridFunction(fes, multidim=len(evecs))
for i in range(len(evecs)):
    gfu.vecs[i].data = evecs[i]

Draw(Norm(gfu), mesh, order=2)
