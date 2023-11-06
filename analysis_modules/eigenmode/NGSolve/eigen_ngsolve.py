from ngsolve import *
from ngsolve.webgui import Draw
from netgen.occ import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy

mu0 = 4 * pi * 1e-7
eps0 = 8.85418782e-12
c0 = 299792458


def eigen2d(geometry_dir):
    # define geometry
    cav_geom = pd.read_csv('D:\Dropbox\multipacting\MPGUI21\geodata.n',
                           header=None, skiprows=3, skipfooter=1, sep='\s+', engine='python')[[1, 0]]

    pnts = list(cav_geom.itertuples(index=False, name=None))
    wp = WorkPlane()
    wp.MoveTo(*pnts[0])
    for p in pnts[1:]:
        wp.LineTo(*p)
    wp.Close().Reverse()
    face = wp.Face()

    # name the boundaries
    face.edges.Max(X).name = "r"
    face.edges.Max(X).col = (1, 0, 0)
    face.edges.Min(X).name = "l"
    face.edges.Min(X).col = (1, 0, 0)
    face.edges.Min(Y).name = "b"
    face.edges.Min(Y).col = (1, 0, 0)

    geo = OCCGeometry(face, dim=2)

    # mesh
    ngmesh = geo.GenerateMesh(maxh=0.01)
    mesh = Mesh(ngmesh)

    # define finite element space
    fes = HCurl(mesh, order=1, dirichlet='default')

    u, v = fes.TnT()

    a = BilinearForm(y * curl(u) * curl(v) * dx).Assemble()
    m = BilinearForm(y * u * v * dx).Assemble()

    apre = BilinearForm(y * curl(u) * curl(v) * dx + y * u * v * dx)
    pre = Preconditioner(apre, "direct", inverse="sparsecholesky")

    with TaskManager():
        a.Assemble()
        m.Assemble()
        apre.Assemble()
        freedof_matrix = a.mat.CreateSmoother(fes.FreeDofs())

        # build gradient matrix as sparse matrix (and corresponding scalar FESpace)
        gradmat, fesh1 = fes.CreateGradient()

        gradmattrans = gradmat.CreateTranspose()  # transpose sparse matrix
        math1 = gradmattrans @ m.mat @ gradmat  # multiply matrices
        math1[0, 0] += 1  # fix the 1-dim kernel
        invh1 = math1.Inverse(inverse="sparsecholesky", freedofs=fesh1.FreeDofs())

        # build the Poisson projector with operator Algebra:
        proj = IdentityMatrix() - gradmat @ invh1 @ gradmattrans @ m.mat

        projpre = proj @ pre.mat
        evals, evecs = solvers.PINVIT(a.mat, m.mat, pre=projpre, num=30, maxit=20, printrates=False);

    # print out eigenvalues
    print("Eigenvalues")
    freq_fes = []
    for i, lam in enumerate(evals):
        freq_fes.append(c0 * np.sqrt(lam) / (2 * np.pi) * 1e-6)
        print(i, lam, 'freq: ', c0 * np.sqrt(lam) / (2 * np.pi) * 1e-6, "MHz")

    # plot results
    gfu = GridFunction(fes, multidim=len(evecs))
    for i in range(len(evecs)):
        gfu.vecs[i].data = evecs[i]

    # # alternative eigenvalue solver
    # u = GridFunction(fes, multidim=30, name='resonances')
    # lamarnoldi = ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(),
    #                         list(u.vecs), shift=300)
    # print(np.sort(c0*np.sqrt(lamarnoldi)/(2*np.pi) * 1e-6))

    evaluate_qois(cav_geom, gfu, mesh, freq_fes)


def eigen3d(geometry_dir):
    # define geometry
    cav_geom = pd.read_csv('D:\Dropbox\multipacting\MPGUI21\geodata.n',
                           header=None, skiprows=3, skipfooter=1, sep='\s+', engine='python')[[1, 0]]

    pnts = list(cav_geom.itertuples(index=False, name=None))
    wp = WorkPlane()
    wp.MoveTo(*pnts[0])
    for p in pnts[1:]:
        wp.LineTo(*p)
    wp.Close().Reverse()
    face = wp.Face()

    # name the boundaries
    face.edges.Max(X).name = "r"
    face.edges.Max(X).col = (1, 0, 0)
    face.edges.Min(X).name = "l"
    face.edges.Min(X).col = (1, 0, 0)
    face.edges.Min(Y).name = "b"
    face.edges.Min(Y).col = (1, 0, 0)

    vol = face.Revolve(Axis((0,0,0), X), 360)
    Draw(vol)
    geo = OCCGeometry(vol, dim=3)

    # mesh
    ngmesh = geo.GenerateMesh(maxh=0.01)
    mesh = Mesh(ngmesh)

    # define finite element space
    fes = HCurl(mesh, order=1, dirichlet='default')

    u, v = fes.TnT()

    a = BilinearForm(y * curl(u) * curl(v) * dx).Assemble()
    m = BilinearForm(y * u * v * dx).Assemble()

    apre = BilinearForm(y * curl(u) * curl(v) * dx + y * u * v * dx)
    pre = Preconditioner(apre, "direct", inverse="sparsecholesky")

    with TaskManager():
        a.Assemble()
        m.Assemble()
        apre.Assemble()
        freedof_matrix = a.mat.CreateSmoother(fes.FreeDofs())

        # build gradient matrix as sparse matrix (and corresponding scalar FESpace)
        gradmat, fesh1 = fes.CreateGradient()

        gradmattrans = gradmat.CreateTranspose()  # transpose sparse matrix
        math1 = gradmattrans @ m.mat @ gradmat  # multiply matrices
        math1[0, 0] += 1  # fix the 1-dim kernel
        invh1 = math1.Inverse(inverse="sparsecholesky", freedofs=fesh1.FreeDofs())

        # build the Poisson projector with operator Algebra:
        proj = IdentityMatrix() - gradmat @ invh1 @ gradmattrans @ m.mat

        projpre = proj @ pre.mat
        evals, evecs = solvers.PINVIT(a.mat, m.mat, pre=projpre, num=30, maxit=20, printrates=False);

    # print out eigenvalues
    print("Eigenvalues")
    freq_fes = []
    for i, lam in enumerate(evals):
        freq_fes.append(c0 * np.sqrt(lam) / (2 * np.pi) * 1e-6)
        print(i, lam, 'freq: ', c0 * np.sqrt(lam) / (2 * np.pi) * 1e-6, "MHz")

    # plot results
    gfu = GridFunction(fes, multidim=len(evecs))
    for i in range(len(evecs)):
        gfu.vecs[i].data = evecs[i]

    # # alternative eigenvalue solver
    # u = GridFunction(fes, multidim=30, name='resonances')
    # lamarnoldi = ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(),
    #                         list(u.vecs), shift=300)
    # print(np.sort(c0*np.sqrt(lamarnoldi)/(2*np.pi) * 1e-6))


def evaluate_qois(cav_geom, gfu, mesh, freq_fes):
    n = 2
    w = 2 * pi * freq_fes[n] * 1e6
    beta = 1

    # calculate Vacc and Eacc
    Vcav = abs(Integrate(gfu.MDComponent(n)[0] * exp(1j * w / (beta * c0) * x), mesh, definedon=mesh.Boundaries('b')))
    Eacc = Vcav / (0.187 * 4)

    # calculate U and R/Q
    U = 2 * pi * 0.5 * eps0 * Integrate(y * InnerProduct(gfu.MDComponent(n), gfu.MDComponent(n)), mesh)
    RoQ = Vcav ** 2 / (w * U)
    print("R/Q: ", RoQ)

    # calculate peak surface fields
    xpnts = cav_geom[(cav_geom[0] > 0) & (cav_geom[1] > min(cav_geom[1])) & (cav_geom[1] < max(cav_geom[1]))]
    vals = [Norm(gfu.MDComponent(n))(mesh(xi, yi)) for indx, (xi, yi) in xpnts.iterrows()]
    Epk = (max(vals))
    print("Epk/Eacc", Epk / Eacc)


eigen2d('d')
