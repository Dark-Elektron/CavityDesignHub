# Author: Prabhu Ramachandran <prabhu@aero.iitb.ac.in>
# Copyright (c) 2007 Prabhu Ramachandran.
# License: BSD Style.
import random

import mayavi
import numpy as np
from tvtk.api import tvtk
from mayavi.scripts import mayavi2


# @mayavi2.standalone
# def main():
#     # Create some random points to view.
#     pd = tvtk.PolyData()
#     pd.points = np.random.random((1000, 3))
#     verts = np.arange(0, 1000, 1)
#     verts.shape = (1000, 1)
#     pd.verts = verts
#     pd.point_data.scalars = np.random.random(1000)
#     pd.point_data.scalars.name = 'scalars'
#
#     # Now visualize it using mayavi2.
#     from mayavi.sources.vtk_data_source import VTKDataSource
#     from mayavi.modules.outline import Outline
#     from mayavi.modules.surface import Surface
#
#     mayavi.new_scene()
#     d = VTKDataSource()
#     d.data = pd
#     mayavi.add_source(d)
#     mayavi.add_module(Outline())
#     s = Surface()
#     mayavi.add_module(s)
#     s.actor.property.trait_set(representation='p', point_size=2)
#     # You could also use glyphs to render the points via the Glyph module.
#
# if __name__ == '__main__':
#     main()

# x = [random.randint(0,10)/10.0 for i in range(100)]
# y = [random.randint(0,10)/10.0 for i in range(100)]
# # x, y = np.meshgrid(x,y)
# # x = [xi for xj in x for xi in xj]
# # y = [yi for yj in y for yi in yj]
# z = [random.randint(0,10)/10.0 for i in range(100)]
# from mayavi import mlab
# s = mlab.points3d(x, y, z, mode = 'sphere', extent=[0,1,0,1,0,1])
# mlab.axes(s, ranges = [min(x), max(x), min(y), max(y), min(z), max(z)])
# mlab.show()
import numpy as np
from scipy.spatial import Delaunay
from mayavi import mlab

X2 = np.array([0, 0, 1, 1])
Y2 = np.array([0.5, 0.45, 1, 0.5])
Z2 = np.array([0, 1, 0.5,0])

# use scipy for delaunay:
p2d = np.vstack([X2,Y2]).T
d2d = Delaunay(p2d)

fig = mlab.figure(1, bgcolor=(1, 0.7, 1), fgcolor=(0.5, 0.5, 0.5))

# Generate triangular Mesh:
tmesh = mlab.triangular_mesh(X2, Y2, Z2, d2d.vertices,
                             scalars=Y2, colormap='jet')

# Simple plot.
mlab.outline(extent=(0,1,0,1,0,1))
mlab.axes(extent=(0,1,0,1,0,1))
mlab.show()