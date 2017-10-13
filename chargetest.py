#!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python


import FDMesh as fdmesh
import PoissonMesh as pmesh
import Poisson     as pois
import numpy  as np
import phys

import scipy.sparse.linalg as la_sparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=1)
np.set_printoptions(linewidth=300)


L = 5.0
W = 6.0

# grid
nx = 50
ny = 50
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

vsp = 2.0
hsp = 1.0

p_mesh = pmesh.PoissonMesh(fdmesh, vsp, hsp, nh=15, nv=15,
							eps_semi=11.9, eps_ox=2.2)

# Find interface nodes
xi = p_mesh.getXInterfaces()
yi = p_mesh.getYInterfaces()
print 'xi: ', xi
print 'yi: ', yi

# Make a Poisson solver object
poisson = pois.Poisson(p_mesh)

#print 'Lap mat:'
pmat =  poisson.Lap_mat
#print pmat.toarray()

# Now build rhs mat
rhs = poisson.buildRHS()
#print 'RHS: ', rhs
print 'RHS shape: ', rhs.shape
print 'Poisson shape: ', pmat.shape

# Solve the resulting system
z = la_sparse.spsolve(pmat, rhs)
#print 'x:', z

# Make some plots
xvals  = p_mesh.getXMesh()
yvals = p_mesh.getYMesh()
X, Y = np.meshgrid(xvals, yvals)
Z = z.reshape(p_mesh.ny, p_mesh.nx)

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z), np.amax(Z), 30)
cp = plt.contourf(X, Y, Z, levels=levels)
plt.contour(X, Y, Z, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	

ax.add_patch(
    patches.Rectangle(
        (0.0, 0.0),
        L,
        W,
        fill=False,      # remove background
        linewidth=2
    )
)


plt.show()	




