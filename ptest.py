#!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python


import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys

import scipy.sparse.linalg as la_sparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=2)
np.set_printoptions(linewidth=300)


L = 3.0
W = 3.0

# grid
nx = 20
ny = 20
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

vsp = 3.0
hsp = 3.0


# Build a mesh using uniform x and y-spacing
# (uniformtity not required)
# Simple zoner function makes left half Silicon, right half Oxide
# Properties dictionary defines permittivities for the two materials
mesh = FVM.FVMMesh( np.linspace(0.0, L, nx), 
					np.linspace(0.0, W, ny),
					lambda x, y: 'Silicon' if x<0.5*L else 'Oxide', 
					{ 'Silicon' : { 'eps' : 11.9 },
					  'Oxide'   : { 'eps' : 3.9  } } ) 
					  
# Assign contacts (Dirichlet BCs) to left and right edges			
mesh.assignContact( lambda x, y: x<1e-8,   'Left' )
mesh.assignContact( lambda x, y: L-x<1e-8, 'Right')
					
print(mesh)

poisson = pois.FVMPoisson(mesh)

#print 'Lap mat:'
pmat =  poisson.Lap_mat
#print pmat.toarray()


# Now build rhs mat
bc_val_dict = { 'Left' : 0.5, 'Right' : 1.0 }
rhs = poisson.buildRHS(bc_val_dict)
#print 'RHS: ', rhs
print 'RHS shape: ', rhs.shape
print 'Poisson shape: ', pmat.shape


# Solve the resulting system
z = la_sparse.spsolve(pmat, rhs)
#print 'z:', z

(X, Y, Z) = poisson.pmesh.mapSolutionToGrid(z)


# Make some plots


fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z), np.amax(Z), 30)
cp = plt.contourf(X, Y, Z, levels=levels)
plt.contour(X, Y, Z, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	

'''
ax.add_patch(
    patches.Rectangle(
        (0.0, 0.0),
        L,
        W,
        fill=False,      # remove background
        linewidth=2
    )
)
'''

plt.show()	




