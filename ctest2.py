#!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

# Solver Poisson
# Charge is assigned to the semiconductor mesh directly
from scipy import interpolate

import Zoner   as zoner
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


L = 5.0
W = 5.0
IL = 0.5

# grid
nx = 30
ny = 30
nh = 20
nv = 15
no = 5
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.0
EOT = 0.85
vsp = 2.0
hsp = IL + (EOT-IL)*Eps_HiK/3.9

pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)


mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			
mesh.assignContact( lambda x, y: x < -hsp + 1e-8,   'Left' )
mesh.assignContact( lambda x, y: x >  L+hsp-1e-8, 'Right')
					
#print(mesh)

#exit()

poisson = pois.FVMPoisson(mesh)

#print 'Lap mat:'
pmat =  poisson.Lap_mat
#print pmat.toarray()


# Now build rhs mat
bc_val_dict = { 'Left' : 0.7, 'Right' : 0.7 }
peak = 8e19
xc = L/2.0
yc = W/2.0
charge_func=lambda x, y : peak * np.exp(-(x-xc)*(x-xc)/(L*L) -(y-yc)*(y-yc)/(W*W))

# Mesh for the semiconductor
xsemi = np.linspace(0.0, L, nx)
ysemi = np.linspace(0.0, W, ny)
XSemi, YSemi = np.meshgrid(xsemi, ysemi)
csemi = charge_func(XSemi, YSemi)
print csemi.shape
csemi_array = np.reshape(csemi, nx*ny)

fdmesh.loadAccumulator(csemi_array)



rhs = poisson.buildRHS(
	bc_val_dict, 
	charge_grid=fdmesh
)
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


ax.add_patch(
    patches.Rectangle(
        (0.0, 0.0),
        L,
        W,
        fill=False,      # remove background
        linewidth=2
    )
)

ax.add_patch(
    patches.Rectangle(
        (-IL, -IL),
        L+2*IL,
        W+2*IL,
        fill=False,      # remove background
        linewidth=2
    )
)


plt.show()	




