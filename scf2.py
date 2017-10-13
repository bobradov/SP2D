#!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import valley as valley
import SP2D as sp2d
import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys

import scipy.sparse.linalg as la_sparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib import colors, ticker, cm

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=300)


#+------------------------------------------------------------+
#|
#| Create Structure
#|
#+------------------------------------------------------------+

L = 5.0
W = 5.0
IL = 0.5

# grid
nx = 30
ny = 30
nh = 20
nv = 15
no = 5

#
# Finite-Difference Mesh: Used by Schroedinger
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.0
EOT = 0.85
vsp = 2.0
hsp = IL + (EOT-IL)*Eps_HiK/3.9

pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)
							
#pgaa_zoner = zoner.finFETZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)
							
#pgaa_zoner = zoner.NWZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)

#
# Finite-Volume Mesh: used by Poisson
mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			

# Assign the contact, call it 'Gate')			
mesh.assignContact( pgaa_zoner.contact_func, 'Gate' )

#+---------------------------------------------------------+
#|
#| Poisson-related Objects
#|
#+---------------------------------------------------------+

# Create Poisson Object
poisson = pois.FVMPoisson(mesh)

# Extract Laplace matrix; fixed dueing simulation
pmat =  poisson.Lap_mat

# Set Boundary Values on contacts
#bc_val_dict = { 'Left' : 0.6, 'Right' : 0.6 }
bc_val_dict = { 'Gate' : 0.6 }


#+---------------------------------------------------------+
#|
#| Schroedinger-related Objects
#|
#+---------------------------------------------------------+

n_eigen = 20  # Number of eigenvalues to use
Ef      = 0.0 # Assign Fermi level to 0.0

# Create Delta4 Valleys
# 'Left' Valleys
mstar = np.array([[0.315, 0.478, 0],[0.478, 0.19, 0],[0, 0, 0.315]])
mult  = 2
Delta4L = valley.valley('Delta4L', L, W, mstar, mult)

# 'Right' Valleys
mstar = np.array([[0.315, -0.478, 0],[-0.478, 0.19, 0],[0, 0, 0.315]])
mult  = 2
Delta4R = valley.valley('Delta4R', L, W, mstar, mult)

# Create Delta2 Valley
mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
mult = 2
Delta2 = valley.valley('Delta2', L, W, mstar, mult)

# Create Schroedinger-Poisson Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Delta4L, Delta4R, Delta2], fdmesh)

#+---------------------------------------------------------+
#|
#| Initial Solver: Poisson with initial guess for charge
#|                 Followed by QM-soltion with guessed potential
#|
#+---------------------------------------------------------+


# Initial guess for charge is a volume-inverted Gaussian
peak = 4e19
xc = L/2.0
yc = W/2.0
rhs_initial = poisson.buildRHS(
	bc_val_dict, 
	charge_func=lambda x, y : peak * np.exp(-(x-xc)*(x-xc)/(L*L) -(y-yc)*(y-yc)/(W*W))
)


# Solve the Poisson system
initial_potential     = la_sparse.spsolve(pmat, rhs_initial)
(X, Y, Z_initial_pot) = poisson.pmesh.mapSolutionToGrid(initial_potential)
prev_potential        = initial_potential           # Initial potential array

# Initial Schroedinger solution
poisson.pmesh.loadAccumulator(initial_potential)
cbfun = lambda x, y: -1.0*poisson.pmesh.nn_interp(x, y)
sp.solveWithCB(cbfun, n_eigen=n_eigen)

# Build QM charge density
wfn_dict          = sp.getSolution(0)
(Xfd,Yfd,Z_eigen) = wfn_dict['wfns']
ZCharge_init      = sp.computeChargeDensity(Ef)
initial_charge    = np.reshape(ZCharge_init, nx*ny) # Initial charge array

#--------------- Done with initial solution



#+------------------------------------------------------------+
#|
#| Self-consistency SP Loop
#|
#+------------------------------------------------------------+
prev_charge    = initial_charge
prev_potential = initial_potential 
max_iter = 200
tol      = 1e-4
alpha    = 0.1

for cur_iter in range(1, max_iter):
	# Prepare for Poisson solve
	# Load previous Schroediner solution into accumulator
	# Get estimated new potential using full previous charge
	# Then "damp" by taking only a fraction of
	# the potential as the update
	fdmesh.loadAccumulator(prev_charge)
	rhs = poisson.buildRHS(
		bc_val_dict, 
		charge_grid=fdmesh
	)
	pred_potential = la_sparse.spsolve(pmat, rhs)
	delta_pot =  pred_potential-prev_potential
	norm = np.sqrt((delta_pot*delta_pot).sum())
	if cur_iter % 5 == 0 or cur_iter == 1:
		print ' >>> Iter: ', cur_iter, ' Norm: ', norm
	
	# Now do Schroedinger solve with damped new potential
	new_potential  = alpha * pred_potential + (1.0-alpha)*prev_potential
	prev_potential = new_potential # Save potential for next iter
	
	poisson.pmesh.loadAccumulator(new_potential)
	cbfun = lambda x, y: -1.0*poisson.pmesh.nn_interp(x, y)
	sp.solveWithCB(cbfun, n_eigen=n_eigen)
	
	wfn_dict          = sp.getSolution(0)
	(Xfd,Yfd,Z_eigen) = wfn_dict['wfns']
	ZCharge           = sp.computeChargeDensity(Ef)
	new_charge        = np.reshape(ZCharge, nx*ny) # Initial charge array 
    
	prev_charge = new_charge
	
	if norm < tol:
		print 'Converged after ', cur_iter, ' iterations.'
		break


# Done with loop, plot results
fig1 = plt.figure(1)

# Potential sub-figure
(Xp, Yp, Z_final_pot) = poisson.pmesh.mapSolutionToGrid(new_potential)

ax = fig1.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z_final_pot), np.amax(Z_final_pot), 30)
cp = plt.contourf(Xp, Yp, Z_final_pot, levels=levels)
plt.contour(X, Y, Z_final_pot, levels=levels, linewidths=0.5, colors='k')
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

# QM Charge sub-figure
fig2 = plt.figure(2)
x2 = fig2.add_subplot(111, aspect='equal')
levels_c = np.linspace(np.amin(ZCharge), np.amax(ZCharge), 30)
#levels_c = np.linspace(5e18, 1e20, 30)
cp_QM = plt.contourf(Xfd, Yfd, ZCharge, levels=levels_c)
plt.colorbar(cp_QM)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')

# Get DOS for each valley
(energy, ret_val) = sp.computeDOS()
D4_DOS = ret_val[0]  # Symmetric D4 valleys, get only one, multiply by 2
D2_DOS = ret_val[2]  # Single D2 valley
plt.figure(3)
plt.plot(energy, 2*D4_DOS, '-r', linewidth=3)
plt.plot(energy, D2_DOS, '-b', linewidth=3)
plt.xlabel('Energy [eV]')
plt.ylabel('DOS')
plt.grid(True)

# Plot eigenfunctions
# Separate plot for each valley
for cur_valley in range(0,3):
	val_dict = sp.getSolution(cur_valley)
	(X,Y,Z) = val_dict['wfns']
	(rows, cols, slabs) = Z.shape
	nrows = 4
	ncols = 3
	div = (slabs-(slabs % ncols))/ncols + 1
	if div > nrows:
		div = nrows
	print 'div: ', div
	fig, axarr = plt.subplots(div, ncols)
	
	for row in range(0,div):
		for col in range(0,ncols):
			cur_slab = row * ncols + col
			if cur_slab == slabs:
			    break
			plotZ = Z[:,:,cur_slab] * Z[:,:,cur_slab]
			
			levels = np.linspace(np.amin(plotZ), np.amax(plotZ), 30)
			axarr[row, col].contourf(X, Y, plotZ, levels=levels)
			axarr[row, col].set_aspect('equal')

plt.show()


