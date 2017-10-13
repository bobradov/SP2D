#!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import valley as valley
import SP2D as sp2d
import SCF  as SCF
import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys


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
W = 9.0
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
Eps_HiK  = 22.2
EOT = 0.85
vsp = 2.4923
hsp = IL + (EOT-IL)*Eps_HiK/3.9
print 'hsp: ', hsp

#pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)
							
#pgaa_zoner = zoner.finFETZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)
							
pgaa_zoner = zoner.NWZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)

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
bc_val_dict = { 'Gate' : 0.7 }


#+---------------------------------------------------------+
#|
#| Schroedinger-related Objects
#|
#+---------------------------------------------------------+


# Create Delta4 Valleys
# 'Left' Valleys
mstar = np.array([[0.315, 0.478, 0],[0.478, 0.19, 0],[0, 0, 0.55]])
mult  = 2
Delta4L = valley.valley('Delta4L', L, W, mstar, mult)

# 'Right' Valleys
mstar = np.array([[0.315, -0.478, 0],[-0.478, 0.19, 0],[0, 0, 0.55]])
mult  = 2
Delta4R = valley.valley('Delta4R', L, W, mstar, mult)

# Create Delta2 Valley
mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
mult = 2
Delta2 = valley.valley('Delta2', L, W, mstar, mult)

# Create Schroedinger Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Delta4R, Delta4L, Delta2], fdmesh)

#+---------------------------------------------------------+
#|
#| SCF Solver: Self-Consistent Solution of QM and Poisson
#|
#+---------------------------------------------------------+

scf = SCF.SCF( sp, poisson )
Ef = 0.0
(iterations, converged) = scf.solve(
								bc_val_dict, 
								max_iter = 200, 
								tol      = 1e-4, 
								alpha    = 0.1,
								n_eigen  = 20,
								Ef       = Ef 
							)

print 'Converged=', converged, ' after ', iterations, ' iterations.'




# Done with loop, plot results
fig1 = plt.figure(1)

# Potential sub-figure
(Xp, Yp, Z_final_pot) = poisson.pmesh.mapSolutionToGrid(scf.final_potential)

ax = fig1.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z_final_pot), np.amax(Z_final_pot), 30)
cp = plt.contourf(Xp, Yp, Z_final_pot, levels=levels)
plt.contour(Xp, Yp, Z_final_pot, levels=levels, linewidths=0.5, colors='k')
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
wfn_dict          = sp.getSolution(0)
(Xfd,Yfd,Z_eigen) = wfn_dict['wfns']
ZCharge           = sp.computeChargeDensity(Ef)
levels_c = np.linspace(np.amin(ZCharge), np.amax(ZCharge), 30)
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
			
			
# Compute charge and velocity

(vel, tot_charge, valley_charge) = sp.computeChargeAndVelocity(Ef)
sidewall_charge = 1e14 * tot_charge / (2.0*(W+vsp))
nw_charge       = 1e14 * tot_charge / (2.0*(W+vsp+L+hsp))
d2_occupancy = valley_charge[2]/(tot_charge)

print 'vel=', vel*1e-7, ' tot_charge=', tot_charge, \
      ' sidewall=', sidewall_charge, ' d2: ', d2_occupancy, \
      ' nw charge: ', nw_charge
#print 'valley charge: ', valley_charge



#sw_total_charge = 1e-14*sw_charge[0]*(2.0*W)
#print 'Standard SW charge: ', sw_charge, ' SW total charge: ', sw_total_charge
plt.show()


