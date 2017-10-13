#!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import ParabolicValley as pvalley
import SP2D as sp2d
import SCF  as SCF
import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys
from   plot_fun import set_plt_defaults
import SuperLattice as super_l


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

# Fin dimensions
L = 5.0
W_target = 45.0
IL = 0.5


# Superlattice
thick     = 0.5
dy_target = 0.25*thick
period    = 5.0*thick
n_periods = np.floor(W_target/period)
W         = n_periods * period
offset    = 0
energy    = 0.0
sl = super_l.SuperLattice(period, thick, energy, offset)




print 'Height=', W, ' period=', period

#dy_target = 5.0/30

ny = int(W/dy_target)
print 'ny: ', ny
#exit()

# grid
nx = 30
#ny = 50
nh = 20
nv = 15
no = 5

#
# Finite-Difference Mesh: Used by Schroedinger
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.2
EOT = 0.85
#vsp = 2.4923
vsp = 2.5
hsp = IL + (EOT-IL)*Eps_HiK/3.9
print 'hsp: ', hsp

'''
pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)
'''
							
#pgaa_zoner = zoner.MLFETZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)
							
pgaa_zoner = zoner.finFETZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)
							
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
bc_val_dict = { 'Gate' : 0.7 }


#+---------------------------------------------------------+
#|
#| Schroedinger-related Objects
#|
#+---------------------------------------------------------+


# Create Delta4 Valleys
# 'Left' Valleys
mstar = np.array([[0.07, 0, 0],[0, 0.07, 0],[0, 0, 0.07]])
mult  = 1
Gamma = pvalley.ParabolicValley('Gamma', L, W, mstar, mult)



# Create Schroedinger Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Gamma], fdmesh)

#+---------------------------------------------------------+
#|
#| SCF Solver: Self-Consistent Solution of QM and Poisson
#|
#+---------------------------------------------------------+

scf = SCF.SCF( sp, poisson, sl )
Ef = 0.0
(iterations, converged) = scf.solve(
								bc_val_dict, 
								max_iter = 200, 
								tol      = 2e-3, 
								alpha    = 0.1,
								n_eigen  = 20,
								Ef       = Ef 
							)

print 'Converged=', converged, ' after ', iterations, ' iterations.'


#set_plt_defaults()
plt.rcParams['figure.facecolor'] = 'w'
plt.rc('xtick', labelsize=18) 
plt.rc('ytick', labelsize=18) 
font = {'family' : 'Bitstream Vera Sans',
	    'weight' : 'bold',
	    'size'   : 16}

plt.rc('font', **font)
plt.rcParams['toolbar'] = 'None'
#plt.tick_params(which='major', length=15, width=2)
#plt.tick_params(which='minor', length=5, width=2)

plt.rcParams['axes.linewidth'] = 2 

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

for cur_rect in pgaa_zoner.plotting_rectangle_tuple():
	start_point = cur_rect[0]
	length      = cur_rect[1]
	height      = cur_rect[2]
	ax.add_patch(
		patches.Rectangle(
		    start_point,
		    length,
		    height,
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
(energy, ret_val) = sp.computeDOS(e_min=-0.1, e_max=0.2)
Gamma_DOS = ret_val[0]  # Symmetric D4 valleys, get only one, multiply by 2
plt.figure(3)
#plt.plot(energy, 2*D4_DOS, '-r', linewidth=3)
plt.plot(energy, Gamma_DOS, '-b', linewidth=3)
plt.xlabel('Energy [eV]')
plt.ylabel('DOS')
plt.grid(True)
#plt.xlim(xmin=-0.1, xmax=0.2)
plt.ylim(ymin=0.0, ymax=100)

# Plot eigenfunctions
# Separate plot for each valley


for cur_valley in range(0,1):
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
			axarr[row, col].xaxis.set_visible(False)
			axarr[row, col].yaxis.set_visible(False)
			
			
# Compute charge and velocity

(vel, tot_charge, valley_charge) = sp.computeChargeAndVelocity(Ef)
sidewall_charge = 1e14 * tot_charge / (2.0*(W+vsp))
nw_charge       = 1e14 * tot_charge / (2.0*(W+vsp+L+hsp))
#d2_occupancy = valley_charge[2]/(tot_charge)

print 'vel=', vel*1e-7, ' tot_charge=', tot_charge
#print 'valley charge: ', valley_charge



#sw_total_charge = 1e-14*sw_charge[0]*(2.0*W)
#print 'Standard SW charge: ', sw_charge, ' SW total charge: ', sw_total_charge
plt.show()


