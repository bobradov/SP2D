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
W = 5.0
IL = 0.5

# grid
nx = 50
ny = 50
nh = 20
nv = 15
no = 5

#
# Finite-Difference Mesh: Used by Schroedinger
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.2
EOT = 0.85
vsp = 2.0
hsp = IL + (EOT-IL)*Eps_HiK/3.9

dev_type='pgaa'

if dev_type=='pgaa':
	pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
								vsp, hsp, nh, nv, no, 
								IL=IL, Eps_HiK=Eps_HiK)
elif dev_type=='fin':
	pgaa_zoner = zoner.finFETZoner(L, W, fdmesh, 
								vsp, hsp, nh, nv, no, 
								IL=IL, Eps_HiK=Eps_HiK)
elif dev_type=='nw':							
	pgaa_zoner = zoner.NWZoner(L, W, fdmesh, 
								vsp, hsp, nh, nv, no, 
								IL=IL, Eps_HiK=Eps_HiK)
else:
	print 'Unknown device type: ', dev_type
	exit()

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
Vg = 0.6
bc_val_dict = { 'Gate' : Vg }


#+---------------------------------------------------------+
#|
#| Schroedinger-related Objects
#|
#+---------------------------------------------------------+


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

# Create Schroedinger Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Delta4L, Delta4R, Delta2], fdmesh)

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

			
# Compute charge and velocity

(vel, tot_charge, valley_charge) = sp.computeChargeAndVelocity(Ef)
sidewall_charge = 1e14 * tot_charge / (2.0*(W+vsp))
nw_charge       = 1e14 * tot_charge / (2.0*(W+vsp+L+hsp))
d2_occupancy = valley_charge[2]/(tot_charge)

# Save everything to 'FOM.dat'
fom = open('FOM.dat', 'w')
fom.write('L ' + str(L) + '\n')
fom.write('W ' + str(W) + '\n')
fom.write('vsp ' + str(vsp) + '\n')
fom.write('Eps ' + str(Eps_HiK) + '\n')
fom.write('Vg  ' + str(Vg)  + '\n')
fom.write('sw_charge ' + str(sidewall_charge) + '\n')
fom.write('nw_charge ' + str(nw_charge) + '\n')
fom.write('tot_charge ' + str(tot_charge) + '\n')
fom.write('d2 ' + str(d2_occupancy) + '\n')
fom.write('vel ' + str(vel) + '\n')
fom.write('converged ' + str(converged) + '\n')

fom.close()






