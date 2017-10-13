#!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import valley as valley
import FDMesh as fdmesh
import Schroedinger as sch

import numpy  as np
import scipy.interpolate as interp
import scipy.integrate   as integral
import scipy.optimize    as optimize
import phys

import scipy.sparse.linalg as la_sparse
import scipy.linalg        as la

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt


# Generate DOS for a single valley
# This is an example of a Delta4 valley in (110) Si
#n_points = 400
L = 6.0
W = 6.0
#mstar = np.array([[0.19, -0.478, 0],[-0.478, 0.315, 0],[0, 0, 0.315]])
mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
mult  = 1

# Create Delta2 Valley
Delta4 = valley.valley('Delta4', L, W, mstar, mult)
# grid
nx = 50
ny = 50
fdmesh = fdmesh.FDMesh(L, W, nx, ny)
SchrD2  = sch.Schroedinger(Delta4, fdmesh)

# Create a CB profile
peak = 0.250
cbfun = lambda x, y : x * (L-x) * peak * 4.0/(L*L)
curvature = 1e18 * peak * 2.0 * 4.0/(L*L)
charge_density = curvature * phys.eps0 * phys.eps_Si * 1e-6 / phys.q
print 'charge density: ', charge_density
total_charge = charge_density * L * W * 1e-7 * 1e-7
surface_charge = 0.5 * total_charge / (W * 1e-7)


# compute surface charge density
print 'surface charge dens: ', surface_charge


(evals, (X,Y,Z)) = SchrD2.solve(n_eigen=12, CB=cbfun)
print 'evals: ', evals

# Make some plots

# DOS
n_energy_points = 300
energy = np.linspace(evals[0]-0.05, evals[0]+0.1, n_energy_points)
DOS = Delta4.computeDOS(energy, evals)

'''
figDOS = plt.figure()
plt.plot(energy, DOS, '-b', linewidth=2)
plt.xlabel('Energy [eV]')
plt.ylabel('DOS [1/eV cm]')
plt.grid(True)
'''

# Charge density
Ef = 0.3
occDOS = DOS * Delta4.occupancy( energy, Ef )
DOS_interp = interp.interp1d( energy, DOS )


fermi_points = 50
fermi_vals  = np.linspace(energy[0]-3.0*phys.kT, energy[-1], fermi_points)
charge_vals = np.zeros(fermi_points)
for i in range(0, fermi_points):
    if i % 5 == 0:
        print 'Working on Fermi level:', fermi_vals[i]
    charge_vals[i] = Delta4.computeSidewallCharge(fermi_vals[i], evals)

'''
# Compute self-consistent Fermi level
RHS = lambda x : Delta4.computeQMCharge( x, evals ) - surface_charge
x0 = 0.15
xf= optimize.fsolve(RHS, x0, xtol=1e-3)
print 'Found Fermi level : ', xf
'''

# Now make some evil plots
plt.figure()
plt.plot(energy, DOS_interp(energy), '-b', linewidth=3)
plt.plot(energy, occDOS, '-r', linewidth=3)
plt.xlabel('Energy [eV]')
plt.ylabel('DOS')
plt.grid(True)


plt.figure()
plt.plot(fermi_vals, charge_vals, '-k', linewidth=3)
plt.xlabel('Fermi Level [eV]')
plt.ylabel('Charge Density [1/cm^2]')
plt.grid(True)

(rows, cols, slabs) = Z.shape
ncols = 3
div = (slabs-(slabs % ncols))/ncols + 1
print 'div: ', div
fig, axarr = plt.subplots(div,ncols)
print axarr.shape
for row in range(0,div):
    for col in range(0,ncols):
        cur_slab = row * ncols + col
        if cur_slab == slabs:
            break
        plotZ = Z[:,:,cur_slab] * Z[:,:,cur_slab]
        #plotZ = Z[:,:,cur_slab]

        levels = np.linspace(np.amin(plotZ), np.amax(plotZ), 30)
        axarr[row, col].contourf(X, Y, plotZ, levels=levels)




plt.show()
