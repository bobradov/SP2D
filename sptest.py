#! /opt/local/bin/python2.7



import valley as valley
import FDMesh as fdmesh
import Schroedinger as sch

import numpy  as np
import scipy.interpolate as interp
import scipy.integrate   as integral
import scipy.optimize    as optimize
import phys

import SP2D as sp2d

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

# Create Delta4 Valley
mstar = np.array([[0.19, 0.478, 0],[0.478, 0.315, 0],[0, 0, 0.315]])
mult  = 4
Delta4 = valley.valley(L, W, mstar, mult)

# Create Delta2 Valley
mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
mult = 2
Delta2 = valley.valley(L, W, mstar, mult)
# grid
nx = 50
ny = 50
fdmesh = fdmesh.FDMesh(L, W, nx, ny)


# Create Schroedinger-Poisson Solver
# Include all valleys
sp = sp2d.SP2D([Delta4, Delta2], fdmesh)


# Create a CB profile
peak = 0.100
cbfun = lambda x, y : x * (L-x) * peak * 4.0/(L*L)
curvature = 1e18 * peak * 2.0 * 4.0/(L*L)
charge_density = curvature * phys.eps0 * phys.eps_Si * 1e-6 / phys.q
print 'charge density: ', charge_density
total_charge = charge_density * L * W * 1e-7 * 1e-7
surface_charge = 0.5 * total_charge / (W * 1e-7)


# compute surface charge density
print 'surface charge dens: ', surface_charge


# Now solve system
sp.solveWithCB(cbfun, n_eigen=12)

res = sp.getSolution(valley_index=0)
evals = res['evals']
print 'evals: ', evals

# Compute DOS
(energy, DOS_array ) = sp.computeDOS(n_energy_points=300)

figDOS = plt.figure()
colors = ['-b','-r','-m']
for index, cur_valley in enumerate(DOS_array):
    plt.plot(energy, cur_valley,
                colors[index], linewidth=3)
plt.plot(energy, DOS_array[0], '-b', linewidth=2)
plt.xlabel('Energy [eV]')
plt.ylabel('DOS [1/eV cm]')
plt.grid(True)


# Charge density
fmax = 0.3
(fermi_vals, charge, valley_charge) = sp.computeChargePlot(fmax)




plt.figure()
plt.plot(fermi_vals, charge, '-k', linewidth=3)
(n_valleys, rows) = valley_charge.shape
print 'Valley shape: ', valley_charge.shape
colors = ['-b','-r','-m']
for cur_valley in range(0, n_valleys):
    plt.plot(fermi_vals, valley_charge[cur_valley,:],
                colors[cur_valley], linewidth=3)
plt.xlabel('Fermi Level [eV]')
plt.ylabel('Charge Density [1/cm^2]')
plt.grid(True)
plt.show()


# Do self-consistent calculation
# Compute self-consistent Fermi level
(charge, valley_charge) = sp.computeCharge( 0.3 )
print 'charge: ', charge, ' valley_charge: ', valley_charge
RHS = lambda x : sp.computeCharge( x )[0] - surface_charge
x0 = 0.15

xf= optimize.fsolve(RHS, x0, xtol=1e-3)
print 'Found Fermi level : ', xf


exit()

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
