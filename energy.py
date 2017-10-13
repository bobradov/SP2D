#!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python



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
L = 5.0
W = 20.0

# Create Delta4 Valley
mstar = np.array([[0.19, 0.478, 0],[0.478, 0.315, 0],[0, 0, 0.315]])
mult  = 4
Delta4 = valley.valley('Delta4', L, W, mstar, mult)

# Create Delta2 Valley
mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
mult = 2
Delta2 = valley.valley('Delta2', L, W, mstar, mult)
# grid
nx = 50
ny = 50
fdmesh = fdmesh.FDMesh(L, W, nx, ny)


# Create Schroedinger-Poisson Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Delta4, Delta2], fdmesh)

# Do self-consistent calculation
# Compute self-consistent Fermi level

peakvals = np.linspace(0.025, 0.2, 10)

charge_array = np.zeros(10)
surf_array   = np.zeros(10)
fermi_array  = np.zeros(10)
valley_charge = np.zeros((2,10))
D4_gs         = np.zeros(10)
D2_gs         = np.zeros(10)
D4_fe         = np.zeros(10)
D2_fe         = np.zeros(10)

print 'Sweeping ...'

for index, cur_peak in enumerate(peakvals):
    # Create a CB profile
    cbfun = lambda x, y : x * (L-x) * cur_peak * 4.0/(L*L)
    curvature = 1e18 * cur_peak * 2.0 * 4.0/(L*L)
    charge_density = curvature * phys.eps0 * phys.eps_Si * 1e-6 / phys.q
    total_charge = charge_density * L * W * 1e-7 * 1e-7
    surface_charge = 0.5 * total_charge / (W * 1e-7)

    # Now solve system
    print 'Starting eigensolve ...'
    sp.solveWithCB(cbfun, n_eigen=40)


    # Ground states
    D4_gs[index] = sp.getSolution(0)['evals'][0]
    D2_gs[index] = sp.getSolution(1)['evals'][0]
    # First excited states
    D4_fe[index] = sp.getSolution(0)['evals'][1]
    D2_fe[index] = sp.getSolution(1)['evals'][1]




plt.figure()
plt.plot(peakvals, D4_gs, '-r', linewidth=3)
plt.plot(peakvals, D2_gs, '-b', linewidth=3)
plt.plot(peakvals, D4_fe, '--r', linewidth=3)
plt.plot(peakvals, D2_fe, '--b', linewidth=3)
plt.grid(True)
plt.xlabel('Band Bending [eV]')
plt.ylabel('Energy Levels')



plt.show()

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
