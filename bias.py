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

# Prepare for some evil plots
def set_plt_defaults():
	plt.rcParams['figure.facecolor'] = 'w'
	plt.rc('xtick', labelsize=18) 
	plt.rc('ytick', labelsize=18) 
	font = {'family' : 'Bitstream Vera Sans',
		    'weight' : 'bold',
		    'size'   : 16}

	plt.rc('font', **font)
	plt.rcParams['toolbar'] = 'None'

	plt.tick_params(which='major', length=15, width=2)
	plt.tick_params(which='minor', length=5, width=2)

	plt.rcParams['axes.linewidth'] = 2 
	

sidewall='110'

sizes = [5.0]
for plot_ind, curW in enumerate(sizes):
	# Generate DOS for a single valley
	# This is an example of a Delta4 valley in (110) Si
	#n_points = 400
	L = 5.0
	W = curW
	
	if sidewall=='110':
		# Create Delta4 Valleys
		mstar = np.array([[0.315, 0.478, 0],[0.478, 0.19, 0],[0, 0, 0.315]])
		mult  = 2
		Delta4L = valley.valley('Delta4L', L, W, mstar, mult)

		mstar = np.array([[0.315, -0.478, 0],[-0.478, 0.19, 0],[0, 0, 0.315]])
		mult  = 2
		Delta4R = valley.valley('Delta4R', L, W, mstar, mult)

		# Create Delta2 Valley
		mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
		mult = 2
		Delta2 = valley.valley('Delta2', L, W, mstar, mult)

	elif sidewall=='100':
		# Create Delta4 Valleys
		mstar = np.array([[0.19, 0, 0],[0, 0.19, 0],[0, 0, 0.92]])
		mult  = 2
		Delta4L = valley.valley('Delta4L', L, W, mstar, mult)

		mstar = np.array([[0.92, 0, 0],[0, 0.19, 0],[0, 0, 0.19]])
		mult  = 2
		Delta4R = valley.valley('Delta4R', L, W, mstar, mult)

		# Create Delta2 Valley
		mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
		mult = 2
		Delta2 = valley.valley('Delta2', L, W, mstar, mult)
	else:
		# Error, unknow sidewall
		print 'Unknow sidewall ', sidewall
		exit()
	
	
	
	
	
	
	# grid
	nx = 50
	ny = 50
	n_eigen = 15
	fd_mesh = fdmesh.FDMesh(L, W, nx, ny)


	# Create Schroedinger-Poisson Solver
	# Include all valleys
	print 'Constructing solver ...'
	sp = sp2d.SP2D([Delta4L, Delta4R, Delta2], fd_mesh)
	#sp = sp2d.SP2D([Delta2], fd_mesh)

	# Do self-consistent calculation
	# Compute self-consistent Fermi level






	# Create a CB profile
	peak = 0.300
	#cbfun = lambda x, y: -peak*np.exp(-(x-2.5)*(x-2.5)/(2.0*2.0))*np.exp(-(y-3.0)*(y-3.0)/(5.0*5.0))
	#cbfun = lambda x, y: -peak*np.exp(-(x-2.5)*(x-2.5)/(3.25*3.25))
	cbfun = lambda x, y : x * (L-x) * peak * 4.0/(L*L)
	#cbfun = lambda x, y : y * (L-y) * peak * 4.0/(L*L)
	curvature = 1e18 * peak * 2.0 * 4.0/(L*L)
	charge_density = curvature * phys.eps0 * phys.eps_Si * 1e-6 / phys.q
	total_charge = charge_density * L * W * 1e-7 * 1e-7
	surface_charge = 0.5 * total_charge / (W * 1e-7)

	# Now solve system
	print 'Starting eigensolve ...'
	sp.solveWithCB(cbfun, n_eigen=n_eigen)

	# RHS for charge self-consistency
	RHS = lambda x : sp.computeSidewallCharge(x)[0] - surface_charge
	# Initial guess
	x0 = 0.15

	print 'Starting Newton solve ...'
	xf= optimize.fsolve(RHS, x0, xtol=1e-5)
	print 'charge=', surface_charge, 'Fermi level : ', xf

	(tot_charge, v_charge) = sp.computeSidewallCharge(xf)

	set_plt_defaults()

	# Get DOS for each valley
	(energy, ret_val) = sp.computeDOS()
	D4_DOS = ret_val[0]
	D2_DOS = ret_val[2]
	plt.figure(1)
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

	
		
	# Build charge density
	ZCharge = sp.computeChargeDensity(xf)
			
	# Now plot the accumulated charge density
	plt.figure()
	levels = np.linspace(np.amin(ZCharge), np.amax(ZCharge), 30)
	cp = plt.contourf(X, Y, ZCharge, levels=levels)
	plt.colorbar(cp)
	plt.xlabel('X [nm]')
	plt.ylabel('Y [nm]')		
		
   

	

	

plt.show()

exit()
