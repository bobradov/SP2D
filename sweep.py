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

# Make some evil plots
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
	


lines = ['-sr','--sb','-Dk','--Dm',':ok']
sizes = [10.0, 8.0, 6.0, 4.0]
#sizes = [20.0]

sidewall='110'

for plot_ind, curW in enumerate(sizes):
	# Generate DOS for a single valley
	# This is an example of a Delta4 valley in (110) Si
	#n_points = 400
	L = curW
	W = 30.0    
	if sidewall=='110':
		# Create Delta4 Valleys
		mstar = np.array([[0.315, 0.478, 0],[0.478, 0.19, 0],[0, 0, 0.315]])
		mult  = 4
		Delta4 = valley.valley('Delta4', L, W, mstar, mult)

		# Create Delta2 Valley
		mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
		mult = 2
		Delta2 = valley.valley('Delta2', L, W, mstar, mult)

	elif sidewall=='100':
		# Create Delta4 Valleys
		mstar = np.array([[0.19, 0, 0],[0, 0.19, 0],[0, 0, 0.92]])
		mult  = 4
		Delta4 = valley.valley('Delta4', L, W, mstar, mult)

		# Create Delta2 Valley
		mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
		mult = 2
		Delta2 = valley.valley('Delta2', L, W, mstar, mult)
	else:
		# Error, unknow sidewall
		print 'Unknow sidewall ', sidewall
		exit()
    
    
    
    
	# grid
	nx = 30
	ny = 60
	fd_mesh = fdmesh.FDMesh(L, W, nx, ny)


	# Create Schroedinger-Poisson Solver
	# Include all valleys
	print 'Constructing solver ...'
	sp = sp2d.SP2D([Delta4, Delta2], fd_mesh)

	# Do self-consistent calculation
	# Compute self-consistent Fermi level

	n_sweeps = 20
	#peakvals = np.linspace(0.005, 0.25, n_sweeps)
	peakvals = np.logspace(-2, -0.5, n_sweeps)

	charge_array  = np.zeros(n_sweeps)
	surf_array    = np.zeros(n_sweeps)
	fermi_array   = np.zeros(n_sweeps)
	valley_charge = np.zeros((2,n_sweeps))
	occ_ratio     = np.zeros(n_sweeps)
	vinj          = np.zeros(n_sweeps)
	vel           = np.zeros(n_sweeps)


	D4_gs         = np.zeros(n_sweeps)
	D2_gs         = np.zeros(n_sweeps)
	D4_fe         = np.zeros(n_sweeps)
	D2_fe         = np.zeros(n_sweeps)

	print 'Sweeping ...'

	for index, cur_peak in enumerate(peakvals):
		# Create a CB profile
		#cbfun = lambda x, y: -cur_peak*np.exp(-(x-2.5)*(x-2.5)/(2.0*2.0))*np.exp(-(y-curW/2.5)*(y-curW/2.5)/(curW*curW))
		cbfun = lambda x, y : x * (L-x) * cur_peak * 4.0/(L*L)
		#cbfun = lambda x, y : y * (L-y) * cur_peak * 4.0/(L*L)
		#cbfun = lambda x, y : x * y * (L-y) * (L-x) * cur_peak * 4.0/(L*L)
		curvature = 1e18 * cur_peak * 2.0 * 4.0/(L*L)
		charge_density = curvature * phys.eps0 * phys.eps_Si * 1e-6 / phys.q
		total_charge = charge_density * L * W * 1e-7 * 1e-7
		surface_charge = 0.5 * total_charge / (W * 1e-7)

		# Now solve system
		print 'Starting eigensolve ...'
		sp.solveWithCB(cbfun, n_eigen=20)

		RHS = lambda x : sp.computeSidewallCharge(x)[0] - surface_charge
		x0 = 0.15

		print 'Starting Newton solve ...'
		xf= optimize.fsolve(RHS, x0, xtol=1e-5)
		print 'charge=', surface_charge, 'Fermi level : ', xf

		(tot_charge, v_charge) = sp.computeSidewallCharge(xf)

		#print 'Shapes: ',tot_charge.shape, ' ', v_charge.shape

		charge_array[index] = total_charge
		surf_array[index]   = surface_charge
		fermi_array[index]  = xf
		valley_charge[:,index] = v_charge
		occ_ratio[index]       = v_charge[1] / tot_charge
		vinj[index]   = np.sqrt(0.315)*v_charge[1]/tot_charge + \
		                np.sqrt(0.19)*v_charge[0]/tot_charge
		                
		vel[index] = sp.computeVelocity(xf) / total_charge
	  

		# Ground states
		D4_gs[index] = sp.getSolution(0)['evals'][0]
		D2_gs[index] = sp.getSolution(1)['evals'][0]
		# First excited states
		D4_fe[index] = sp.getSolution(0)['evals'][1]
		D2_fe[index] = sp.getSolution(1)['evals'][1]

	set_plt_defaults()

	'''    
	plt.figure()
	plt.plot(surf_array, fermi_array, '-k', linewidth=3)
	plt.plot(surf_array, D4_gs, '-r', linewidth=3)
	plt.plot(surf_array, D2_gs, '-b', linewidth=3)
	plt.plot(surf_array, D4_fe, '--r', linewidth=3)
	plt.plot(surf_array, D2_fe, '--b', linewidth=3)
	plt.grid(True)
	plt.xlabel('Surface Sheet Density')
	plt.ylabel('Energy Levels')


	plt.figure()
	plt.plot(surf_array, valley_charge[0,:], '-r', linewidth=3)
	plt.plot(surf_array, valley_charge[1,:], '-b', linewidth=3)
	plt.grid(True)
	plt.xlabel('Surface Sheet Density')
	plt.ylabel('Valley Sheet Density')
	'''    


	plt.figure(1)
	plt.semilogx(surf_array, occ_ratio, lines[plot_ind], linewidth=3)
	plt.grid(True)
	plt.xlabel('Surface Sheet Density')
	plt.ylabel('Fraction in D2 valleys')

	plt.figure(2)
	plt.semilogx(surf_array, vinj, lines[plot_ind], linewidth=3)
	plt.grid(True)
	plt.xlabel('Surface Sheet Density')
	plt.ylabel('Relative injection velocity')


	plt.figure(3)
	plt.semilogx(surf_array, vel, lines[plot_ind], linewidth=3)
	plt.grid(True)
	plt.xlabel('Surface Sheet Density')
	plt.ylabel('Absolute injection velocity')

plt.show()

exit()
