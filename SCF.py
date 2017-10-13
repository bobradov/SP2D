
import numpy as np
import scipy.sparse.linalg as la_sparse
import SP2D as sp2d
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
from   SuperLattice import SuperLattice

class SCF:
	def __init__(self, schr, poiss, super_lattice='none'):
		self.schr  = schr
		self.poiss = poiss
		self.super_lattice = super_lattice

	def solve(self, bc_val_dict, max_iter=200, tol=1e-4, alpha=0.1,
								 Ef=0.0, n_eigen=20):


		#+---------------------------------------------------------+
		#|
		#| Initial Solver: Poisson with initial guess for charge
		#|                 Followed by QM-soltion with guessed potential
		#|
		#+---------------------------------------------------------+


		# Initial guess for charge is a volume-inverted Gaussian
		peak = 4e19
		L  = self.schr.fdmesh.L
		W  = self.schr.fdmesh.W
		nx = self.schr.fdmesh.nx
		ny = self.schr.fdmesh.ny
		xc = L/2.0
		yc = W/2.0
		rhs_initial = self.poiss.buildRHS(
			bc_val_dict,
			charge_func=lambda x, y : peak * np.exp(
			                        -(x-xc)*(x-xc)/(L*L)
			                        -(y-yc)*(y-yc)/(W*W))
		)


		# Solve the Poisson system
		initial_potential = la_sparse.spsolve(self.poiss.Lap_mat, rhs_initial)
		# Push solution to grid
		(X, Y, Z_initial_pot) = self.poiss.pmesh.mapSolutionToGrid(
													initial_potential
													)
		# Initial potential array
		prev_potential = initial_potential

		# Initial Schroedinger solution
		self.poiss.pmesh.loadAccumulator(initial_potential)

		# Use a conduction band which is simply the negative of the potenatial
		cbfun = lambda x, y: -1.0*self.poiss.pmesh.nn_interp(x, y)

		# Now solve eiogensystem
		self.schr.solveWithCB(cbfun, n_eigen=n_eigen)

		# Build QM charge density
		wfn_dict          = self.schr.getSolution(0)
		(Xfd,Yfd,Z_eigen) = wfn_dict['wfns']
		ZCharge_init      = self.schr.computeChargeDensity(Ef)
		initial_charge    = np.reshape(ZCharge_init, nx*ny)



		#+------------------------------------------------------------+
		#|
		#| Self-consistency SP Loop
		#|
		#+------------------------------------------------------------+
		prev_charge    = initial_charge
		prev_potential = initial_potential
		converged      = False

		print('+--------------+------------------------+')
		print('|    Iter      |    Potential Norm      |')
		print('+--------------+------------------------+')

		for cur_iter in range(1, max_iter):
			# Prepare for Poisson solve
			# Load previous Schroediner solution into accumulator
			# Get estimated new potential using full previous charge
			# Then "damp" by taking only a fraction of
			# the potential as the update
			self.schr.fdmesh.loadAccumulator(prev_charge)
			rhs = self.poiss.buildRHS(
				bc_val_dict,
				charge_grid=self.schr.fdmesh
			)
			pred_potential = la_sparse.spsolve(self.poiss.Lap_mat, rhs)

			# Compute update error
			delta_pot =  pred_potential-prev_potential
			norm = np.sqrt((delta_pot*delta_pot).sum())

			if cur_iter % 5 == 0 or cur_iter == 1:
				#print ' >>> Iter: ', cur_iter, ' Norm: ', norm
				print('| {0:5d}        |    {1:8.4f}            |'.format(cur_iter, norm))
			# Now do Schroedinger solve with damped new potential
			new_potential  = alpha*pred_potential + (1.0-alpha)*prev_potential

			# Save potential for next iter
			prev_potential = new_potential

			self.poiss.pmesh.loadAccumulator(new_potential)

			# Construct conduction band:
			# Start with negative of electrostatic potential
			cbfun   = lambda x, y: -1.0*self.poiss.pmesh.nn_interp(x, y)

			# Add a superlattice contribution to the conduction band
			if self.super_lattice != 'none':
				sl_fun  = lambda x, y: self.super_lattice.evaluate(y)
				tot_fun = lambda x, y: cbfun(x, y) + sl_fun(x, y)
			else:
				tot_fun = cbfun

			# Solve eigenproblem with completed conduction band
			self.schr.solveWithCB(tot_fun, n_eigen=n_eigen)

			wfn_dict          = self.schr.getSolution(0)
			(Xfd,Yfd,Z_eigen) = wfn_dict['wfns']
			ZCharge           = self.schr.computeChargeDensity(Ef)
			new_charge        = np.reshape(ZCharge, nx*ny) # Initial charge array

			prev_charge = new_charge

			if norm <= tol:
				converged = True
				break

		print('+--------------+------------------------+')

		# Save potential from final iteration of Poisson
		self.final_potential = new_potential

		return (cur_iter, converged)
