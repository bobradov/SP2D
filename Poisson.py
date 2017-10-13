
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as la_sparse

import phys as phys
import PoissonMesh as pmesh


class Poisson:
    def __init__(self, p_mesh, bc_dict):
        self.pmesh = p_mesh
        self.bc_dict = bc_dict

        self.Lap_mat = self.buildLaplaceMat()
        # Make sure Lap-matrix is Hermitian
        assert (z in self.Lap_mat for z in np.transpose(self.Lap_mat))
	'''
    def solve(self, n_eigen=4, CB='none'):
        # Kinetic energy operator is fixed and assembled
        # Build Hamiltonian by adding potential energy
        H = self.buildHamiltonian(CB)

        # Now solve eigenproblem
        vals, vecs = la_sparse.eigsh(H, k=n_eigen, sigma=0.1)

        # eigenvalues are energies relative to
        # CB minimum
        # Eigenvectors are in global matrix space
        # Must be mapped onto 2-D grid to obtain wfns
        wfns = self.fdmesh.map2D(vecs)

        return (vals / self.energy_scale, wfns)
	'''

    def mapper(self, i, j):
		return self.pmesh.nx * i + j

    def buildRHS(self, bc_val_dict, charge='none'):
		rows = self.pmesh.nx * self.pmesh.ny
		rhs  = np.zeros(rows)
		
		# Add Dirichlet on left
		if self.bc_dict['left']=='Dirichlet':
			j = 0
			for i in range(0, self.pmesh.ny):
				row   = self.mapper(i, j)
				rhs[row] = bc_val_dict['left']
			
		# Add Dirichlet on right
		if self.bc_dict['right']=='Dirichlet':
			j = self.pmesh.nx-1
			for i in range(0, self.pmesh.ny):
				row   = self.mapper(i, j)
				rhs[row] = bc_val_dict['right']
		
		# Add Dirichlet on bot
		if self.bc_dict['bot']=='Dirichlet':
			i = 0
			for j in range(0, self.pmesh.nx):
				row   = self.mapper(i, j)
				rhs[row] = bc_val_dict['bot']
			
		# Add Dirichlet on top,
		if self.bc_dict['top']=='Dirichlet':
			i = self.pmesh.ny-1
			for j in range(0, self.pmesh.nx):
				row   = self.mapper(i, j)
				rhs[row] = bc_val_dict['top']		
		
			
		if charge !='none':
			# Use charge function to add charge contribution
			# Add charge only to semiconductor nodes
			# Assume that charge density is zero elsewhere
			(ixl, ixh) = self.pmesh.getXInterfaces()
			(iyl, iyh) = self.pmesh.getYInterfaces()
	
			dx = self.pmesh.schr_mesh.dx
			dy = self.pmesh.schr_mesh.dy
	
			i_start = iyl
			i_end   = iyh+1
			j_start = ixl
			j_end   = ixh+1
			
			#print 'dx=', dx, ' j_start=', j_start, ' j_end=', j_end
			#print 'dy=', dx, ' i_start=', i_start, ' i_end=', i_end
			
			for i in range(i_start, i_end):
				y = self.pmesh.schr_mesh.W-(i-i_start)*dy
				#print 'y=', y
				for j in range(j_start, j_end):
					x = (j-j_start)*dx
					row   = self.mapper(i, j)
					val   = -1.0*100*charge(x, y)*dx*dy*1e-14*phys.q/phys.eps0
					#print 'x=', x, ' val=', val
					rhs[row] = val
			
			
		return rhs
		

    def buildLaplaceMat(self):
        rows = self.pmesh.nx * self.pmesh.ny
        cols = rows
        Lap_mat_dok = sparse.dok_matrix((rows, rows))
      
      	print 'Laplace rows: ', rows
      	print 'Laplace nx:   ', self.pmesh.nx
      	print 'Laplace ny:   ', self.pmesh.ny
      
        # Assemble Laplace operator
        # Loop over nodes, add coupling
        # to neighbors
        
        
        
        (ixl, ixh) = self.pmesh.getXInterfaces()
        (iyl, iyh) = self.pmesh.getYInterfaces()
        
        print 'x-Interior nodes: ', (ixl, ixh)
        print 'y-Interior nodes: ', (iyl, iyh)
        
        # Coupling constants for interior nodes
        tx = self.pmesh.eps_semi*self.pmesh.schr_mesh.dy/self.pmesh.schr_mesh.dx
        ty = self.pmesh.eps_semi*self.pmesh.schr_mesh.dx/self.pmesh.schr_mesh.dy
         
        # Coupling constants for interior dielectric nodes
        (ox_dx, ox_dy) = self.pmesh.getMeshSpacing()
        ox_tx = self.pmesh.eps_ox * ox_dy / ox_dx
        ox_ty = self.pmesh.eps_ox * ox_dx / ox_dy
        
       
        
        #-----------------------------------------------
        #
        # Start with semiconductor interior nodes
        # These are nodes at least one meshline away
        # from interface nodes
        #
        #-----------------------------------------------
        
        for i in range(iyl+1, iyh-1+1):
            for j in range(ixl+1, ixh-1+1):
                # on-site
                row   = self.mapper(i, j)
                col   = row
                Lap_mat_dok[row,col] = -2*tx - 2*ty

                # Coupled right
                col_r   = self.mapper(i, j+1)
                Lap_mat_dok[row, col_r] = tx

                # Coupled left
                col_l   = self.mapper(i, j-1)
                Lap_mat_dok[row, col_l] = tx

                # Coupled up
                col_u   = self.mapper(i-1, j)
                Lap_mat_dok[row, col_u] = ty

                # Couple down
                col_d   = self.mapper(i+1, j)
                Lap_mat_dok[row, col_d] = ty

        #-------------------------------------------------------
        # Handle Interface nodes
        # These are the nodes at the semiconductor-dielectric
        # interface
        # Coupling at interface nodes more complicated
        # Construct on the fly using the following helper
        # variables:
        eps   = self.pmesh.eps_semi
        epsox = self.pmesh.eps_ox
        dx    = self.pmesh.schr_mesh.dx
        dy    = self.pmesh.schr_mesh.dy 
        #-------------------------------------------------------
        
        # Left boundary
        
        '''
        # Top left corner
        row = self.mapper(iyl, ixl)
        col = row
        Lap_mat_dok[row, col] = 0.0
        # right
        cr = self.mapper(iyl, ixl+1)
        Lap_mat_dok[row, cr] = 
        '''
        
        for i in range(iyl, iyh+1):
			j = ixl

			# Resolve material properties
			# on both sides of the boundary

			# on-site
			row   = self.mapper(i, j)
			col   = row
			Lap_mat_dok[row, col] = 0.0

			# Coupled right - right material properties
			col_r   = self.mapper(i, j+1)
			Lap_mat_dok[row, col_r] = eps*dy/dx
			Lap_mat_dok[row, col]  -= eps*dy/dx

			# Coupled left - left material properties
			col_l   = self.mapper(i, j-1)
			Lap_mat_dok[row, col_l] = epsox*dy/ox_dx
			Lap_mat_dok[row, col]  -= epsox*dy/ox_dx

			# Coupled up - right and left material props
			# Mesh density is same on both sides
			col_u   = self.mapper(i-1, j)
			Lap_mat_dok[row, col_u] = 0.5*eps*dx/dy + 0.5*epsox*ox_dx/dy
			Lap_mat_dok[row, col]  -= 0.5*eps*dx/dy + 0.5*epsox*ox_dx/dy

			# Couple down - right and left material props
			# Mesh density is the same on both sides
			col_d   = self.mapper(i+1, j)
			Lap_mat_dok[row, col_d] = 0.5*eps*dx/dy + 0.5*epsox*ox_dx/dy
			Lap_mat_dok[row, col]  -= 0.5*eps*dx/dy + 0.5*epsox*ox_dx/dy     

		# Right boundary
        
        for i in range(iyl, iyh+1):
			j = ixh

			# Resolve material properties
			# on both sides of the boundary

			# on-site
			row   = self.mapper(i, j)
			col   = row
			Lap_mat_dok[row,col] = 0.0

			# Coupled right - right material properties
			col_r   = self.mapper(i, j+1)
			Lap_mat_dok[row, col_r] = epsox*dy/ox_dx
			Lap_mat_dok[row, col]  -= epsox*dy/ox_dx

			# Coupled left - left material properties
			col_l   = self.mapper(i, j-1)
			Lap_mat_dok[row, col_l] = eps*dy/dx
			Lap_mat_dok[row, col]  -= eps*dy/dx

			# Coupled up - left and right material props
			col_u   = self.mapper(i-1, j)
			Lap_mat_dok[row, col_u] = 0.5*eps*dx/dy + 0.5*epsox*ox_dx/dy
			Lap_mat_dok[row, col]  -= 0.5*eps*dx/dy + 0.5*epsox*ox_dx/dy

			# Couple down - left and right material props
			col_d   = self.mapper(i+1, j)
			Lap_mat_dok[row, col_d] = 0.5*eps*dx/dy + 0.5*epsox*ox_dx/dy
			Lap_mat_dok[row, col]  -= 0.5*eps*dx/dy + 0.5*epsox*ox_dx/dy
			
		# Top boundary
        
        for j in range(ixl, ixh+1):
			i = iyl

			# Resolve material properties
			# on both sides of the boundary

			# on-site
			row   = self.mapper(i, j)
			col   = row
			Lap_mat_dok[row,col] = 0.0

			# Coupled right - combined material properties
			col_r   = self.mapper(i, j+1)
			Lap_mat_dok[row, col_r] = 0.5*eps*dy/dx + 0.5*epsox*ox_dy/dx
			Lap_mat_dok[row, col]  -= 0.5*eps*dy/dx + 0.5*epsox*ox_dy/dx

			# Coupled left - left material properties
			col_l   = self.mapper(i, j-1)
			Lap_mat_dok[row, col_l] = 0.5*eps*dy/dx + 0.5*epsox*ox_dy/dx
			Lap_mat_dok[row, col]  -= 0.5*eps*dy/dx + 0.5*epsox*ox_dy/dx

			# Coupled up - top material properties
			col_u   = self.mapper(i-1, j)
			Lap_mat_dok[row, col_u] = epsox*dx/ox_dy
			Lap_mat_dok[row, col]  -= epsox*dx/ox_dy

			# Couple down - right material properties
			col_d   = self.mapper(i+1, j)
			Lap_mat_dok[row, col_d] = eps*dx/dy	       	  
			Lap_mat_dok[row, col]  -= eps*dx/dy

		# Bottom boundary
        
        for j in range(ixl, ixh+1):
			i = iyh

			# Resolve material properties
			# on both sides of the boundary

			# on-site
			row   = self.mapper(i, j)
			col   = row
			Lap_mat_dok[row,col] = 0.0

			# Coupled right - combined material properties
			col_r   = self.mapper(i, j+1)
			Lap_mat_dok[row, col_r] = 0.5*eps*dy/dx + 0.5*epsox*ox_dy/dx
			Lap_mat_dok[row, col]  -= 0.5*eps*dy/dx + 0.5*epsox*ox_dy/dx

			# Coupled left - left material properties
			col_l   = self.mapper(i, j-1)
			Lap_mat_dok[row, col_l] = 0.5*eps*dy/dx + 0.5*epsox*ox_dy/dx
			Lap_mat_dok[row, col]  -= 0.5*eps*dy/dx + 0.5*epsox*ox_dy/dx

			# Coupled up - top material properties
			col_u   = self.mapper(i-1, j)
			Lap_mat_dok[row, col_u] = eps*dx/dy
			Lap_mat_dok[row, col]  -= eps*dx/dy

			# Couple down - right material properties
			col_d   = self.mapper(i+1, j)
			Lap_mat_dok[row, col_d] = epsox*dx/ox_dy
			Lap_mat_dok[row, col]  -= epsox*dx/ox_dy
		
		#---------------------------------------------------------	
		# Interior elements of insulator	       	 
		# Loop over all non-boundary nodes, skipping everything
		# that touches the semiconductor
		# Must also skip all Dirichlet boundaty nodes
		#---------------------------------------------------------
		
		# Determine boundary nodes
        if self.bc_dict['bot']=='Dirichlet':
			start_i = 1
        else:
			start_i = 0
			
        if self.bc_dict['top']=='Dirichlet':
			end_i = self.pmesh.ny-1
        else:
			end_i = self.pmesh.ny
			
        if self.bc_dict['left']=='Dirichlet':
			start_j = 1
        else:
			start_j = 0
			
        if self.bc_dict['right']=='Dirichlet':
			end_j = self.pmesh.nx-1
        else:
			end_j = self.pmesh.nx
		
        for i in range(start_i, end_i):
            for j in range(start_j, end_j):

				# Skip if semiconductor node
				if i>=ixl and i<=ixh and j>=iyl and j<=iyh:
					continue
	
				# Interior dielectric node
				# Need to overwrite boundary nodes

				# on-site
				row   = self.mapper(i, j)
				col   = row
				Lap_mat_dok[row,col] = 0.0

				# Coupled right
				if j<self.pmesh.nx-1:
					col_r   = self.mapper(i, j+1)
					Lap_mat_dok[row, col_r] = ox_tx
					Lap_mat_dok[row, col]  -= ox_tx

				# Coupled left
				if j>0:
					col_l   = self.mapper(i, j-1)
					Lap_mat_dok[row, col_l] = ox_tx
					Lap_mat_dok[row, col]  -= ox_tx

				# Couple up
				if i>0:
					col_u   = self.mapper(i-1, j)
					Lap_mat_dok[row, col_u] = ox_ty
					Lap_mat_dok[row, col]  -= ox_ty

				# Couple down
				if i<self.pmesh.ny-1:
					col_d   = self.mapper(i+1, j)
					Lap_mat_dok[row, col_d] = ox_ty
					Lap_mat_dok[row, col]  -= ox_ty
			
		#-----------------------------------------------
		#
		# Add boundary conditions
		# Not adding anything gives a Neumann condition
		#
		#-----------------------------------------------
        
        # Minus sign on the "1" is necessary
        # due to the final sign flip at the end
        
        # Add Dirichlet on left
        if self.bc_dict['left']=='Dirichlet':
			j = 0
			for i in range(0, self.pmesh.ny):
				row   = self.mapper(i, j)
				Lap_mat_dok[row, row] = -1.0
			
		# Add Dirichlet on right
        if self.bc_dict['right']=='Dirichlet':
			j = self.pmesh.nx-1
			for i in range(0, self.pmesh.ny):
				row   = self.mapper(i, j)
				Lap_mat_dok[row, row] = -1.0
		
		# Add Dirichlet on bot
        if self.bc_dict['bot']=='Dirichlet':
			i = 0
			for j in range(0, self.pmesh.nx):
				row   = self.mapper(i, j)
				Lap_mat_dok[row, row] = -1.0
			
		# Add Dirichlet on top,
        if self.bc_dict['top']=='Dirichlet':
			i = self.pmesh.ny-1
			for j in range(0, self.pmesh.nx):
				row   = self.mapper(i, j)
				Lap_mat_dok[row, row] = -1.0		
        
        
			
			
        return -Lap_mat_dok.tocsr()
