
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as la_sparse

import phys as phys
import valley as valley


class Schroedinger:
    def __init__(self, valley, fdmesh):
        self.valley = valley
        self.fdmesh = fdmesh
        self.energy_scale = phys.q * 2.0 * phys.m0 * \
            phys.nm_to_m * phys.nm_to_m / (phys.hbar * phys.hbar)

        self.Tmat = self.buildKineticEnergyMat()
        # Make sure T-matrix is Hermitian
        assert (z in self.Tmat for z in np.transpose(self.Tmat))

    def solve(self, n_eigen=4, CB='none', sigma=0.0):
        # Kinetic energy operator is fixed and assembled
        # Build Hamiltonian by adding potential energy
        H = self.buildHamiltonian(CB)

        # Now solve eigenproblem
        vals, vecs = la_sparse.eigsh(H, k=n_eigen, sigma=sigma)

        # eigenvalues are energies relative to
        # CB minimum
        # Eigenvectors are in global matrix space
        # Must be mapped onto 2-D grid to obtain wfns
        wfns = self.fdmesh.map2D(vecs)

        return (vals / self.energy_scale + self.valley.shift, wfns)


    # Adds the potential energy to the diagonal
    # of the Kinetic Energy to obtain the total
    # Hamiltonian
    def buildHamiltonian(self, CB):
        if CB == 'none':
           return self.Tmat
		# A conduction-band lambda has been passed into the function
		# Need to map it only Hamiltonian matrix
        '''
        Vmat = self.fdmesh.mapToGlobaleMat(CB)
        #print self.energy_scale * np.diag(Vmat)
        return self.Tmat + self.energy_scale * np.diag(Vmat)
        '''
		# Loop over nodes and add contribution to Hamiltonian
		# Store values in a 1-D array, then add to diagonal
        Varray = np.zeros(self.fdmesh.nx * self.fdmesh.ny)
        mapper = lambda i,j : self.fdmesh.nx * i + j

        for i in range(0, self.fdmesh.ny):
            y = self.fdmesh.dy * i
            for j in range(0, self.fdmesh.nx):
                gb_index = mapper(i, j)
                x = self.fdmesh.dx * j
                Varray[gb_index] = CB(x, y)

        return self.Tmat + self.energy_scale * np.diag(Varray)



    def buildKineticEnergyMat(self):
        rows = self.fdmesh.nx * self.fdmesh.ny
        cols = rows
        Tmat_dok = sparse.dok_matrix((rows, rows))
        #dx   = self.fdmesh.L / (self.fdmesh.nx-1)
        #dy   = self.fdmesh.W / (self.fdmesh.ny-1)
        dx   = self.fdmesh.dx
        dy   = self.fdmesh.dy
        mxx  = self.valley.mstar[0,0]
        myy  = self.valley.mstar[1,1]
        mxy  = self.valley.mstar[0,1]
        tx   = 1.0/(mxx * dx * dx)
        ty   = 1.0/(myy * dy * dy)


        # Assemble kinetic energy operator
        # The mapper lambda computes indices
        # of the global matrix based on mesh
        # indices (i,j)
        mapper = lambda i,j : self.fdmesh.nx * i + j

        # Loop over nodes, add coupling
        # to neighbors
        for i in range(0, self.fdmesh.ny):
            for j in range(0, self.fdmesh.nx):
                # on-site
                row   = mapper(i, j)
                col   = row
                Tmat_dok[row,col] = -2*tx -2*ty

                # Coupled right
                if j<self.fdmesh.nx-1:
                    col_r   = mapper(i, j+1)
                    Tmat_dok[row, col_r] = tx

                # Coupled left
                if j>0:
                    col_l   = mapper(i, j-1)
                    Tmat_dok[row, col_l] = tx

                # Couple up
                if i>0:
                    col_u   = mapper(i-1, j)
                    Tmat_dok[row, col_u] = ty

                # Couple down
                if i<self.fdmesh.ny-1:
                    col_d   = mapper(i+1, j)
                    Tmat_dok[row, col_d] = ty

                # Mixed partial derivatives
                if mxy != 0.0:
                    txy  = 1.0/(2.0 * mxy * dx * dy)
                    # Upper left
                    if j>0 and i>0:
                        col_ul = mapper(i-1, j-1)
                        Tmat_dok[row, col_ul] = txy
                    # Lower right
                    if j<self.fdmesh.nx-1 and i<self.fdmesh.ny-1:
                        col_lr = mapper(i+1, j+1)
                        Tmat_dok[row, col_lr] = txy
                    # Upper right
                    if i>0 and j<self.fdmesh.nx-1:
                        col_ur = mapper(i-1, j+1)
                        Tmat_dok[row, col_ur] = -txy
                    # Lower left
                    if i<self.fdmesh.ny-1 and j>0:
                        col_ll = mapper(i+1, j-1)
                        Tmat_dok[row, col_ll] = -txy

        return -Tmat_dok.tocsr()
