import numpy as np


class FDMesh:

    def __init__(self, L, W, nx, ny):
        self.L = L
        self.W = W
        self.nx = nx
        self.ny = ny
        self.dx = L/(nx-1.0)
        self.dy = W/(ny-1.0)


	# Create "plottable" arrays of a set
	# of arrays on meshes.
	# For example, a set of eigenvectors

    def map2D(self, vecs):
        # Create the X,Y arrays
        X = np.linspace(0.0, self.L, self.nx)
        Y = np.linspace(0.0, self.W, self.ny)
        X, Y = np.meshgrid(X, Y)

        # Now create the "Z" arrays, i.e. wfns
        # Do just the ground state first

        (rows, num_eigenvecs) = vecs.shape
        #print 'rows:', rows, ' num_eigenvecs:', num_eigenvecs
        #print 'nx=', self.nx, ' ny=', self.ny
        #print vecs
        Z = np.zeros((self.ny, self.nx, num_eigenvecs))
        for cur_col in range(0,num_eigenvecs):
          cur_wf = vecs[:,cur_col]
          cur_wf_vec = cur_wf[:]
          Z[:,:,cur_col] = np.reshape(cur_wf_vec,
                                     (self.ny, self.nx),
                                     order='C')

        return (X, Y, Z)

    # mesh_func is a lambda(x,y) which is evaluated
    # on the mesh, the applied to the diagonal of a matrix
    def mapToGlobaleMat(self, mesh_func):
        X = np.linspace(0.0, self.L, self.nx)
        Y = np.linspace(0.0, self.W, self.ny)
        X, Y = np.meshgrid(X, Y)
        Z = mesh_func(X,Y)
        #print Z
        return np.reshape(Z, (self.nx * self.ny), order='C')
        
    def loadAccumulator(self, z):
    	self.acc = z
    	
    def	nn_interp(self, x, y):
    	# Find indices nearest to x, y
    	# Evaluate from array in accumulator
    	j = int(round(x/self.dx))
    	i = int(round(y/self.dy))
    	return self.acc[self.nx*i+j]
    	
    	
