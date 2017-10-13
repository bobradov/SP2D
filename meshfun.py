import numpy as np

class FDMesh:

    def __init__(self, L, W, nx, ny):
        self.L = L
        self.W = W
        self.nx = nx
        self.ny = ny

    def map2D(self, vecs):
          # Create the X,Y arrays
          X = np.linspace(0.0, self.L, self.nx)
          Y = np.linspace(0.0, self.W, self.ny)
          X, Y = np.meshgrid(X, Y)

          # Now create the "Z" arrays, i.e. wfns
          # Do just the ground state first

          (rows, num_eigenvecs) = vecs.shape
          #print 'rows:', rows, ' num_eigenvecs:', num_eigenvecs
          Z = np.zeros((ny, nx, num_eigenvecs))
          for cur_col in range(0,num_eigenvecs):
              cur_wf = vecs[:,cur_col]
              cur_wf_vec = cur_wf[:]
              Z[:,:,cur_col] = np.reshape(cur_wf_vec,
                                         (self.ny, self.nx),
                                         order='C')
