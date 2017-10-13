#! /opt/local/bin/python2.7

import numpy as np
import phys as phys
#from scipy.sparse import dok_matrix
import scipy.sparse as sparse
import scipy.sparse.linalg as la_sparse
import scipy.linalg        as la

n = 10

S = sparse.dok_matrix((n,n))
for i in range(0,n):
    S[i,i] = 2.0
    if i>0:
        S[i,i-1] = -1.0
    if i<n-1:
        S[i,i+1] = -1.0

#print S

csr_mat = S.tocsr()

print 'In CSR format:'
print csr_mat.toarray()

# Do an eigensolve
vals, vecs = la_sparse.eigsh(csr_mat, k=3, sigma=0.1)
print 'Eigenvalues: ', vals

dvals, dvecs = la.eig(csr_mat.toarray())
print 'Dense solve:', dvals
