bias.py                                                                                             0000755 0002735 0004704 00000011171 13027250002 013111  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python



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
                                                                                                                                                                                                                                                                                                                                                                                                       charge.py                                                                                           0000755 0002735 0004704 00000004122 13025746215 013440  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python


import FDMesh as fdmesh
import PoissonMesh as pmesh
import Poisson     as pois
import numpy  as np
import phys

import scipy.sparse.linalg as la_sparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=1)
np.set_printoptions(linewidth=300)


L = 3.0
W = 3.0

# grid
nx = 10
ny = 10
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

vsp = 15.0
hsp = 6.0
#vsp = 3.0
#hsp = 3.0

p_mesh = pmesh.PoissonMesh(fdmesh, vsp, hsp, nh=5, nv=5,
							eps_semi=1, eps_ox=1)

# Find interface nodes
xi = p_mesh.getXInterfaces()
yi = p_mesh.getYInterfaces()
print 'xi: ', xi
print 'yi: ', yi

# Make a Poisson solver object

bc_dict = { 'top'   : 'Neumann', 
			'bot'   : 'Neumann', 
			'left'  : 'Dirichlet', 
			'right' : 'Dirichlet' }
			

'''
bc_dict = { 'top'   : 'Dirichlet', 
			'bot'   : 'Dirichlet', 
			'left'  : 'Neumann', 
			'right' : 'Neumann' }			
'''			
poisson = pois.Poisson(p_mesh, bc_dict)

#print 'Lap mat:'
pmat =  poisson.Lap_mat
#print pmat.toarray()

# Now build rhs mat
bc_val_dict = { 'left'  : 0.0, 
				'right' : 1.0,
				'top'   : 0.0,
				'bot'   : 1.0 }
				
peak = 1e20*0
xc = L/2.0
yc = W/2.0
rhs = poisson.buildRHS(bc_val_dict, charge = 
	lambda x, y : peak * np.exp(-(x-xc)*(x-xc)/(L*L) -(y-yc)*(y-yc)/(W*W))
	)
#print 'RHS: ', rhs
print 'RHS shape: ', rhs.shape
print 'Poisson shape: ', pmat.shape

# Solve the resulting system
z = la_sparse.spsolve(pmat, rhs)
#print 'x:', z

# Make some plots
xvals  = p_mesh.getXMesh()
yvals = p_mesh.getYMesh()
X, Y = np.meshgrid(xvals, yvals)
Z = z.reshape(p_mesh.ny, p_mesh.nx)

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z), np.amax(Z), 30)
cp = plt.contourf(X, Y, Z, levels=levels)
plt.contour(X, Y, Z, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	

ax.add_patch(
    patches.Rectangle(
        (0.0, 0.0),
        L,
        W,
        fill=False,      # remove background
        linewidth=2
    )
)


plt.show()	




                                                                                                                                                                                                                                                                                                                                                                                                                                              chargetest.py                                                                                       0000755 0002735 0004704 00000003140 13025274167 014341  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python


import FDMesh as fdmesh
import PoissonMesh as pmesh
import Poisson     as pois
import numpy  as np
import phys

import scipy.sparse.linalg as la_sparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=1)
np.set_printoptions(linewidth=300)


L = 5.0
W = 6.0

# grid
nx = 50
ny = 50
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

vsp = 2.0
hsp = 1.0

p_mesh = pmesh.PoissonMesh(fdmesh, vsp, hsp, nh=15, nv=15,
							eps_semi=11.9, eps_ox=2.2)

# Find interface nodes
xi = p_mesh.getXInterfaces()
yi = p_mesh.getYInterfaces()
print 'xi: ', xi
print 'yi: ', yi

# Make a Poisson solver object
poisson = pois.Poisson(p_mesh)

#print 'Lap mat:'
pmat =  poisson.Lap_mat
#print pmat.toarray()

# Now build rhs mat
rhs = poisson.buildRHS()
#print 'RHS: ', rhs
print 'RHS shape: ', rhs.shape
print 'Poisson shape: ', pmat.shape

# Solve the resulting system
z = la_sparse.spsolve(pmat, rhs)
#print 'x:', z

# Make some plots
xvals  = p_mesh.getXMesh()
yvals = p_mesh.getYMesh()
X, Y = np.meshgrid(xvals, yvals)
Z = z.reshape(p_mesh.ny, p_mesh.nx)

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z), np.amax(Z), 30)
cp = plt.contourf(X, Y, Z, levels=levels)
plt.contour(X, Y, Z, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	

ax.add_patch(
    patches.Rectangle(
        (0.0, 0.0),
        L,
        W,
        fill=False,      # remove background
        linewidth=2
    )
)


plt.show()	




                                                                                                                                                                                                                                                                                                                                                                                                                                CombinatorialLoop.py                                                                                0000755 0002735 0004704 00000011730 13035447135 015627  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python


# Class for combinatorial looping


class CombinatorialLoop(object):
    def __init__(self, DictOfLists):
    	'Constructor for the iterator, takes a dictionary of lists as input'
    	# Create a ListOfLists which stores the possible values for each
    	# dimension. The iterator loops over the outer 
    	# product of all the dimensions.
    	# Heterogeneous list types are possible, no type
    	# checking is done.
    	#
    	# The .self value at each iteration is a list
    	# of elemements, equal in length to the size of 
    	# the ListOfLists. Each iteration produces a
    	# different element from the outer product.
    	#
    	
    	self.vars = DictOfLists.keys()
    	ListOfLists = []
    	for cur_var in self.vars:
    		ListOfLists.append(DictOfLists[cur_var])
    	
    	self.data   = ListOfLists
    	self.NLists = len( ListOfLists )
    	self.first  = True
    	
    	# How many items in each list?
    	self.ListLengths = []
    	for i in range(0,self.NLists):
    		self.ListLengths.append( len(ListOfLists[i]) )
    	
    	# Initialize counter for each list to 0	
        self.CurrentIndex = [0]*self.NLists
        
        # Initialize current state
        self.current = []
        for i in range(0,self.NLists):
        	self.current.append( ListOfLists[i][0] )
        
        # Find the final element of each list 
        self.high = []
        for i in range(0,self.NLists):
        	self.high.append( ListOfLists[i][ self.ListLengths[i]-1 ] )
        	
        	
        	
    
    def __iter__(self):
        'Returns itself as an iterator object'
        return self
        
        
        
        
	# The next function finds the next suitable n-tuple from the
	# outer product. 
	# The return value is an array (not a tuple) containing 
	# an element from the outer product.
	# The increment algorithm is a couunter, with the first list
	# representing the least-significant digit, the next list
	# represeting the next digit etc.
	# The 'counter' algorithm proceeds until all list elements are 
	# incremented to max values.

    def next(self):
        'Returns the next value until no more items remain'
        
        # Is this the very first use? If so, return
        # initial state.
        # Incrementing happens with next use of 'next'
        if self.first:
        	self.first = False
        	ret_dict = {}
        	for index, cur_var in enumerate(self.vars):
        		ret_dict[cur_var] = self.current[index]
        	return ret_dict
        
        # Increment lowest index; if already at max, reset 
        # and incement next index
        incremented = False
        
        for digit in range(0,self.NLists):
        	# Try to increment this index
        	self.CurrentIndex[digit] += 1
        	                
        	if self.CurrentIndex[digit] == self.ListLengths[digit]:
        		# Carry
        		# The current digit is not incremented
        		# Instead, it's reset to the zero index
        		# 'incremented' is not set to true, allowing
        		# the next digit to be incremented in the 
        		# next iteration of the loop
        		self.CurrentIndex[digit] = 0
        		self.current[digit] = \
        			self.data[digit][self.CurrentIndex[digit]]
        		continue
        	else:
        		# No carry, actually increment this digit
        		# Done with 'next', no more looping over digits
        		# Declare 'incremented' True and break out of loop
        		self.current[digit] = self.data[digit][self.CurrentIndex[digit]]
        		incremented  = True
        		break
        		
        # Termination condition: no digits could be incremented,
        # every digit is at max index
        # If termination condition is not met, return the current
        # state instead of raising exception	 
        if not incremented:
            raise StopIteration
        else:
            # Build a dictionary which contains the current
            # values of each variable
            ret_dict = {}
            for index, cur_var in enumerate(self.vars):
        		ret_dict[cur_var] = self.current[index]
            return ret_dict
            
            
    # The reset function is used for restarting the iterator
    # This is useful if one wishes to re-use the iterator, rather
    # than creating a new one.              
            
    def reset(self):
    	# Initialize counter for each list to 0	
        self.CurrentIndex = [0]*self.NLists
        
        # Initialize current state
        self.current = []
        for i in range(0,self.NLists):
        	self.current.append( ListOfLists[i][0] )
       	
       	# Restart count 	
       	self.first  = True
              
      
# Testsuite
      
if __name__ == "__main__":      
     # Do some testcases       
	counter = CombinatorialLoop( { 'alpha1' : ['a','b','c'],
	                               'int1'   : [1,2],
	                               'alpha2' : ['x','y','z'] } )             
		        
	for curVal in counter:
		print 'curVal=', curVal
                                        ctest2.py                                                                                           0000755 0002735 0004704 00000005003 13027243172 013406  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

# Solver Poisson
# Charge is assigned to the semiconductor mesh directly
from scipy import interpolate

import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys

import scipy.sparse.linalg as la_sparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=2)
np.set_printoptions(linewidth=300)


L = 5.0
W = 5.0
IL = 0.5

# grid
nx = 30
ny = 30
nh = 20
nv = 15
no = 5
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.0
EOT = 0.85
vsp = 2.0
hsp = IL + (EOT-IL)*Eps_HiK/3.9

pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)


mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			
mesh.assignContact( lambda x, y: x < -hsp + 1e-8,   'Left' )
mesh.assignContact( lambda x, y: x >  L+hsp-1e-8, 'Right')
					
#print(mesh)

#exit()

poisson = pois.FVMPoisson(mesh)

#print 'Lap mat:'
pmat =  poisson.Lap_mat
#print pmat.toarray()


# Now build rhs mat
bc_val_dict = { 'Left' : 0.7, 'Right' : 0.7 }
peak = 8e19
xc = L/2.0
yc = W/2.0
charge_func=lambda x, y : peak * np.exp(-(x-xc)*(x-xc)/(L*L) -(y-yc)*(y-yc)/(W*W))

# Mesh for the semiconductor
xsemi = np.linspace(0.0, L, nx)
ysemi = np.linspace(0.0, W, ny)
XSemi, YSemi = np.meshgrid(xsemi, ysemi)
csemi = charge_func(XSemi, YSemi)
print csemi.shape
csemi_array = np.reshape(csemi, nx*ny)

fdmesh.loadAccumulator(csemi_array)



rhs = poisson.buildRHS(
	bc_val_dict, 
	charge_grid=fdmesh
)
#print 'RHS: ', rhs
print 'RHS shape: ', rhs.shape
print 'Poisson shape: ', pmat.shape


# Solve the resulting system
z = la_sparse.spsolve(pmat, rhs)
#print 'z:', z

(X, Y, Z) = poisson.pmesh.mapSolutionToGrid(z)


# Make some plots


fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z), np.amax(Z), 30)
cp = plt.contourf(X, Y, Z, levels=levels)
plt.contour(X, Y, Z, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	


ax.add_patch(
    patches.Rectangle(
        (0.0, 0.0),
        L,
        W,
        fill=False,      # remove background
        linewidth=2
    )
)

ax.add_patch(
    patches.Rectangle(
        (-IL, -IL),
        L+2*IL,
        W+2*IL,
        fill=False,      # remove background
        linewidth=2
    )
)


plt.show()	




                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             ctest.py                                                                                            0000755 0002735 0004704 00000004176 13027030200 013320  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys

import scipy.sparse.linalg as la_sparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=2)
np.set_printoptions(linewidth=300)


L = 5.0
W = 5.0
IL = 0.5

# grid
nx = 30
ny = 30
nh = 20
nv = 15
no = 5
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.0
EOT = 0.85
vsp = 2.0
hsp = IL + (EOT-IL)*Eps_HiK/3.9

pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)


mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			
mesh.assignContact( lambda x, y: x < -hsp + 1e-8,   'Left' )
mesh.assignContact( lambda x, y: x >  L+hsp-1e-8, 'Right')
					
#print(mesh)

#exit()

poisson = pois.FVMPoisson(mesh)

#print 'Lap mat:'
pmat =  poisson.Lap_mat
#print pmat.toarray()


# Now build rhs mat
bc_val_dict = { 'Left' : 0.7, 'Right' : 0.7 }
peak = 8e19
xc = L/2.0
yc = W/2.0
rhs = poisson.buildRHS(
	bc_val_dict, 
	charge_func=lambda x, y : peak * np.exp(-(x-xc)*(x-xc)/(L*L) -(y-yc)*(y-yc)/(W*W))
)
#print 'RHS: ', rhs
print 'RHS shape: ', rhs.shape
print 'Poisson shape: ', pmat.shape


# Solve the resulting system
z = la_sparse.spsolve(pmat, rhs)
#print 'z:', z

(X, Y, Z) = poisson.pmesh.mapSolutionToGrid(z)


# Make some plots


fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z), np.amax(Z), 30)
cp = plt.contourf(X, Y, Z, levels=levels)
plt.contour(X, Y, Z, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	


ax.add_patch(
    patches.Rectangle(
        (0.0, 0.0),
        L,
        W,
        fill=False,      # remove background
        linewidth=2
    )
)

ax.add_patch(
    patches.Rectangle(
        (-IL, -IL),
        L+2*IL,
        W+2*IL,
        fill=False,      # remove background
        linewidth=2
    )
)


plt.show()	




                                                                                                                                                                                                                                                                                                                                                                                                  energy.py                                                                                           0000755 0002735 0004704 00000005650 13020153236 013475  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python



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
                                                                                        FDMesh.py                                                                                           0000644 0002735 0004704 00000003366 13027031116 013310  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               import numpy as np


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
    	
    	
                                                                                                                                                                                                                                                                          FOMtoPandas.py                                                                                      0000644 0002735 0004704 00000002435 13032034362 014313  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               import pandas as pnd
import numpy  as np
import glob
import os
import os.path



# Read a set of data files, collate into a single Pandas DataFrame.
# The data files are all expected to have the same set of variables.
# Return value: Pandas DataFrame.

def FOMtoPandas(root='Split_', separator='_', fom_name='FOM.dat'):
		
	# Helper function to determine if a string represents a number		
	def is_number(s):
		try:
		    float(s)
		    return True
		except ValueError:
		    return False	
				

	# 1. Find all directories named 'Split_*'
	split_dirs = glob.glob(root+'*')
	split_dirs.sort()
	
	
	# Save all information into a dictionary, convert to Pandas DF
	table_dict = {}

	for cur_split in split_dirs:
		# Get split number
		(dummy, split_num) = cur_split.split(separator)
	
		# Now parse the FOM file
		fom_file = cur_split + '/' + fom_name
		if not os.path.isfile(fom_file):
			continue 
		fom = open(fom_file, 'r')
		for cur_line in fom:
			(name, val) = cur_line.split()
			# Add to table
			# Is this the first time we see this variable?
			if not name in table_dict:
				table_dict[name] = []
			# Is it an expected to be non-numerical?
			if is_number(val):
				val = float(val)
			table_dict[name].append(val)
		
		fom.close()

	# Create a pandas table
	df = pnd.DataFrame(table_dict)
	return df
                                                                                                                                                                                                                                   FVMMesh.py                                                                                          0000644 0002735 0004704 00000014760 13027252243 013455  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               import numpy as np

class FVMVertex:
	def __init__(self, index, x, y, ind_i, ind_j, boundary=False):
		self.index = index
		self.x     = float(x)
		self.y     = float(y)
		self.ind_i = ind_i
		self.ind_j = ind_j
		self.boundary = boundary
		self.contact  = 'none'
		
	def __str__(self):
		tstr = '((' + str(self.index) +  ') '
		tstr += '(' + str(self.x) + ', ' + str(self.y) + ') '
		tstr += str(self.boundary)
		if self.boundary:
			tstr += ' ' + self.contact
		return tstr
		


class FVMElement:
	def __init__(self, Index, Material, VertexList):
		self.index      = Index
		self.material   = Material
		self.vertexList = VertexList
		
	def __str__(self):
		tstr = 'Element: ' + str(self.index) + ' '
		tstr += str(self.material) + ' '
		tstr += str(self.vertexList)
		return tstr
		
	def getVertexIter(self, mesh):
		return ElemVertexIter(self, mesh)
		
#+-----------------------------------------------------------------------+
#|
#| Finite Volume Mesh Class                                        
#|
#| Supports cartesian meshes with arbitrary line distributions
#| Element properties are mapped onto the mesh via a zoner function
#| A dictionary linking material to properties is included in prop_dict
#|
#+-----------------------------------------------------------------------+		
		

class FVMMesh:
	def __init__(self, xlines, ylines, zoner_func, prop_dict):
		# Check to make sure there are no duplicate lines
		xdpl = self.duplicates(xlines)
		if xdpl !=-1:
			print 'Duplicated lines found in xlines:'
			print 'Duplicated line:', xlines[xdpl]
			raise ValueError('Duplicated xline')
			
		ydpl = self.duplicates(ylines)
		if ydpl !=-1:
			print 'Duplicated lines found in ylines:'
			print 'Duplicated line:', ylines[ydpl]
			raise ValueError('Duplicated yline')
			
		self.xlines = xlines   # x-coordinates of lines
		self.ylines = ylines   # y-coordinates of lines
		self.nxl    = len(xlines)
		self.nyl    = len(ylines)
		self.prop_dict = prop_dict
		self.nVert     = self.nxl * self.nyl
		self.nElem     = (self.nxl-1)*(self.nyl-1)
		self.zoner     = zoner_func
		
		# Allocate mesh arrays
		self.vertex_list = []
		self.elem_list   = []
		self.contacts    = {}
		
		self.buildElements()
		
	# Determine of any of the lines are duplicated
	# If so, return the index of the duplicated line
	# Its duplicate is at index+1
	# Assumes the input array is sorted
	def duplicates(self, lines):
		for index in range(0, len(lines)-1):
			if lines[index]==lines[index+1]:
				return index
		return -1
		
	def getElementIter(self):
		return ElementIter(self)
		
	def mapSolutionToGrid(self, z):
		X, Y = np.meshgrid(self.xlines, self.ylines)
		Z    = z.reshape(len(self.ylines), len(self.xlines))
		return (X, Y, Z)
		
	def mapFunctionToGrid(self, func):
		X, Y = np.meshgrid(self.xlines, self.ylines)
		Z    = func(X, Y)
		return (X, Y, Z)
		
	def mapper(self, i, j):
		return i*self.nxl + j
		
	
	def loadAccumulator(self, solution):
		self.solution = solution
		
	# Nearest-neighbor interpolation for x, y
	def nn_interp(self, x, y):
		# find the xline which most closely matches x
		index_x      = (np.abs(self.xlines-x)).argmin()
		index_y      = (np.abs(self.ylines-y)).argmin()
		global_index = self.mapper(index_y, index_x)
		
		return self.solution[global_index]
		
		
		
	def __str__(self):
		# Header
		tstr = 'FVM Mesh:\n'
		tstr += 'nxl: ' + str(self.nxl) + '\n'
		tstr += 'nyl: ' + str(self.nyl) + '\n'
		tstr += 'nElem: ' + str(self.nElem) + '\n'
		
		# All elements
		for cur_elem in self.elem_list:
			tstr += cur_elem.__str__() + '    '
			# Get vertex coordinates
			for curv_index in cur_elem.vertexList:
				tstr += self.vertex_list[curv_index].__str__() + ' '
			tstr += '\n'
			
		# Contacts
		for cur_con in self.contacts.keys():
			tstr += cur_con + ' : '
			for cur_ind in self.contacts[cur_con]:
				tstr += str(cur_ind) + ' '
			tstr += '\n'
 		return tstr 	
		
	def buildElements(self):
		# Loop over lines, create vertices, edges, elements
		# Vertical lines are left edge of element
		# 
		v_index = 0
		for ind_i, yi in enumerate(self.ylines):
			for ind_j, xi in enumerate(self.xlines):
				# Add a new vertex
				self.vertex_list.append(FVMVertex(v_index, 
												  xi, yi,
												  ind_i, ind_j))
				v_index += 1
				
		
		# Add elements
		# Loop over vertices, add element which has the vertex for
		# its top-left corner
		v_index    = 0
		elem_index = 0
		for yi in range(0, self.nyl-1):
			for xi in range(0, self.nxl-1):
				# Create nodelist
				# Local node numbering: counter-clockwise, start with bot left
				# 3 -- 2
				# |    |
				# 0 -- 1
				#
				node_list = [v_index, 
				         v_index+1, 
				         v_index+self.nxl+1,
				         v_index+self.nxl]
				         
				# Determine material based on zoner function
				# Sample at the center of the element
				
				#print 'doing v_index=', v_index, ' nxl=', self.nxl, ' nyl=', self.nyl
				
				mat = self.zoner( 0.5*(self.vertex_list[v_index].x +
				                       self.vertex_list[v_index+1].x),
				                  0.5*(self.vertex_list[v_index].y +
				                  	   self.vertex_list[v_index+self.nxl].y ) )
				                  	   
				self.elem_list.append( FVMElement( elem_index,
												   mat, 
												   node_list ) )
				v_index    += 1
				elem_index += 1
			# Completed row
			v_index += 1 # Move to next row
			
	def assignContact(self, cont_func, name):
		# Find all the nodes that belong to this contact
		self.contacts[name] = []
		
		for cur_vert in self.vertex_list:
			if cont_func(cur_vert.x, cur_vert.y):
				cur_vert.boundary = True
				cur_vert.contact  = name
				self.contacts[name].append(cur_vert.index)
					
#+---------------------------------------------------------------
#|
#| Iterators for FVMMesh
#|
#+---------------------------------------------------------------

# Return elements of mesh
# Returns mesh object, not index
class ElementIter:
	def __init__(self, mesh):
		self.mesh = mesh
		self.i    = 0
		
	def __iter__(self):
		return self
		
	def next(self):
		if self.i<self.mesh.nElem:
			self.i += 1
			return self.mesh.elem_list[self.i-1]
		else:
			raise StopIteration()
			
# Return vertices of an element
# Returns a tuple:
# (global_index_of_vertex, vertex_object)			
class ElemVertexIter:
	def __init__(self, elem, mesh):
		self.mesh = mesh
		self.elem = elem
		self.i    = 0
		
	def __iter__(self):
		return self
		
	def next(self):
		if self.i < 4:
			self.i += 1
			global_node = self.elem.vertexList[self.i-1]
			return (global_node, self.mesh.vertex_list[global_node])
		else:
			raise StopIteration()			
						
				
		
	
		
	
                FVMPoisson.py                                                                                       0000644 0002735 0004704 00000011520 13027027531 014202  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as la_sparse

import phys as phys
import FVMMesh as mesh


class FVMPoisson:
    def __init__(self, p_mesh):
        self.pmesh = p_mesh
     

        self.Lap_mat = self.buildLaplaceMat()
        # Make sure Lap-matrix is Hermitian
        assert (z in self.Lap_mat for z in np.transpose(self.Lap_mat))
	
   
   
    def buildRHS(self, bc_dict, charge_func='none', charge_grid='none'):
		rows = self.pmesh.nVert
		rhs  = np.zeros(rows)

		# Add Dirichlet BCs
		for cur_cont in self.pmesh.contacts.keys():
			val = bc_dict[cur_cont]

			# Loop over Dirichlet nodes
			for cur_node_index in self.pmesh.contacts[cur_cont]:
				rhs[cur_node_index] = val
				
		#
		# Add charge to all nodes
		# Two options: 
		#             1. Use function (x, y)-->charge_dens (charge_func)
		#             2. Map directly from a grid of semiconductor nodes
		#
		# Can be both, so charge is accumulated from function to grid.
		# If both options are 'none', Laplace solution is obtained
		#
		
		# Add charge using function
		if charge_func != 'none':
			charge_dens_func = charge_func
		elif charge_grid !='none':
			charge_dens_func = lambda x, y: charge_grid.nn_interp(x, y)
		
		# Loop over all elements
		# If a semiconductor element, add charge to all nodes of element
		for cur_elem in self.pmesh.getElementIter():
			if cur_elem.material != 'Silicon':
				continue
			
			# Shortcut references for local vertices	
			vIter = cur_elem.getVertexIter(self.pmesh)
			(gb0, v0) = vIter.next()
			(gb1, v1) = vIter.next()
			(gb2, v2) = vIter.next()
			(gb3, v3) = vIter.next()	
					
	  		# Get spacings for the element
	  		dx = v1.x-v0.x
	  		dy = v3.y-v0.y	
			
			# Loop over all nodes of this element
			# Only a quarter of the nodal area of each node is 
			# in the element
			for gb_index, vertex in cur_elem.getVertexIter(self.pmesh):
				#vertex = self.pmesh.vertex_list[global_index]
				x    = vertex.x
				y    = vertex.y 
				dA   = 0.25 * dx * dy # Node gets one quarter of element
				qval = -1.0*100*charge_dens_func(x, y)*dA*1e-14*phys.q/phys.eps0
				rhs[gb_index] += qval
							
				
		return rhs
   		
   		
   
    def buildLaplaceMat(self):
    	rows  = self.pmesh.nVert
      	Lap_mat_dok = sparse.dok_matrix((rows, rows))
      	
      	# Loop over all elements
        for cur_elem in self.pmesh.getElementIter():
      		# All vertex indices in element
      		#print 'Working on element: ', cur_elem
      		
      		
			# Shortcut references for local vertices	
			vIter = cur_elem.getVertexIter(self.pmesh)
			(gb0, v0) = vIter.next()
			(gb1, v1) = vIter.next()
			(gb2, v2) = vIter.next()
			(gb3, v3) = vIter.next()	
      		
      		# Get spacings for the element
			dx = v1.x-v0.x
			dy = v3.y-v0.y
      		
      		# Coupling coefficients
			eps = self.pmesh.prop_dict[cur_elem.material]['eps']
			tx = eps * 0.5 * dy / dx
			ty = eps * 0.5 * dx / dy
      		
			
			for loc_index, (global_index, cur_vert) in enumerate(
										cur_elem.getVertexIter(self.pmesh)):
      			# Sequence is:
      			# 3 -- 2
				# |    |
				# 0 -- 1
				# 
				# Construct partial control volume for each node
				
				#--------------------------------------------------
        		# Global row is equal to the global vertex index
				#--------------------------------------------------
        
        		# Is it a Dirichlet node?
        		# If so, simply assign a 1.0 into the matrix
				if cur_vert.boundary:
					Lap_mat_dok[global_index, global_index] = 1.0

				# Not a Dirchlet node	
				else:
					if loc_index == 0:
						# Couple right
						row_r = v1.index
						Lap_mat_dok[global_index, global_index] += tx
						Lap_mat_dok[global_index, row_r]        -= tx

						# Couple up
						row_u = v3.index
						Lap_mat_dok[global_index, global_index] += ty
						Lap_mat_dok[global_index, row_u]        -= ty

					if loc_index == 1:
						# Couple up
						row_u = v2.index
						Lap_mat_dok[global_index, global_index] += ty
						Lap_mat_dok[global_index, row_u]        -= ty

						# Couple left
						row_l = v0.index
						Lap_mat_dok[global_index, global_index] += tx
						Lap_mat_dok[global_index, row_l]        -= tx

					if loc_index == 2:
						# Couple down
						row_d = v1.index
						Lap_mat_dok[global_index, global_index] += ty
						Lap_mat_dok[global_index, row_d]        -= ty

						# Couple left
						row_l = v3.index
						Lap_mat_dok[global_index, global_index] += tx
						Lap_mat_dok[global_index, row_l]        -= tx

					if loc_index == 3:
						# Couple down
						row_d = v0.index
						Lap_mat_dok[global_index, global_index] += ty
						Lap_mat_dok[global_index, row_d]        -= ty

						# Couple right
						row_r = v2.index
						Lap_mat_dok[global_index, global_index] += tx
						Lap_mat_dok[global_index, row_r]        -= tx
        				
        				
        					
			
        return Lap_mat_dok.tocsr()
                                                                                                                                                                                ingaas.py                                                                                           0000755 0002735 0004704 00000014666 13167755426 013502  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import ParabolicValley as pvalley
import SP2D as sp2d
import SCF  as SCF
import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys
from   plot_fun import set_plt_defaults
import SuperLattice as super_l


import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib import colors, ticker, cm

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=300)


#+------------------------------------------------------------+
#|
#| Create Structure
#|
#+------------------------------------------------------------+

# Fin dimensions
L = 5.0
W_target = 45.0
IL = 0.5


# Superlattice
thick     = 0.5
dy_target = 0.25*thick
period    = 5.0*thick
n_periods = np.floor(W_target/period)
W         = n_periods * period
offset    = 0
energy    = 0.0
sl = super_l.SuperLattice(period, thick, energy, offset)




print 'Height=', W, ' period=', period

#dy_target = 5.0/30

ny = int(W/dy_target)
print 'ny: ', ny
#exit()

# grid
nx = 30
#ny = 50
nh = 20
nv = 15
no = 5

#
# Finite-Difference Mesh: Used by Schroedinger
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.2
EOT = 0.85
#vsp = 2.4923
vsp = 2.5
hsp = IL + (EOT-IL)*Eps_HiK/3.9
print 'hsp: ', hsp

'''
pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)
'''
							
#pgaa_zoner = zoner.MLFETZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)
							
pgaa_zoner = zoner.finFETZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)
							
#pgaa_zoner = zoner.NWZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)

#
# Finite-Volume Mesh: used by Poisson
mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			

# Assign the contact, call it 'Gate')			
mesh.assignContact( pgaa_zoner.contact_func, 'Gate' )

#+---------------------------------------------------------+
#|
#| Poisson-related Objects
#|
#+---------------------------------------------------------+

# Create Poisson Object
poisson = pois.FVMPoisson(mesh)

# Extract Laplace matrix; fixed dueing simulation
pmat =  poisson.Lap_mat

# Set Boundary Values on contacts
bc_val_dict = { 'Gate' : 0.7 }


#+---------------------------------------------------------+
#|
#| Schroedinger-related Objects
#|
#+---------------------------------------------------------+


# Create Delta4 Valleys
# 'Left' Valleys
mstar = np.array([[0.07, 0, 0],[0, 0.07, 0],[0, 0, 0.07]])
mult  = 1
Gamma = pvalley.ParabolicValley('Gamma', L, W, mstar, mult)



# Create Schroedinger Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Gamma], fdmesh)

#+---------------------------------------------------------+
#|
#| SCF Solver: Self-Consistent Solution of QM and Poisson
#|
#+---------------------------------------------------------+

scf = SCF.SCF( sp, poisson, sl )
Ef = 0.0
(iterations, converged) = scf.solve(
								bc_val_dict, 
								max_iter = 200, 
								tol      = 2e-3, 
								alpha    = 0.1,
								n_eigen  = 20,
								Ef       = Ef 
							)

print 'Converged=', converged, ' after ', iterations, ' iterations.'


#set_plt_defaults()
plt.rcParams['figure.facecolor'] = 'w'
plt.rc('xtick', labelsize=18) 
plt.rc('ytick', labelsize=18) 
font = {'family' : 'Bitstream Vera Sans',
	    'weight' : 'bold',
	    'size'   : 16}

plt.rc('font', **font)
plt.rcParams['toolbar'] = 'None'
#plt.tick_params(which='major', length=15, width=2)
#plt.tick_params(which='minor', length=5, width=2)

plt.rcParams['axes.linewidth'] = 2 

# Done with loop, plot results
fig1 = plt.figure(1)

# Potential sub-figure
(Xp, Yp, Z_final_pot) = poisson.pmesh.mapSolutionToGrid(scf.final_potential)

ax = fig1.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z_final_pot), np.amax(Z_final_pot), 30)
cp = plt.contourf(Xp, Yp, Z_final_pot, levels=levels)
plt.contour(Xp, Yp, Z_final_pot, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	

for cur_rect in pgaa_zoner.plotting_rectangle_tuple():
	start_point = cur_rect[0]
	length      = cur_rect[1]
	height      = cur_rect[2]
	ax.add_patch(
		patches.Rectangle(
		    start_point,
		    length,
		    height,
		    fill=False,      # remove background
		    linewidth=2
		)
	)



# QM Charge sub-figure
fig2 = plt.figure(2)
x2 = fig2.add_subplot(111, aspect='equal')
wfn_dict          = sp.getSolution(0)
(Xfd,Yfd,Z_eigen) = wfn_dict['wfns']
ZCharge           = sp.computeChargeDensity(Ef)
levels_c = np.linspace(np.amin(ZCharge), np.amax(ZCharge), 30)
cp_QM = plt.contourf(Xfd, Yfd, ZCharge, levels=levels_c)
plt.colorbar(cp_QM)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')

# Get DOS for each valley
(energy, ret_val) = sp.computeDOS(e_min=-0.1, e_max=0.2)
Gamma_DOS = ret_val[0]  # Symmetric D4 valleys, get only one, multiply by 2
plt.figure(3)
#plt.plot(energy, 2*D4_DOS, '-r', linewidth=3)
plt.plot(energy, Gamma_DOS, '-b', linewidth=3)
plt.xlabel('Energy [eV]')
plt.ylabel('DOS')
plt.grid(True)
#plt.xlim(xmin=-0.1, xmax=0.2)
plt.ylim(ymin=0.0, ymax=100)

# Plot eigenfunctions
# Separate plot for each valley


for cur_valley in range(0,1):
	val_dict = sp.getSolution(cur_valley)
	(X,Y,Z) = val_dict['wfns']
	(rows, cols, slabs) = Z.shape
	nrows = 4
	ncols = 3
	div = (slabs-(slabs % ncols))/ncols + 1
	if div > nrows:
		div = nrows
	print 'div: ', div
	fig, axarr = plt.subplots(div, ncols)
	
	for row in range(0,div):
		for col in range(0,ncols):
			cur_slab = row * ncols + col
			if cur_slab == slabs:
			    break
			plotZ = Z[:,:,cur_slab] * Z[:,:,cur_slab]
			
			levels = np.linspace(np.amin(plotZ), np.amax(plotZ), 30)
			axarr[row, col].contourf(X, Y, plotZ, levels=levels)
			axarr[row, col].set_aspect('equal')
			axarr[row, col].xaxis.set_visible(False)
			axarr[row, col].yaxis.set_visible(False)
			
			
# Compute charge and velocity

(vel, tot_charge, valley_charge) = sp.computeChargeAndVelocity(Ef)
sidewall_charge = 1e14 * tot_charge / (2.0*(W+vsp))
nw_charge       = 1e14 * tot_charge / (2.0*(W+vsp+L+hsp))
#d2_occupancy = valley_charge[2]/(tot_charge)

print 'vel=', vel*1e-7, ' tot_charge=', tot_charge
#print 'valley charge: ', valley_charge



#sw_total_charge = 1e-14*sw_charge[0]*(2.0*W)
#print 'Standard SW charge: ', sw_charge, ' SW total charge: ', sw_total_charge
plt.show()


                                                                          meshfun.py                                                                                          0000644 0002735 0004704 00000001572 13020130274 013642  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               import numpy as np

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
                                                                                                                                      nw.py                                                                                               0000755 0002735 0004704 00000006170 13023560035 012630  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import valley as valley
import FDMesh as fdmesh
import Schroedinger as sch

import numpy  as np
import scipy.interpolate as interp
import scipy.integrate   as integral
import scipy.optimize    as optimize
import phys

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
#mstar = np.array([[0.19, -0.478, 0],[-0.478, 0.315, 0],[0, 0, 0.315]])
mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
mult  = 1

# Create Delta2 Valley
Delta4 = valley.valley('Delta4', L, W, mstar, mult)
# grid
nx = 50
ny = 50
fdmesh = fdmesh.FDMesh(L, W, nx, ny)
SchrD2  = sch.Schroedinger(Delta4, fdmesh)

# Create a CB profile
peak = 0.250
cbfun = lambda x, y : x * (L-x) * peak * 4.0/(L*L)
curvature = 1e18 * peak * 2.0 * 4.0/(L*L)
charge_density = curvature * phys.eps0 * phys.eps_Si * 1e-6 / phys.q
print 'charge density: ', charge_density
total_charge = charge_density * L * W * 1e-7 * 1e-7
surface_charge = 0.5 * total_charge / (W * 1e-7)


# compute surface charge density
print 'surface charge dens: ', surface_charge


(evals, (X,Y,Z)) = SchrD2.solve(n_eigen=12, CB=cbfun)
print 'evals: ', evals

# Make some plots

# DOS
n_energy_points = 300
energy = np.linspace(evals[0]-0.05, evals[0]+0.1, n_energy_points)
DOS = Delta4.computeDOS(energy, evals)

'''
figDOS = plt.figure()
plt.plot(energy, DOS, '-b', linewidth=2)
plt.xlabel('Energy [eV]')
plt.ylabel('DOS [1/eV cm]')
plt.grid(True)
'''

# Charge density
Ef = 0.3
occDOS = DOS * Delta4.occupancy( energy, Ef )
DOS_interp = interp.interp1d( energy, DOS )


fermi_points = 50
fermi_vals  = np.linspace(energy[0]-3.0*phys.kT, energy[-1], fermi_points)
charge_vals = np.zeros(fermi_points)
for i in range(0, fermi_points):
    if i % 5 == 0:
        print 'Working on Fermi level:', fermi_vals[i]
    charge_vals[i] = Delta4.computeSidewallCharge(fermi_vals[i], evals)

'''
# Compute self-consistent Fermi level
RHS = lambda x : Delta4.computeQMCharge( x, evals ) - surface_charge
x0 = 0.15
xf= optimize.fsolve(RHS, x0, xtol=1e-3)
print 'Found Fermi level : ', xf
'''

# Now make some evil plots
plt.figure()
plt.plot(energy, DOS_interp(energy), '-b', linewidth=3)
plt.plot(energy, occDOS, '-r', linewidth=3)
plt.xlabel('Energy [eV]')
plt.ylabel('DOS')
plt.grid(True)


plt.figure()
plt.plot(fermi_vals, charge_vals, '-k', linewidth=3)
plt.xlabel('Fermi Level [eV]')
plt.ylabel('Charge Density [1/cm^2]')
plt.grid(True)

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
                                                                                                                                                                                                                                                                                                                                                                                                        ParabolicValley.py                                                                                  0000644 0002735 0004704 00000004042 13047110236 015246  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               import numpy as np
import scipy.interpolate as interp
import scipy.integrate   as integral
import phys as phys
from   valley import *

class ParabolicValley(valley):
   
    def __init__(self, name, L, W, mstar, mult, shift=0.0):
        """
        Constructor calls base class for base things, but
        adds custom construction explicity
        """
        super(ParabolicValley, self).__init__(name, L, W, 
                                              mstar, mult, shift)
        
        
        
    def singleStateDOS(self, E, Ec):
        """
        Compute single-state DOS at energy E
        Usese bound-state energy Ec
        """
        prefac = 2.0 / (phys.pi * phys.hbar)
        root   = np.sqrt(phys.m0 * self.mstarDOS/(2.0 * phys.q * (E-Ec)))
        return 1e-9 * prefac * root * self.mult * phys.q
        
    
    
    def singleStateVel(self, Ef, Ec, limit=300, epsrel=1e-5):
        """
        Computes the velocity integral for the single bound state
        Returns both the velocity integral and the associated charge
        Overall vel is assembled from calls to singleStateVel:
        
        vel = sum(vel_i)/sum(q_i)
        
        where vel_i is the velocity integral from the single state, and
        q_i is the charge integral from the single state.
        
        NOTE: vel_i actually has units of velocity * charge
        """
        v_integrand = lambda x : \
                      1e-7 * 1e-4 * 1e21 * self.singleStateDOS(x, Ec) * \
                      self.occupancy(x, Ef) * \
                      np.sqrt((2.0*phys.kT*phys.q /
                              (phys.pi*self.mstar[2, 2]*phys.m0)) * (x-Ec))
                               

        (velnum, abserr) = integral.quad( v_integrand,
                                            Ec, np.inf,
                                            limit=limit,
                                            epsrel=epsrel )
                                            
                                            
                                            
        return velnum
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              parm_subs.py                                                                                        0000644 0002735 0004704 00000001767 13032517643 014212  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #+-------------------------------------------------------------
#|
#| Function for parameter substitution
#| Takes a template file name (infile), output file 
#| name (outfile), and a dictionary of name, value pairs
#| (val_dict).
#| Expects the template to use the following substitution 
#| format:
#| text @keyword@ text ...
#| with the result:
#| text keyword_value text ...
#|
#+-------------------------------------------------------------- 



def parm_subs(infile, outfile, val_dict):

	inf  = open(infile, 'r')
	outf = open(outfile, 'w')
	
	for cur_line in inf:
		#print 'Got: ', cur_line
		
		# Does the line contain substitution targets?
		if '@' not in cur_line:
			outf.write(cur_line)
		else:
			# There may be substitution targets
			for cur_word in val_dict.keys():
				if '@'+cur_word+'@' in cur_line:
					#print 'Substitution for ', cur_word
					cur_line = cur_line.replace('@'+cur_word+'@', 
							    str(val_dict[cur_word]))
			outf.write(cur_line)
											
		
	outf.close()
	inf.close()
         phys.py                                                                                             0000644 0002735 0004704 00000000307 13020130274 013153  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               # Phsyical constants


h       = 6.626e-34
m0      = 9.10938356e-31
pi      = 3.14159265
hbar    = h/(2.0*pi)
q       = 1.6e-19
kT      = 26e-3
nm_to_m = 1e-9
eps0    = 8.85418782e-12
eps_Si  = 11.9
                                                                                                                                                                                                                                                                                                                         plotall.py                                                                                          0000755 0002735 0004704 00000004167 13047432547 013673  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python 

import pandas as pnd
import numpy  as np
import glob
import os
import matplotlib.pyplot as plt
from   FOMtoPandas import FOMtoPandas
import plot_fun as pf

# Create a new table, add a row for each split
# Use data from FOM.dat in each split

plot_types = ['pgaa', 'nw', 'fin']

df = FOMtoPandas()
print df.head()
df = df[df['converged']=='True']
df = df[df['vsp'] !=3.0 ]

# Average charge in channel slice
df['average_charge'] = 1e21*df['tot_charge']/(df['L']*df['W'])

# Charge per unit height
df['charge_height'] = 1e14*df['tot_charge']/(df['W']+ 2.0*df['vsp'])

# Create sidewall charge
df['hsp']    = 2.49
#df['fin_ns'] = df['tot_charge']*1e14/(2.0*df['W'] + df['L'] +
#										2.0*df['hsp'] + 2.0*df['vsp'])
										
df['fin_ns'] = df['tot_charge']*1e14/(2.0*df['W'] + df['L'])
										
df['sw_ns']  = df['tot_charge']*1e14/(2.0*(df['W']+2.0*df['vsp']))

df['gaa_ns'] = df['tot_charge']*1e14/(2.0*(df['W']+df['L'] + 
										   2.0*df['vsp'] + 2.0*df['hsp']))
										 
df['gate_ns'] = 0.0
df.ix[(df['dev_type']=='fin'), 'gate_ns'] = df.ix[ (df['dev_type']=='fin'), 'fin_ns']
df.ix[(df['dev_type']=='pgaa'),'gate_ns'] = df.ix[ (df['dev_type']=='pgaa'),'sw_ns'] 
df.ix[(df['dev_type']=='nw'), 'gate_ns'] = df.ix[ (df['dev_type']=='nw'), 'gaa_ns']  
		
df['vel'] = df['vel']*1.05								   
						 

df.to_csv(path_or_buf='done.csv')
print df.head()


# Make some evil plots
lines   = ['-','--','-.',':','-o', '-D']
colors  = ['r','b','k','m']
symbols = ['-o','D','s','^','D']
markersize = 7
lw  = 3

# Plotting
pf.set_plt_defaults()

cat_var = 'vsp'

# Create a list of data frames, one for each plot type
dev_df = [ df[df['dev_type']==x] for x in plot_types ]

for dev_index, cur_dev_df in enumerate(dev_df):
	cur_col = colors[dev_index]
	vsp_val = cur_dev_df['vsp'].unique().max()
	cat_df  = cur_dev_df[cur_dev_df['vsp']==vsp_val]
	
	
	for index, cur_height in enumerate(np.sort(cat_df['W'].unique())):
		sdf = (cat_df[cat_df['W']==cur_height]).sort(columns=['Vg'], ascending=True)
		cur_line = lines[index] + cur_col	
		
		pf.plot_fun(sdf, cur_line, lw)
	
plt.show()
	
	
	
	


		
	
                                                                                                                                                                                                                                                                                                                                                                                                         plot_fun.py                                                                                         0000644 0002735 0004704 00000013167 13033473624 014043  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               import matplotlib.pyplot as plt
from   matplotlib.ticker import FormatStrFormatter
import numpy             as np
import scipy.interpolate as sp_interp
import scipy.optimize    as optimize


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



def plot_fun(sdf, cur_line, lw):		
		
		
		
		fignum = 1	
	
		# Charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['gate_ns']
		plt.semilogy(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg [V]')
		plt.ylabel('Gate ns [1/cm^2]')
		plt.grid(True)
	
	

		# Gate Capacitance
		plt.figure(fignum)
		fignum += 1
		vvals = sdf['Vg'].values
		qvals = sdf['gate_ns']
			# Interpolate onto sparse mesh for smoothing
		print 'interpolating: vvals:', vvals.shape, ' qvals:', qvals.shape
		#print sdf 
		q_interpolator = sp_interp.interp1d(vvals, qvals, kind='cubic')
		v_i = np.linspace(vvals[0], vvals[-1], 10)
		q_i = q_interpolator(v_i)
		dV  = v_i[1]-v_i[0]
		cvals = 1e6*1.6e-19*np.gradient(q_i, dV)
			# Now interpolate back to original grid
		c_interpolator = sp_interp.interp1d(v_i, cvals, kind='cubic')
		cvals_final = c_interpolator(vvals)
		plt.plot(vvals, cvals_final, cur_line, linewidth=lw)
		plt.xlabel('Vg [V]')
		plt.ylabel('Gate Capacitance [uF/cm^2]')
		plt.grid(True)
		plt.ylim(ymin=0.0, ymax=3.5)
		
		# Gate Capacitance for total charge
		plt.figure(fignum)
		fignum += 1
		vvals = sdf['Vg'].values
		qvals_tot = sdf['tot_charge']*1e14
			# Interpolate onto sparse mesh for smoothing
		
		q_interpolator_tot = sp_interp.interp1d(vvals, qvals_tot, kind='cubic')
		v_i_tot = np.linspace(vvals[0], vvals[-1], 10)
		q_i_tot = q_interpolator_tot(v_i_tot)
		dV_tot  = v_i[1]-v_i[0]
		cvals_tot = 1e6*1.6e-19*np.gradient(q_i_tot, dV_tot)
			# Now interpolate back to original grid
		c_interpolator_tot = sp_interp.interp1d(v_i_tot, cvals_tot, kind='cubic')
		cvals_final_tot = c_interpolator_tot(vvals)
		if sdf['dev_type'].values[0] != 'fin':
			plt.plot(vvals, cvals_final_tot, cur_line, linewidth=lw)
		plt.xlabel('Vg [V]')
		plt.ylabel('Gate Capacitance [uF/cm]')
		plt.grid(True)
		plt.ylim(ymin=0.0, ymax=100)
	
	
	
		# Gate capacitance with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg'].values
		dV    = xvals[1]-xvals[0]
		qvals = sdf['gate_ns']
		cvals = 1e6*1.6e-19*np.gradient(qvals, dV)
			# Interpolant for capacitance
		#c_interp = sp_interp.interp1d(xvals, cvals)
		rhs = lambda x : c_interpolator(x)-0.15
			# Now find shift require to put cval=1.5 to Vg=0.0
		v0 = -0.05
		vf= optimize.fsolve(rhs, v0, xtol=1e-6)
		print 'Got vf=', vf
		yvals = cvals_final
		plt.plot(vvals-vf, cvals_final, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt [V]')
		plt.ylabel('Gate Capacitance [uF/cm^2]')
		plt.grid(True)
		plt.ylim(ymin=0.0, ymax=3.5)
		
	
		# Charge with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['gate_ns']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt [V]')
		plt.ylabel('Gate ns [1/cm^2]')
		plt.grid(True)
		
		# Charge with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['charge_height']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt [V]')
		plt.ylabel('Charge / height')
		plt.grid(True)
		
		# Total Charge with Vt shift
		# Don't plot for FinFET
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = 1e7*sdf['tot_charge']
		if sdf['dev_type'].values[0] != 'fin':
			plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt [V]')
		plt.ylabel('Total Charge [1/nm]')
		plt.grid(True)
		
		# Average Charge with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['average_charge']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt [V]')
		plt.ylabel('Average Charge Conc. [1/cm^3]')
		plt.grid(True)
	
	
		# Delta2 Fraction vs Vg-Vt
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['d2']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt, [V]')
		plt.ylabel('Delta 2 Fraction')
		plt.grid(True)
	
		# Delta2 Fraction vs sidewall charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['gate_ns']
		yvals = sdf['d2']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Gate ns [1/cm^2]')
		plt.ylabel('Delta 2 Fraction')
		plt.grid(True)
		plt.xlim(xmin=2e11, xmax=2e13)


		# Velocity vs. Vg-Vt
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['vel']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt [V]')
		plt.ylabel('Injection Velocity [cm/s]')
		plt.grid(True)

		# Velocity vs. sidewall charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['gate_ns']
		yvals = sdf['vel']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Gate ns [1/cm^2]')
		plt.ylabel('Injection Velocity [cm/s]')
		plt.grid(True)
		plt.xlim(xmin=2e11, xmax=2e13)
		
		# Velocity vs. average charge conc.
		
		
		
		plt.figure(fignum)
		plt.tick_params(axis='x', which='minor')
		ax = plt.gca()
		#ax.set_yscale('log')
		plt.tick_params(axis='x', which='minor')
		ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
		fignum += 1
		xvals = sdf['average_charge']
		yvals = sdf['vel']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Average Charge Conc. [1/cm^3]')
		plt.ylabel('Injection Velocity [cm/s]')
		plt.grid(True)
		plt.xlim(xmin=2e18, xmax=2e20)
                                                                                                                                                                                                                                                                                                                                                                                                         plot_nw.py                                                                                          0000755 0002735 0004704 00000010540 13032252544 013665  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python 

import pandas as pnd
import numpy  as np
import glob
import os
import matplotlib.pyplot as plt
import scipy.interpolate as sp_interp
import scipy.optimize    as optimize
from   FOMtoPandas import FOMtoPandas

# Create a new table, add a row for each split
# Use data from FOM.dat in each split

df = FOMtoPandas()
print df.head()
df = df[df['converged']=='True']
df = df[ (df['vsp']==2.4923) | (df['vsp']==2.5)]
df.to_csv(path_or_buf='done.csv')
print df.head()

# Make some plots
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

lines   = ['-','--','-.',':','-o', '-D']
colors  = ['r','b','k','m']
symbols = ['-o','D','s','^','D']
markersize = 7
lw  = 3

# Plotting
set_plt_defaults()

cat_var = 'vsp'

non_fin_df = df[df['dev_type'] != 'fin']
fin_df     = df[df['dev_type'] == 'fin']

for cat_index, cat_val in enumerate(np.sort(non_fin_df[cat_var].unique())):
	cat_df = df[ (df[cat_var]==cat_val ) & (df['dev_type'] !='fin') ]
	cur_col = colors[cat_index]

	print 'dev_type: ', df['dev_type']

	for index, cur_height in enumerate(np.sort(non_fin_df['W'].unique())):
		sdf = (cat_df[cat_df['W']==cur_height]).sort(columns=['Vg'], ascending=True)
		cur_line = lines[index] + cur_col	
		#cur_line = symbols[index] + 'b'
		fignum = 1	
		
		# Charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['sw_charge']
		plt.semilogy(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg')
		plt.ylabel('ns')
		plt.grid(True)
		
		
	
		# Gate Capacitance
		plt.figure(fignum)
		fignum += 1
		vvals = sdf['Vg'].values
		qvals = sdf['sw_charge']
			# Interpolate onto sparse mesh for smoothing
		print 'interpolating: vvals:', vvals.shape, ' qvals:', qvals.shape
		#print sdf 
		q_interpolator = sp_interp.interp1d(vvals, qvals, kind='cubic')
		v_i = np.linspace(vvals[0], vvals[-1], 10)
		q_i = q_interpolator(v_i)
		dV  = v_i[1]-v_i[0]
		cvals = 1e6*1.6e-19*np.gradient(q_i, dV)
			# Now interpolate back to original grid
		c_interpolator = sp_interp.interp1d(v_i, cvals, kind='cubic')
		cvals_final = c_interpolator(vvals)
		plt.plot(vvals, cvals_final, cur_line, linewidth=lw)
		plt.xlabel('Vg')
		plt.ylabel('Gate Capacitance')
		plt.grid(True)
		
		
		
		# Gate capacitance with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg'].values
		dV    = xvals[1]-xvals[0]
		qvals = sdf['sw_charge']
		cvals = 1e6*1.6e-19*np.gradient(qvals, dV)
			# Interpolant for capacitance
		c_interp = sp_interp.interp1d(xvals, cvals)
		rhs = lambda x : c_interp(x)-1.0
			# Now find shift require to put cval=1.5 to Vg=0.0
		v0 = 0.0
		vf= optimize.fsolve(rhs, v0, xtol=1e-5)
		print 'Got vf=', vf
		yvals = cvals_final
		plt.plot(vvals-vf, cvals_final, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('Gate Capacitance')
		plt.grid(True)
		
		# Charge with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['sw_charge']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('ns')
		plt.grid(True)
		
		
		# Delta2 Fraction vs Vg-Vt
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['d2']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('d2_frac')
		plt.grid(True)
		
		# Delta2 Fraction vs sidewall charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['sw_charge']
		yvals = sdf['d2']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg')
		plt.ylabel('d2_frac')
		plt.grid(True)
		plt.xlim(xmin=1e11, xmax=1e14)
	
	
		# Velocity vs. Vg-Vt
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['vel']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('velocity')
		plt.grid(True)
	
		# Velocity vs. sidewall charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['sw_charge']
		yvals = sdf['vel']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('ns')
		plt.ylabel('velocity')
		plt.grid(True)
		plt.xlim(xmin=1e11, xmax=1e14)
	
plt.show()
	
	
	
	


		
	
                                                                                                                                                                plot_pgaa.py                                                                                        0000755 0002735 0004704 00000010121 13032001120 014122  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python 

import pandas as pnd
import numpy  as np
import glob
import os
import matplotlib.pyplot as plt
import scipy.interpolate as sp_interp
import scipy.optimize    as optimize
from   FOMtoPandas import FOMtoPandas

# Create a new table, add a row for each split
# Use data from FOM.dat in each split

df = FOMtoPandas()
print df.head()
df = df[df['converged']=='True']
df = df[df['vsp'] !=2.4923]
df.to_csv(path_or_buf='done.csv')
print df.head()

# Make some plots
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

lines   = ['-','--','-.',':','-o', '-D']
colors  = ['r','b','k','m']
symbols = ['-o','D','s','^','D']
markersize = 7
lw  = 3

# Plotting
set_plt_defaults()

cat_var = 'vsp'


for cat_index, cat_val in enumerate(np.sort(df[cat_var].unique())):
	cat_df = df[df[cat_var]==cat_val]
	cur_col = colors[cat_index]

	for index, cur_height in enumerate(np.sort(df['W'].unique())):
		sdf = (cat_df[cat_df['W']==cur_height]).sort(columns=['Vg'], ascending=True)
		cur_line = lines[index] + cur_col	
		#cur_line = symbols[index] + 'b'
		fignum = 1	
		
		# Charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['sw_charge']
		plt.semilogy(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg')
		plt.ylabel('ns')
		plt.grid(True)
		
		
	
		# Gate Capacitance
		plt.figure(fignum)
		fignum += 1
		vvals = sdf['Vg'].values
		qvals = sdf['sw_charge']
			# Interpolate onto sparse mesh for smoothing
		q_interpolator = sp_interp.interp1d(vvals, qvals, kind='cubic')
		v_i = np.linspace(vvals[0], vvals[-1], 10)
		q_i = q_interpolator(v_i)
		dV  = v_i[1]-v_i[0]
		cvals = 1e6*1.6e-19*np.gradient(q_i, dV)
			# Now interpolate back to original grid
		c_interpolator = sp_interp.interp1d(v_i, cvals, kind='cubic')
		cvals_final = c_interpolator(vvals)
		plt.plot(vvals, cvals_final, cur_line, linewidth=lw)
		plt.xlabel('Vg')
		plt.ylabel('Gate Capacitance')
		plt.grid(True)
		
		
		
		# Gate capacitance with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg'].values
		dV    = xvals[1]-xvals[0]
		qvals = sdf['sw_charge']
		cvals = 1e6*1.6e-19*np.gradient(qvals, dV)
			# Interpolant for capacitance
		c_interp = sp_interp.interp1d(xvals, cvals)
		rhs = lambda x : c_interp(x)-1.0
			# Now find shift require to put cval=1.5 to Vg=0.0
		v0 = 0.0
		vf= optimize.fsolve(rhs, v0, xtol=1e-5)
		print 'Got vf=', vf
		yvals = cvals_final
		plt.plot(vvals-vf, cvals_final, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('Gate Capacitance')
		plt.grid(True)
		
		# Charge with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['sw_charge']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('ns')
		plt.grid(True)
		
		
		# Delta2 Fraction vs Vg-Vt
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['d2']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('d2_frac')
		plt.grid(True)
		
		# Delta2 Fraction vs sidewall charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['sw_charge']
		yvals = sdf['d2']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg')
		plt.ylabel('d2_frac')
		plt.grid(True)
		plt.xlim(xmin=1e11, xmax=1e14)
	
	
		# Velocity vs. Vg-Vt
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['vel']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('velocity')
		plt.grid(True)
	
		# Velocity vs. sidewall charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['sw_charge']
		yvals = sdf['vel']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('ns')
		plt.ylabel('velocity')
		plt.grid(True)
		plt.xlim(xmin=1e11, xmax=1e14)
	
plt.show()
	
	
	
	


		
	
                                                                                                                                                                                                                                                                                                                                                                                                                                               plot-stress.py                                                                                      0000755 0002735 0004704 00000004472 13051374400 014506  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python 

import pandas as pnd
import numpy  as np
import glob
import os
import matplotlib.pyplot as plt
from   FOMtoPandas import FOMtoPandas
import plot_fun as pf

# Create a new table, add a row for each split
# Use data from FOM.dat in each split

plot_types = ['pgaa', 'fin']

df = FOMtoPandas()
print df
#exit()
df = df[df['converged']=='True']
#df = df[df['vsp'] !=3.0 ]

# Average charge in channel slice
df['average_charge'] = 1e21*df['tot_charge']/(df['L']*df['W'])

# Charge per unit height
df['charge_height'] = 1e14*df['tot_charge']/(df['W']+ 2.0*df['vsp'])

# Create sidewall charge
df['hsp']    = 2.49
#df['fin_ns'] = df['tot_charge']*1e14/(2.0*df['W'] + df['L'] +
#										2.0*df['hsp'] + 2.0*df['vsp'])
										
df['fin_ns'] = df['tot_charge']*1e14/(2.0*df['W'] + df['L'])
										
df['sw_ns']  = df['tot_charge']*1e14/(2.0*(df['W']+2.0*df['vsp']))

df['gaa_ns'] = df['tot_charge']*1e14/(2.0*(df['W']+df['L'] + 
										   2.0*df['vsp'] + 2.0*df['hsp']))
										 
df['gate_ns'] = 0.0
df.ix[(df['dev_type']=='fin'), 'gate_ns'] = df.ix[ (df['dev_type']=='fin'), 'fin_ns']
df.ix[(df['dev_type']=='pgaa'),'gate_ns'] = df.ix[ (df['dev_type']=='pgaa'),'sw_ns'] 
df.ix[(df['dev_type']=='nw'), 'gate_ns'] = df.ix[ (df['dev_type']=='nw'), 'gaa_ns']  
		
df['vel'] = df['vel']*1.05								   
						 

df.to_csv(path_or_buf='done.csv')
print df.head()


# Make some evil plots
lines   = ['-','--','-.',':','-o', '-D']
colors  = ['r','k','b','m']
symbols = ['o', '']
markersize = 7
lw  = 3

# Plotting
pf.set_plt_defaults()

cat_var = 'd2_shift'

print df['d2_shift'].unique()
#exit()

# Create a list of data frames, one for each plot type
dev_df = [ df[df['dev_type']==x] for x in plot_types ]

for dev_index, cur_dev_df in enumerate(dev_df):
	cur_col = colors[dev_index]
	
	for d2_index, d2_val in enumerate(sorted(cur_dev_df['d2_shift'].unique())):
		print ' >>> Working on d2 = ', d2_val
		cat_df  = cur_dev_df[cur_dev_df['d2_shift']==d2_val]
		cur_symbol = symbols[d2_index]
		if cur_col == 'k':
			cur_symbol = ''
	
		for index, cur_height in enumerate(np.sort(cat_df['W'].unique())):
			sdf = (cat_df[cat_df['W']==cur_height]).sort(columns=['Vg'], ascending=True)
			cur_line = lines[index] + cur_symbol + cur_col	
		
			pf.plot_fun(sdf, cur_line, lw)
	
plt.show()
	
	
	
	


		
	
                                                                                                                                                                                                      PoissonMesh.py                                                                                      0000644 0002735 0004704 00000004011 13024575263 014452  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               # Mesh for Poisson
# Superset of the Schroedinger mesh

import FDMesh as fdmesh
import numpy  as np

class PoissonMesh:
	
	def __init__(self, schr_mesh, vsp, hsp, eps_ox=3.9, 
							eps_semi=11.9, nh=3, nv=3):
		self.schr_mesh    = schr_mesh
		self.eps_ox       = eps_ox
		self.eps_semi     = eps_semi
		self.vsp          = vsp
		self.hsp          = hsp
		self.nh           = nh
		self.nv           = nv
		if self.nh <2:
			self.nh = 2
		if self.nv <2:
			self.nv = 2
			
		self.nx = self.schr_mesh.nx + 2*(self.nh-1)
		self.ny = self.schr_mesh.ny + 2*(self.nv-1) 
		
		self.buildMesh()
		
	def getXMesh(self):
		return self.xmesh	
		
	def getYMesh(self):
		return self.ymesh
		
	def getXInterfaces(self):
		return (self.nh-1, self.nh-1+self.schr_mesh.nx-1)
		
	def getYInterfaces(self):
		return (self.nv-1, self.nv-1+self.schr_mesh.ny-1)
		
	def getMeshSpacing(self):
		return (self.dx, self.dy)
		
	
		
	# Build the arrays of x and y coordinates
	# This mesh is non-uniform
	# The "center" part of the mesh comes from th Schroedinger Mesh
	# The outer portions are insulator regions, used only for Poisson	
	def buildMesh(self):
	
		#------ Horizontal grid
		self.dx = self.hsp / (self.nh-1)
		
		# Additional vertical grids on left side
		xl = np.linspace(-self.hsp, 0-self.dx, self.nh-1)
		
		# Additional vertical grids on right side
		xr = np.linspace(self.schr_mesh.L+self.dx, 
						 self.schr_mesh.L+self.hsp, self.nh-1)
		
		# FDMesh horizontal grid
		xc = np.linspace(0.0, self.schr_mesh.L, self.schr_mesh.nx)
		
		self.xmesh = np.concatenate((xl, xc, xr))
		print 'xmesh: ', self.xmesh 
		
		#------ Vertical grid
		self.dy = self.vsp / (self.nv-1)
		
		# Additional vertical grids on bottom
		yb = np.linspace(-self.vsp, 0-self.dy, self.nv-1)
		
		# Additional vertical grids on top
		yt = np.linspace(self.schr_mesh.W+self.dy, 
						 self.schr_mesh.W+self.vsp, self.nv-1)
		
		# FDMesh vertical grid
		yc = np.linspace(0.0, self.schr_mesh.W, self.schr_mesh.ny)
		
		self.ymesh = np.concatenate((yb, yc, yt))
		print 'ymesh: ', self.ymesh
	
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       Poisson.py                                                                                          0000644 0002735 0004704 00000030265 13026004325 013633  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               
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
                                                                                                                                                                                                                                                                                                                                           ptest.py                                                                                            0000755 0002735 0004704 00000003707 13026311075 013347  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python


import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys

import scipy.sparse.linalg as la_sparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=2)
np.set_printoptions(linewidth=300)


L = 3.0
W = 3.0

# grid
nx = 20
ny = 20
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

vsp = 3.0
hsp = 3.0


# Build a mesh using uniform x and y-spacing
# (uniformtity not required)
# Simple zoner function makes left half Silicon, right half Oxide
# Properties dictionary defines permittivities for the two materials
mesh = FVM.FVMMesh( np.linspace(0.0, L, nx), 
					np.linspace(0.0, W, ny),
					lambda x, y: 'Silicon' if x<0.5*L else 'Oxide', 
					{ 'Silicon' : { 'eps' : 11.9 },
					  'Oxide'   : { 'eps' : 3.9  } } ) 
					  
# Assign contacts (Dirichlet BCs) to left and right edges			
mesh.assignContact( lambda x, y: x<1e-8,   'Left' )
mesh.assignContact( lambda x, y: L-x<1e-8, 'Right')
					
print(mesh)

poisson = pois.FVMPoisson(mesh)

#print 'Lap mat:'
pmat =  poisson.Lap_mat
#print pmat.toarray()


# Now build rhs mat
bc_val_dict = { 'Left' : 0.5, 'Right' : 1.0 }
rhs = poisson.buildRHS(bc_val_dict)
#print 'RHS: ', rhs
print 'RHS shape: ', rhs.shape
print 'Poisson shape: ', pmat.shape


# Solve the resulting system
z = la_sparse.spsolve(pmat, rhs)
#print 'z:', z

(X, Y, Z) = poisson.pmesh.mapSolutionToGrid(z)


# Make some plots


fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z), np.amax(Z), 30)
cp = plt.contourf(X, Y, Z, levels=levels)
plt.contour(X, Y, Z, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	

'''
ax.add_patch(
    patches.Rectangle(
        (0.0, 0.0),
        L,
        W,
        fill=False,      # remove background
        linewidth=2
    )
)
'''

plt.show()	




                                                         report.py                                                                                           0000755 0002735 0004704 00000010121 13032000747 013505  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python 

import pandas as pnd
import numpy  as np
import glob
import os
import matplotlib.pyplot as plt
import scipy.interpolate as sp_interp
import scipy.optimize    as optimize
from   FOMtoPandas import FOMtoPandas

# Create a new table, add a row for each split
# Use data from FOM.dat in each split

df = FOMtoPandas()
print df.head()
df = df[df['converged']=='True']
df = df[df['vsp'] ==2.4923]
df.to_csv(path_or_buf='done.csv')
print df.head()

# Make some plots
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

lines   = ['-','--','-.',':','-o', '-D']
colors  = ['r','b','k','m']
symbols = ['-o','D','s','^','D']
markersize = 7
lw  = 3

# Plotting
set_plt_defaults()

cat_var = 'vsp'


for cat_index, cat_val in enumerate(np.sort(df[cat_var].unique())):
	cat_df = df[df[cat_var]==cat_val]
	cur_col = colors[cat_index]

	for index, cur_height in enumerate(np.sort(df['W'].unique())):
		sdf = (cat_df[cat_df['W']==cur_height]).sort(columns=['Vg'], ascending=True)
		cur_line = lines[index] + cur_col	
		#cur_line = symbols[index] + 'b'
		fignum = 1	
		
		# Charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['sw_charge']
		plt.semilogy(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg')
		plt.ylabel('ns')
		plt.grid(True)
		
		
	
		# Gate Capacitance
		plt.figure(fignum)
		fignum += 1
		vvals = sdf['Vg'].values
		qvals = sdf['sw_charge']
			# Interpolate onto sparse mesh for smoothing
		q_interpolator = sp_interp.interp1d(vvals, qvals, kind='cubic')
		v_i = np.linspace(vvals[0], vvals[-1], 10)
		q_i = q_interpolator(v_i)
		dV  = v_i[1]-v_i[0]
		cvals = 1e6*1.6e-19*np.gradient(q_i, dV)
			# Now interpolate back to original grid
		c_interpolator = sp_interp.interp1d(v_i, cvals, kind='cubic')
		cvals_final = c_interpolator(vvals)
		plt.plot(vvals, cvals_final, cur_line, linewidth=lw)
		plt.xlabel('Vg')
		plt.ylabel('Gate Capacitance')
		plt.grid(True)
		
		
		
		# Gate capacitance with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg'].values
		dV    = xvals[1]-xvals[0]
		qvals = sdf['sw_charge']
		cvals = 1e6*1.6e-19*np.gradient(qvals, dV)
			# Interpolant for capacitance
		c_interp = sp_interp.interp1d(xvals, cvals)
		rhs = lambda x : c_interp(x)-1.0
			# Now find shift require to put cval=1.5 to Vg=0.0
		v0 = 0.0
		vf= optimize.fsolve(rhs, v0, xtol=1e-5)
		print 'Got vf=', vf
		yvals = cvals_final
		plt.plot(vvals-vf, cvals_final, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('Gate Capacitance')
		plt.grid(True)
		
		# Charge with Vt shift
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['sw_charge']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('ns')
		plt.grid(True)
		
		
		# Delta2 Fraction vs Vg-Vt
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['d2']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('d2_frac')
		plt.grid(True)
		
		# Delta2 Fraction vs sidewall charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['sw_charge']
		yvals = sdf['d2']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg')
		plt.ylabel('d2_frac')
		plt.grid(True)
		plt.xlim(xmin=1e11, xmax=1e14)
	
	
		# Velocity vs. Vg-Vt
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['Vg']
		yvals = sdf['vel']
		plt.plot(xvals-vf, yvals, cur_line, linewidth=lw)
		plt.xlabel('Vg-Vt')
		plt.ylabel('velocity')
		plt.grid(True)
	
		# Velocity vs. sidewall charge
		plt.figure(fignum)
		fignum += 1
		xvals = sdf['sw_charge']
		yvals = sdf['vel']
		plt.semilogx(xvals, yvals, cur_line, linewidth=lw)
		plt.xlabel('ns')
		plt.ylabel('velocity')
		plt.grid(True)
		plt.xlim(xmin=1e11, xmax=1e14)
	
plt.show()
	
	
	
	


		
	
                                                                                                                                                                                                                                                                                                                                                                                                                                               runfile.py                                                                                          0000755 0002735 0004704 00000007701 13031267517 013661  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import valley as valley
import SP2D as sp2d
import SCF  as SCF
import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys


import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib import colors, ticker, cm

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=300)


#+------------------------------------------------------------+
#|
#| Create Structure
#|
#+------------------------------------------------------------+

L = 5.0
W = 5.0
IL = 0.5

# grid
nx = 50
ny = 50
nh = 20
nv = 15
no = 5

#
# Finite-Difference Mesh: Used by Schroedinger
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.2
EOT = 0.85
vsp = 2.0
hsp = IL + (EOT-IL)*Eps_HiK/3.9

dev_type='pgaa'

if dev_type=='pgaa':
	pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
								vsp, hsp, nh, nv, no, 
								IL=IL, Eps_HiK=Eps_HiK)
elif dev_type=='fin':
	pgaa_zoner = zoner.finFETZoner(L, W, fdmesh, 
								vsp, hsp, nh, nv, no, 
								IL=IL, Eps_HiK=Eps_HiK)
elif dev_type=='nw':							
	pgaa_zoner = zoner.NWZoner(L, W, fdmesh, 
								vsp, hsp, nh, nv, no, 
								IL=IL, Eps_HiK=Eps_HiK)
else:
	print 'Unknown device type: ', dev_type
	exit()

#
# Finite-Volume Mesh: used by Poisson
mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			

# Assign the contact, call it 'Gate')			
mesh.assignContact( pgaa_zoner.contact_func, 'Gate' )

#+---------------------------------------------------------+
#|
#| Poisson-related Objects
#|
#+---------------------------------------------------------+

# Create Poisson Object
poisson = pois.FVMPoisson(mesh)

# Extract Laplace matrix; fixed dueing simulation
pmat =  poisson.Lap_mat

# Set Boundary Values on contacts
Vg = 0.6
bc_val_dict = { 'Gate' : Vg }


#+---------------------------------------------------------+
#|
#| Schroedinger-related Objects
#|
#+---------------------------------------------------------+


# Create Delta4 Valleys
# 'Left' Valleys
mstar = np.array([[0.315, 0.478, 0],[0.478, 0.19, 0],[0, 0, 0.315]])
mult  = 2
Delta4L = valley.valley('Delta4L', L, W, mstar, mult)

# 'Right' Valleys
mstar = np.array([[0.315, -0.478, 0],[-0.478, 0.19, 0],[0, 0, 0.315]])
mult  = 2
Delta4R = valley.valley('Delta4R', L, W, mstar, mult)

# Create Delta2 Valley
mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
mult = 2
Delta2 = valley.valley('Delta2', L, W, mstar, mult)

# Create Schroedinger Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Delta4L, Delta4R, Delta2], fdmesh)

#+---------------------------------------------------------+
#|
#| SCF Solver: Self-Consistent Solution of QM and Poisson
#|
#+---------------------------------------------------------+

scf = SCF.SCF( sp, poisson )
Ef = 0.0
(iterations, converged) = scf.solve(
								bc_val_dict, 
								max_iter = 200, 
								tol      = 1e-4, 
								alpha    = 0.1,
								n_eigen  = 20,
								Ef       = Ef 
							)

print 'Converged=', converged, ' after ', iterations, ' iterations.'

			
# Compute charge and velocity

(vel, tot_charge, valley_charge) = sp.computeChargeAndVelocity(Ef)
sidewall_charge = 1e14 * tot_charge / (2.0*(W+vsp))
nw_charge       = 1e14 * tot_charge / (2.0*(W+vsp+L+hsp))
d2_occupancy = valley_charge[2]/(tot_charge)

# Save everything to 'FOM.dat'
fom = open('FOM.dat', 'w')
fom.write('L ' + str(L) + '\n')
fom.write('W ' + str(W) + '\n')
fom.write('vsp ' + str(vsp) + '\n')
fom.write('Eps ' + str(Eps_HiK) + '\n')
fom.write('Vg  ' + str(Vg)  + '\n')
fom.write('sw_charge ' + str(sidewall_charge) + '\n')
fom.write('nw_charge ' + str(nw_charge) + '\n')
fom.write('tot_charge ' + str(tot_charge) + '\n')
fom.write('d2 ' + str(d2_occupancy) + '\n')
fom.write('vel ' + str(vel) + '\n')
fom.write('converged ' + str(converged) + '\n')

fom.close()






                                                               run.py                                                                                              0000755 0002735 0004704 00000006264 13051372354 013022  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python 



# Run script

import CombinatorialLoop as cl
import SplitTable        as st
import numpy as np
import os
from subprocess import call
from parm_subs  import parm_subs

	
# Default values
default_dict = { 'L' : 5.0, 
			 'W' : 5.0, 
			 'vsp' : 2.4923, 
			 'Vg': 0.6, 
			 'Eps': 22.2,
			 'dev_type' : 'fin',
			 'lattice_thick'            : 0.5,
			 'period_thickness_ratio'   : 7.0,
			 'barrier'                  : 0.0,
			 'd2_shift'                 : 0.0
			  }
			 
# Make a combinatorial loop over a few parameters	
c_pgaa = cl.CombinatorialLoop( { 'W'   : [3.0, 5.0, 7.0],
                                 'Vg'  : np.linspace(-0.2, 0.8, 30),
                                 'd2_shift' : [0,  -0.05],
                                 'dev_type' : ['pgaa']
                                 })             
                                 
c_nw = cl.CombinatorialLoop( { 'W'        : [3.0, 4.0, 5.0, 7.0, 11.0],
                               'Vg'       : np.linspace(-0.2, 0.8, 30),
                               'dev_type' : ['nw']
                             })    
                                 
c_fin = cl.CombinatorialLoop( { 'W'        : [30.0],
                                'Vg'       : np.linspace(-0.2, 0.8, 30),
                                'dev_type' : ['fin']
                              })     
                               
split_table = st.SplitTable( [c_pgaa, c_fin], default_dict ) 
n_nodes  = 60
n_splits = 0
offset   = 0



# Create run directories       
for split, cur_split in enumerate(split_table):

	print cur_split
	
	# Make a run directiory
	split_dir = 'Split_' + str(split+offset)
	if not os.path.exists(split_dir):
		os.makedirs(split_dir)
    	
    # Put substituted version of template into run directory
	parm_subs('template.py', split_dir + '/' + 'sim.py', cur_split)
	
	n_splits += 1
	
	
	
# Create jobfiles
for cur_jobfile_index in range(0, n_nodes):
	jf_name = 'jobfile_' + str(cur_jobfile_index)
	jobf = open(jf_name, 'w')
	
	# Header
	jobf.write('# Job Name\n' ) 
	jobf.write('#BSUB -J pGAA_SP\n') 
	jobf.write('# merge stdout and stderror\n')
	jobf.write(
	'#BSUB -oo /project/SAS_xfer/data_hub/bobradovic/5nmCMOS/pGAA/pGAA_Template/OutLog.o\n')
	jobf.write(
	'#BSUB -eo /project/SAS_xfer/data_hub/bobradovic/5nmCMOS/pGAA/pGAA_Template/OutLog.e\n')
	jobf.write('#BSUB -q "advmpi.q advsmp.q"\n')
	jobf.write('#BSUB -n 16\n') 
	jobf.write('#BSUB -R "span[hosts=1]"\n\n')
	
	# Now loop over all split directories which belong to this node
	for cur_split in range(0, n_splits):
		if cur_split % n_nodes == cur_jobfile_index:
			jobf.write('cd Split_' + str(cur_split+offset) + '\n')
			jobf.write('if [ -e "FOM.dat" ]\n')
			jobf.write('then\n')
			jobf.write('echo "Split already finished." > job.stdout\n')
			jobf.write('cd ..\n')
			jobf.write('else\n')
			jobf.write(
			'/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python sim.py > job.stdout'
			+ '\n')
			jobf.write('cd ..' + '\n')
			jobf.write('fi\n')
	jobf.write('exit' + '\n')
			
	jobf.close()
	
# Now run all the jobs
for job in range(0, n_nodes):
	print 'Submitting job: ', job
	jobname = 'jobfile_' + str(job)
	os.system("bsub < " + jobname)
	


                                                                                                                                                                                                                                                                                                                                            scf1.py                                                                                             0000755 0002735 0004704 00000011170 13027257370 013045  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import valley as valley
import SP2D as sp2d
import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys

import scipy.sparse.linalg as la_sparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=300)


L = 5.0
W = 5.0
IL = 0.5

# grid
nx = 30
ny = 30
nh = 20
nv = 15
no = 5
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.0
EOT = 0.85
vsp = 2.0
hsp = IL + (EOT-IL)*Eps_HiK/3.9

pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)


mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			
mesh.assignContact( lambda x, y: x < -hsp + 1e-8,   'Left' )
mesh.assignContact( lambda x, y: x >  L+hsp-1e-8, 'Right')
					
#print(mesh)

#exit()

poisson = pois.FVMPoisson(mesh)

#print 'Lap mat:'
pmat =  poisson.Lap_mat
#print pmat.toarray()


# Now build rhs mat
bc_val_dict = { 'Left' : 0.7, 'Right' : 0.7 }
peak = 8e19
xc = L/2.0
yc = W/2.0
rhs = poisson.buildRHS(
	bc_val_dict, 
	charge_func=lambda x, y : peak * np.exp(-(x-xc)*(x-xc)/(L*L) -(y-yc)*(y-yc)/(W*W))
)
#print 'RHS: ', rhs
print 'RHS shape: ', rhs.shape
print 'Poisson shape: ', pmat.shape


# Solve the resulting system
potential_solution = la_sparse.spsolve(pmat, rhs)

(X, Y, Z) = poisson.pmesh.mapSolutionToGrid(potential_solution)


# Make some plots


fig = plt.figure()
ax = fig.add_subplot(121, aspect='equal')
levels = np.linspace(np.amin(Z), np.amax(Z), 30)
cp = plt.contourf(X, Y, Z, levels=levels)
plt.contour(X, Y, Z, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	


ax.add_patch(
    patches.Rectangle(
        (0.0, 0.0),
        L,
        W,
        fill=False,      # remove background
        linewidth=2
    )
)

ax.add_patch(
    patches.Rectangle(
        (-IL, -IL),
        L+2*IL,
        W+2*IL,
        fill=False,      # remove background
        linewidth=2
    )
)


#plt.show()	




#
# Schroedinger part
#
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

# Create Schroedinger-Poisson Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Delta4L, Delta4R, Delta2], fdmesh)


n_eigen = 20

# Now solve system
# Need to map the potential_solution to a CB function
poisson.pmesh.loadAccumulator(potential_solution)
cbfun = lambda x, y: -1.0*poisson.pmesh.nn_interp(x, y)
print 'Starting eigensolve ...'
sp.solveWithCB(cbfun, n_eigen=n_eigen)

Ef = 0.0

# Build charge density
val_dict = sp.getSolution(0)
(X,Y,Z) = val_dict['wfns']
ZCharge = sp.computeChargeDensity(Ef)
		
# Now plot the accumulated charge density
'''
plt.figure()
levels = np.linspace(np.amin(ZCharge), np.amax(ZCharge), 30)
cp = plt.contourf(X, Y, ZCharge, levels=levels)
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	
'''
ax2 = fig.add_subplot(122, aspect='equal')
levels_c = np.linspace(np.amin(ZCharge), np.amax(ZCharge), 30)
cp = plt.contourf(X, Y, ZCharge, levels=levels_c)
#plt.contour(X, Y, ZCharge, levels=levels_c, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')


#plt.show()

print 'Re-solving system:'

# Now re-solve Poisson with the new charge
csemi_array = np.reshape(ZCharge, nx*ny)
fdmesh.loadAccumulator(csemi_array)
rhs = poisson.buildRHS(
	bc_val_dict, 
	charge_grid=fdmesh
)

# Solve the resulting system
z = la_sparse.spsolve(pmat, rhs)
(X, Y, Z) = poisson.pmesh.mapSolutionToGrid(z)
print 'Done re-solving.'

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z), np.amax(Z), 30)
cp = plt.contourf(X, Y, Z, levels=levels)
plt.contour(X, Y, Z, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	


ax.add_patch(
    patches.Rectangle(
        (0.0, 0.0),
        L,
        W,
        fill=False,      # remove background
        linewidth=2
    )
)

ax.add_patch(
    patches.Rectangle(
        (-IL, -IL),
        L+2*IL,
        W+2*IL,
        fill=False,      # remove background
        linewidth=2
    )
)

plt.show()




                                                                                                                                                                                                                                                                                                                                                                                                        scf2.py                                                                                             0000755 0002735 0004704 00000016574 13030545431 013053  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import valley as valley
import SP2D as sp2d
import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys

import scipy.sparse.linalg as la_sparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib import colors, ticker, cm

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=300)


#+------------------------------------------------------------+
#|
#| Create Structure
#|
#+------------------------------------------------------------+

L = 5.0
W = 5.0
IL = 0.5

# grid
nx = 30
ny = 30
nh = 20
nv = 15
no = 5

#
# Finite-Difference Mesh: Used by Schroedinger
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.0
EOT = 0.85
vsp = 2.0
hsp = IL + (EOT-IL)*Eps_HiK/3.9

pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)
							
#pgaa_zoner = zoner.finFETZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)
							
#pgaa_zoner = zoner.NWZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)

#
# Finite-Volume Mesh: used by Poisson
mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			

# Assign the contact, call it 'Gate')			
mesh.assignContact( pgaa_zoner.contact_func, 'Gate' )

#+---------------------------------------------------------+
#|
#| Poisson-related Objects
#|
#+---------------------------------------------------------+

# Create Poisson Object
poisson = pois.FVMPoisson(mesh)

# Extract Laplace matrix; fixed dueing simulation
pmat =  poisson.Lap_mat

# Set Boundary Values on contacts
#bc_val_dict = { 'Left' : 0.6, 'Right' : 0.6 }
bc_val_dict = { 'Gate' : 0.6 }


#+---------------------------------------------------------+
#|
#| Schroedinger-related Objects
#|
#+---------------------------------------------------------+

n_eigen = 20  # Number of eigenvalues to use
Ef      = 0.0 # Assign Fermi level to 0.0

# Create Delta4 Valleys
# 'Left' Valleys
mstar = np.array([[0.315, 0.478, 0],[0.478, 0.19, 0],[0, 0, 0.315]])
mult  = 2
Delta4L = valley.valley('Delta4L', L, W, mstar, mult)

# 'Right' Valleys
mstar = np.array([[0.315, -0.478, 0],[-0.478, 0.19, 0],[0, 0, 0.315]])
mult  = 2
Delta4R = valley.valley('Delta4R', L, W, mstar, mult)

# Create Delta2 Valley
mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
mult = 2
Delta2 = valley.valley('Delta2', L, W, mstar, mult)

# Create Schroedinger-Poisson Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Delta4L, Delta4R, Delta2], fdmesh)

#+---------------------------------------------------------+
#|
#| Initial Solver: Poisson with initial guess for charge
#|                 Followed by QM-soltion with guessed potential
#|
#+---------------------------------------------------------+


# Initial guess for charge is a volume-inverted Gaussian
peak = 4e19
xc = L/2.0
yc = W/2.0
rhs_initial = poisson.buildRHS(
	bc_val_dict, 
	charge_func=lambda x, y : peak * np.exp(-(x-xc)*(x-xc)/(L*L) -(y-yc)*(y-yc)/(W*W))
)


# Solve the Poisson system
initial_potential     = la_sparse.spsolve(pmat, rhs_initial)
(X, Y, Z_initial_pot) = poisson.pmesh.mapSolutionToGrid(initial_potential)
prev_potential        = initial_potential           # Initial potential array

# Initial Schroedinger solution
poisson.pmesh.loadAccumulator(initial_potential)
cbfun = lambda x, y: -1.0*poisson.pmesh.nn_interp(x, y)
sp.solveWithCB(cbfun, n_eigen=n_eigen)

# Build QM charge density
wfn_dict          = sp.getSolution(0)
(Xfd,Yfd,Z_eigen) = wfn_dict['wfns']
ZCharge_init      = sp.computeChargeDensity(Ef)
initial_charge    = np.reshape(ZCharge_init, nx*ny) # Initial charge array

#--------------- Done with initial solution



#+------------------------------------------------------------+
#|
#| Self-consistency SP Loop
#|
#+------------------------------------------------------------+
prev_charge    = initial_charge
prev_potential = initial_potential 
max_iter = 200
tol      = 1e-4
alpha    = 0.1

for cur_iter in range(1, max_iter):
	# Prepare for Poisson solve
	# Load previous Schroediner solution into accumulator
	# Get estimated new potential using full previous charge
	# Then "damp" by taking only a fraction of
	# the potential as the update
	fdmesh.loadAccumulator(prev_charge)
	rhs = poisson.buildRHS(
		bc_val_dict, 
		charge_grid=fdmesh
	)
	pred_potential = la_sparse.spsolve(pmat, rhs)
	delta_pot =  pred_potential-prev_potential
	norm = np.sqrt((delta_pot*delta_pot).sum())
	if cur_iter % 5 == 0 or cur_iter == 1:
		print ' >>> Iter: ', cur_iter, ' Norm: ', norm
	
	# Now do Schroedinger solve with damped new potential
	new_potential  = alpha * pred_potential + (1.0-alpha)*prev_potential
	prev_potential = new_potential # Save potential for next iter
	
	poisson.pmesh.loadAccumulator(new_potential)
	cbfun = lambda x, y: -1.0*poisson.pmesh.nn_interp(x, y)
	sp.solveWithCB(cbfun, n_eigen=n_eigen)
	
	wfn_dict          = sp.getSolution(0)
	(Xfd,Yfd,Z_eigen) = wfn_dict['wfns']
	ZCharge           = sp.computeChargeDensity(Ef)
	new_charge        = np.reshape(ZCharge, nx*ny) # Initial charge array 
    
	prev_charge = new_charge
	
	if norm < tol:
		print 'Converged after ', cur_iter, ' iterations.'
		break


# Done with loop, plot results
fig1 = plt.figure(1)

# Potential sub-figure
(Xp, Yp, Z_final_pot) = poisson.pmesh.mapSolutionToGrid(new_potential)

ax = fig1.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z_final_pot), np.amax(Z_final_pot), 30)
cp = plt.contourf(Xp, Yp, Z_final_pot, levels=levels)
plt.contour(X, Y, Z_final_pot, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	

ax.add_patch(
    patches.Rectangle(
        (0.0, 0.0),
        L,
        W,
        fill=False,      # remove background
        linewidth=2
    )
)

ax.add_patch(
    patches.Rectangle(
        (-IL, -IL),
        L+2*IL,
        W+2*IL,
        fill=False,      # remove background
        linewidth=2
    )
)

# QM Charge sub-figure
fig2 = plt.figure(2)
x2 = fig2.add_subplot(111, aspect='equal')
levels_c = np.linspace(np.amin(ZCharge), np.amax(ZCharge), 30)
#levels_c = np.linspace(5e18, 1e20, 30)
cp_QM = plt.contourf(Xfd, Yfd, ZCharge, levels=levels_c)
plt.colorbar(cp_QM)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')

# Get DOS for each valley
(energy, ret_val) = sp.computeDOS()
D4_DOS = ret_val[0]  # Symmetric D4 valleys, get only one, multiply by 2
D2_DOS = ret_val[2]  # Single D2 valley
plt.figure(3)
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
	nrows = 4
	ncols = 3
	div = (slabs-(slabs % ncols))/ncols + 1
	if div > nrows:
		div = nrows
	print 'div: ', div
	fig, axarr = plt.subplots(div, ncols)
	
	for row in range(0,div):
		for col in range(0,ncols):
			cur_slab = row * ncols + col
			if cur_slab == slabs:
			    break
			plotZ = Z[:,:,cur_slab] * Z[:,:,cur_slab]
			
			levels = np.linspace(np.amin(plotZ), np.amax(plotZ), 30)
			axarr[row, col].contourf(X, Y, plotZ, levels=levels)
			axarr[row, col].set_aspect('equal')

plt.show()


                                                                                                                                    scf3.py                                                                                             0000755 0002735 0004704 00000014670 13034463167 013060  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import valley as valley
import SP2D as sp2d
import SCF  as SCF
import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys
from   plot_fun import set_plt_defaults


import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib import colors, ticker, cm

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=300)


#+------------------------------------------------------------+
#|
#| Create Structure
#|
#+------------------------------------------------------------+

L = 5.0
W = 20.0
IL = 0.5

dy_target = 5.0/30
ny = int(W/dy_target)
print 'ny: ', ny
#exit()

# grid
nx = 30
#ny = 50
nh = 20
nv = 15
no = 5

#
# Finite-Difference Mesh: Used by Schroedinger
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.2
EOT = 0.85
#vsp = 2.4923
vsp = 2.5
hsp = IL + (EOT-IL)*Eps_HiK/3.9
print 'hsp: ', hsp

#pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)
							
#pgaa_zoner = zoner.MLFETZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)
							
pgaa_zoner = zoner.finFETZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)
							
#pgaa_zoner = zoner.NWZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)

#
# Finite-Volume Mesh: used by Poisson
mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			

# Assign the contact, call it 'Gate')			
mesh.assignContact( pgaa_zoner.contact_func, 'Gate' )

#+---------------------------------------------------------+
#|
#| Poisson-related Objects
#|
#+---------------------------------------------------------+

# Create Poisson Object
poisson = pois.FVMPoisson(mesh)

# Extract Laplace matrix; fixed dueing simulation
pmat =  poisson.Lap_mat

# Set Boundary Values on contacts
bc_val_dict = { 'Gate' : 0.7 }


#+---------------------------------------------------------+
#|
#| Schroedinger-related Objects
#|
#+---------------------------------------------------------+


# Create Delta4 Valleys
# 'Left' Valleys
mstar = np.array([[0.315, 0.478, 0],[0.478, 0.19, 0],[0, 0, 0.55]])
mult  = 2
Delta4L = valley.valley('Delta4L', L, W, mstar, mult)

# 'Right' Valleys
mstar = np.array([[0.315, -0.478, 0],[-0.478, 0.19, 0],[0, 0, 0.55]])
mult  = 2
Delta4R = valley.valley('Delta4R', L, W, mstar, mult)

# Create Delta2 Valley
mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
mult = 2
Delta2 = valley.valley('Delta2', L, W, mstar, mult)

# Create Schroedinger Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Delta4L, Delta4R, Delta2], fdmesh)

#+---------------------------------------------------------+
#|
#| SCF Solver: Self-Consistent Solution of QM and Poisson
#|
#+---------------------------------------------------------+

scf = SCF.SCF( sp, poisson )
Ef = 0.0
(iterations, converged) = scf.solve(
								bc_val_dict, 
								max_iter = 200, 
								tol      = 1e-4, 
								alpha    = 0.1,
								n_eigen  = 20,
								Ef       = Ef 
							)

print 'Converged=', converged, ' after ', iterations, ' iterations.'


#set_plt_defaults()
plt.rcParams['figure.facecolor'] = 'w'
plt.rc('xtick', labelsize=18) 
plt.rc('ytick', labelsize=18) 
font = {'family' : 'Bitstream Vera Sans',
	    'weight' : 'bold',
	    'size'   : 16}

plt.rc('font', **font)
plt.rcParams['toolbar'] = 'None'
#plt.tick_params(which='major', length=15, width=2)
#plt.tick_params(which='minor', length=5, width=2)

plt.rcParams['axes.linewidth'] = 2 

# Done with loop, plot results
fig1 = plt.figure(1)

# Potential sub-figure
(Xp, Yp, Z_final_pot) = poisson.pmesh.mapSolutionToGrid(scf.final_potential)

ax = fig1.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z_final_pot), np.amax(Z_final_pot), 30)
cp = plt.contourf(Xp, Yp, Z_final_pot, levels=levels)
plt.contour(Xp, Yp, Z_final_pot, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	

for cur_rect in pgaa_zoner.plotting_rectangle_tuple():
	start_point = cur_rect[0]
	length      = cur_rect[1]
	height      = cur_rect[2]
	ax.add_patch(
		patches.Rectangle(
		    start_point,
		    length,
		    height,
		    fill=False,      # remove background
		    linewidth=2
		)
	)



# QM Charge sub-figure
fig2 = plt.figure(2)
x2 = fig2.add_subplot(111, aspect='equal')
wfn_dict          = sp.getSolution(0)
(Xfd,Yfd,Z_eigen) = wfn_dict['wfns']
ZCharge           = sp.computeChargeDensity(Ef)
levels_c = np.linspace(np.amin(ZCharge), np.amax(ZCharge), 30)
cp_QM = plt.contourf(Xfd, Yfd, ZCharge, levels=levels_c)
plt.colorbar(cp_QM)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')

# Get DOS for each valley
(energy, ret_val) = sp.computeDOS(e_min=-0.1, e_max=0.2)
D4_DOS = ret_val[0]  # Symmetric D4 valleys, get only one, multiply by 2
D2_DOS = ret_val[2]  # Single D2 valley
plt.figure(3)
plt.plot(energy, 2*D4_DOS, '-r', linewidth=3)
plt.plot(energy, D2_DOS, '-b', linewidth=3)
plt.xlabel('Energy [eV]')
plt.ylabel('DOS')
plt.grid(True)
plt.xlim(xmin=-0.1, xmax=0.2)
plt.ylim(ymin=0.0, ymax=500)

# Plot eigenfunctions
# Separate plot for each valley
for cur_valley in range(0,3):
	val_dict = sp.getSolution(cur_valley)
	(X,Y,Z) = val_dict['wfns']
	(rows, cols, slabs) = Z.shape
	nrows = 4
	ncols = 3
	div = (slabs-(slabs % ncols))/ncols + 1
	if div > nrows:
		div = nrows
	print 'div: ', div
	fig, axarr = plt.subplots(div, ncols)
	
	for row in range(0,div):
		for col in range(0,ncols):
			cur_slab = row * ncols + col
			if cur_slab == slabs:
			    break
			plotZ = Z[:,:,cur_slab] * Z[:,:,cur_slab]
			
			levels = np.linspace(np.amin(plotZ), np.amax(plotZ), 30)
			axarr[row, col].contourf(X, Y, plotZ, levels=levels)
			axarr[row, col].set_aspect('equal')
			
			
# Compute charge and velocity

(vel, tot_charge, valley_charge) = sp.computeChargeAndVelocity(Ef)
sidewall_charge = 1e14 * tot_charge / (2.0*(W+vsp))
nw_charge       = 1e14 * tot_charge / (2.0*(W+vsp+L+hsp))
d2_occupancy = valley_charge[2]/(tot_charge)

print 'vel=', vel*1e-7, ' tot_charge=', tot_charge, \
      ' sidewall=', sidewall_charge, ' d2: ', d2_occupancy, \
      ' nw charge: ', nw_charge
#print 'valley charge: ', valley_charge



#sw_total_charge = 1e-14*sw_charge[0]*(2.0*W)
#print 'Standard SW charge: ', sw_charge, ' SW total charge: ', sw_total_charge
plt.show()


                                                                        scf4.py                                                                                             0000755 0002735 0004704 00000016005 13167751322 013053  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import ParabolicValley as pvalley
import SP2D as sp2d
import SCF  as SCF
import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys
from   plot_fun import set_plt_defaults
import SuperLattice as super_l


import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib import colors, ticker, cm

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=300)


#+------------------------------------------------------------+
#|
#| Create Structure
#|
#+------------------------------------------------------------+

# Fin dimensions
L = 5.0 
W_target = 3.0
IL = 0.5


# Superlattice
thick     = 1.0
dy_target = 0.0625*thick
period    = 2.5*thick
n_periods = np.floor(W_target/period)
W         = n_periods * period
offset    = -2.0*thick
energy    = 0.0
#energy    = 0.135
sl = super_l.SuperLattice(period, thick, energy, offset)


W = W_target


print 'Height=', W_target, ' period=', period

#dy_target = 5.0/30

ny = int(W/dy_target)
print 'ny: ', ny
#exit()

# grid
nx = 30
#ny = 50
nh = 20
nv = 15
no = 5

#
# Finite-Difference Mesh: Used by Schroedinger
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.2
EOT = 0.85
vsp = 2.4923
#vsp = 2.0 
hsp = IL + (EOT-IL)*Eps_HiK/3.9
print 'hsp: ', hsp

pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)
							
#pgaa_zoner = zoner.MLFETZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)
							
#pgaa_zoner = zoner.finFETZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)

'''							
pgaa_zoner = zoner.NWZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)
'''
#
# Finite-Volume Mesh: used by Poisson
mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			

# Assign the contact, call it 'Gate')			
mesh.assignContact( pgaa_zoner.contact_func, 'Gate' )

#+---------------------------------------------------------+
#|
#| Poisson-related Objects
#|
#+---------------------------------------------------------+

# Create Poisson Object
poisson = pois.FVMPoisson(mesh)

# Extract Laplace matrix; fixed dueing simulation
pmat =  poisson.Lap_mat

# Set Boundary Values on contacts
bc_val_dict = { 'Gate' : 1.0 }


#+---------------------------------------------------------+
#|
#| Schroedinger-related Objects
#|
#+---------------------------------------------------------+


# Create Delta4 Valleys
# 'Left' Valleys
mstar = np.array([[0.315, 0.478, 0],[0.478, 0.19, 0],[0, 0, 0.55]])
mult  = 2
Delta4L = pvalley.ParabolicValley('Delta4L', L, W, mstar, mult)

# 'Right' Valleys
mstar = np.array([[0.315, -0.478, 0],[-0.478, 0.19, 0],[0, 0, 0.55]])
mult  = 2
Delta4R = pvalley.ParabolicValley('Delta4R', L, W, mstar, mult)

# Create Delta2 Valley
# Also add energy offset due to strain
mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
mult = 2
shift = -0.05
Delta2 = pvalley.ParabolicValley('Delta2', L, W, mstar, mult, shift=shift)

# Create Schroedinger Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Delta4L, Delta4R, Delta2], fdmesh)

#+---------------------------------------------------------+
#|
#| SCF Solver: Self-Consistent Solution of QM and Poisson
#|
#+---------------------------------------------------------+

scf = SCF.SCF( sp, poisson, sl )
Ef = 0.0
(iterations, converged) = scf.solve(
								bc_val_dict, 
								max_iter = 200, 
								tol      = 2e-3, 
								alpha    = 0.1,
								n_eigen  = 20,
								Ef       = Ef 
							)

print 'Converged=', converged, ' after ', iterations, ' iterations.'


#set_plt_defaults()
plt.rcParams['figure.facecolor'] = 'w'
plt.rc('xtick', labelsize=18) 
plt.rc('ytick', labelsize=18) 
font = {'family' : 'Bitstream Vera Sans',
	    'weight' : 'bold',
	    'size'   : 16}

plt.rc('font', **font)
plt.rcParams['toolbar'] = 'None'
#plt.tick_params(which='major', length=15, width=2)
#plt.tick_params(which='minor', length=5, width=2)

plt.rcParams['axes.linewidth'] = 2 

# Done with loop, plot results
fig1 = plt.figure(1)

# Potential sub-figure
(Xp, Yp, Z_final_pot) = poisson.pmesh.mapSolutionToGrid(scf.final_potential)

ax = fig1.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z_final_pot), np.amax(Z_final_pot), 30)
cp = plt.contourf(Xp, Yp, Z_final_pot, levels=levels)
plt.contour(Xp, Yp, Z_final_pot, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	

for cur_rect in pgaa_zoner.plotting_rectangle_tuple():
	start_point = cur_rect[0]
	length      = cur_rect[1]
	height      = cur_rect[2]
	ax.add_patch(
		patches.Rectangle(
		    start_point,
		    length,
		    height,
		    fill=False,      # remove background
		    linewidth=2
		)
	)



# QM Charge sub-figure
fig2 = plt.figure(2)
x2 = fig2.add_subplot(111, aspect='equal')
wfn_dict          = sp.getSolution(0)
(Xfd,Yfd,Z_eigen) = wfn_dict['wfns']
ZCharge           = sp.computeChargeDensity(Ef)
levels_c = np.linspace(np.amin(ZCharge), np.amax(ZCharge), 30)
cp_QM = plt.contourf(Xfd, Yfd, ZCharge, levels=levels_c)
plt.colorbar(cp_QM)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')

# Get DOS for each valley
(energy, ret_val) = sp.computeDOS(e_min=-0.1, e_max=0.2)
D4_DOS = ret_val[0]  # Symmetric D4 valleys, get only one, multiply by 2
D2_DOS = ret_val[2]  # Single D2 valley
plt.figure(3)
plt.plot(energy, 2*D4_DOS, '-r', linewidth=3)
plt.plot(energy, D2_DOS, '-b', linewidth=3)
plt.xlabel('Energy [eV]')
plt.ylabel('DOS')
plt.grid(True)
plt.xlim(xmin=-0.1, xmax=0.2)
plt.ylim(ymin=0.0, ymax=500)

# Plot eigenfunctions
# Separate plot for each valley


for cur_valley in range(0,3):
	val_dict = sp.getSolution(cur_valley)
	(X,Y,Z) = val_dict['wfns']
	(rows, cols, slabs) = Z.shape
	nrows = 4
	ncols = 3
	div = (slabs-(slabs % ncols))/ncols + 1
	if div > nrows:
		div = nrows
	print 'div: ', div
	fig, axarr = plt.subplots(div, ncols)
	
	for row in range(0,div):
		for col in range(0,ncols):
			cur_slab = row * ncols + col
			if cur_slab == slabs:
			    break
			plotZ = Z[:,:,cur_slab] * Z[:,:,cur_slab]
			
			levels = np.linspace(np.amin(plotZ), np.amax(plotZ), 30)
			axarr[row, col].contourf(X, Y, plotZ, levels=levels)
			axarr[row, col].set_aspect('equal')
			axarr[row, col].xaxis.set_visible(False)
			axarr[row, col].yaxis.set_visible(False)
			
			
# Compute charge and velocity

(vel, tot_charge, valley_charge) = sp.computeChargeAndVelocity(Ef)
sidewall_charge = 1e14 * tot_charge / (2.0*(W+vsp))
nw_charge       = 1e14 * tot_charge / (2.0*(W+vsp+L+hsp))
d2_occupancy = valley_charge[2]/(tot_charge)

print 'vel=', vel*1e-7, ' tot_charge=', tot_charge, \
      ' sidewall=', sidewall_charge, ' d2: ', d2_occupancy, \
      ' nw charge: ', nw_charge
#print 'valley charge: ', valley_charge



#sw_total_charge = 1e-14*sw_charge[0]*(2.0*W)
#print 'Standard SW charge: ', sw_charge, ' SW total charge: ', sw_total_charge
plt.show()


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           SCF.py                                                                                              0000644 0002735 0004704 00000010446 13034464733 012627  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               
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
	
		print '+--------------+------------------------+'
		print '|    Iter      |    Potential Norm      |'
		print '+--------------+------------------------+'
	
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
				print '| {0:5d}        |    {1:8.4f}            |'.format(cur_iter, norm)
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
		
		print '+--------------+------------------------+'
		
		# Save potential from final iteration of Poisson		
		self.final_potential = new_potential
				
		return (cur_iter, converged)
                                                                                                                                                                                                                          Schroedinger.py                                                                                     0000644 0002735 0004704 00000011120 13045205607 014611  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               
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
                                                                                                                                                                                                                                                                                                                                                                                                                                                sim.py                                                                                              0000644 0002735 0004704 00000007701 13031303111 012757  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import valley as valley
import SP2D as sp2d
import SCF  as SCF
import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys


import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib import colors, ticker, cm

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=300)


#+------------------------------------------------------------+
#|
#| Create Structure
#|
#+------------------------------------------------------------+

L = 5.0
W = 9.0
IL = 0.5

# grid
nx = 50
ny = 50
nh = 20
nv = 15
no = 5

#
# Finite-Difference Mesh: Used by Schroedinger
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.2
EOT = 0.85
vsp = 2.0
hsp = IL + (EOT-IL)*Eps_HiK/3.9

dev_type='pgaa'

if dev_type=='pgaa':
	pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
								vsp, hsp, nh, nv, no, 
								IL=IL, Eps_HiK=Eps_HiK)
elif dev_type=='fin':
	pgaa_zoner = zoner.finFETZoner(L, W, fdmesh, 
								vsp, hsp, nh, nv, no, 
								IL=IL, Eps_HiK=Eps_HiK)
elif dev_type=='nw':							
	pgaa_zoner = zoner.NWZoner(L, W, fdmesh, 
								vsp, hsp, nh, nv, no, 
								IL=IL, Eps_HiK=Eps_HiK)
else:
	print 'Unknown device type: ', dev_type
	exit()

#
# Finite-Volume Mesh: used by Poisson
mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			

# Assign the contact, call it 'Gate')			
mesh.assignContact( pgaa_zoner.contact_func, 'Gate' )

#+---------------------------------------------------------+
#|
#| Poisson-related Objects
#|
#+---------------------------------------------------------+

# Create Poisson Object
poisson = pois.FVMPoisson(mesh)

# Extract Laplace matrix; fixed dueing simulation
pmat =  poisson.Lap_mat

# Set Boundary Values on contacts
Vg = 0.7
bc_val_dict = { 'Gate' : Vg }


#+---------------------------------------------------------+
#|
#| Schroedinger-related Objects
#|
#+---------------------------------------------------------+


# Create Delta4 Valleys
# 'Left' Valleys
mstar = np.array([[0.315, 0.478, 0],[0.478, 0.19, 0],[0, 0, 0.315]])
mult  = 2
Delta4L = valley.valley('Delta4L', L, W, mstar, mult)

# 'Right' Valleys
mstar = np.array([[0.315, -0.478, 0],[-0.478, 0.19, 0],[0, 0, 0.315]])
mult  = 2
Delta4R = valley.valley('Delta4R', L, W, mstar, mult)

# Create Delta2 Valley
mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
mult = 2
Delta2 = valley.valley('Delta2', L, W, mstar, mult)

# Create Schroedinger Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Delta4L, Delta4R, Delta2], fdmesh)

#+---------------------------------------------------------+
#|
#| SCF Solver: Self-Consistent Solution of QM and Poisson
#|
#+---------------------------------------------------------+

scf = SCF.SCF( sp, poisson )
Ef = 0.0
(iterations, converged) = scf.solve(
								bc_val_dict, 
								max_iter = 200, 
								tol      = 1e-4, 
								alpha    = 0.1,
								n_eigen  = 20,
								Ef       = Ef 
							)

print 'Converged=', converged, ' after ', iterations, ' iterations.'

			
# Compute charge and velocity

(vel, tot_charge, valley_charge) = sp.computeChargeAndVelocity(Ef)
sidewall_charge = 1e14 * tot_charge / (2.0*(W+vsp))
nw_charge       = 1e14 * tot_charge / (2.0*(W+vsp+L+hsp))
d2_occupancy = valley_charge[2]/(tot_charge)

# Save everything to 'FOM.dat'
fom = open('FOM.dat', 'w')
fom.write('L ' + str(L) + '\n')
fom.write('W ' + str(W) + '\n')
fom.write('vsp ' + str(vsp) + '\n')
fom.write('Eps ' + str(Eps_HiK) + '\n')
fom.write('Vg  ' + str(Vg)  + '\n')
fom.write('sw_charge ' + str(sidewall_charge) + '\n')
fom.write('nw_charge ' + str(nw_charge) + '\n')
fom.write('tot_charge ' + str(tot_charge) + '\n')
fom.write('d2 ' + str(d2_occupancy) + '\n')
fom.write('vel ' + str(vel) + '\n')
fom.write('converged ' + str(converged) + '\n')

fom.close()






                                                               SP2D.py                                                                                             0000644 0002735 0004704 00000017364 13045205623 012723  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               

import valley as valley
import FDMesh as fdmesh
import Schroedinger as sch

import numpy  as np
import scipy.interpolate as interp
import scipy.integrate   as integral
import scipy.optimize    as optimize
import phys



class SP2D:

    def __init__(self, valley_array, fdmesh):
        self.valley_array = valley_array
        self.fdmesh       = fdmesh
        self.Sch_array    = []
        self.nvalleys     = 0
        self.eigen_solution = []
        self.min_energy   = 0.0
        self.max_energy   = 0.0

        # Build Schroedinger solvers, one for each valley
        for cur_valley in self.valley_array:
            self.Sch_array.append( sch.Schroedinger(cur_valley, self.fdmesh) )
            self.nvalleys += 1
            # Each valley stores the "eigen_solution"
            # The eigen_solution is a dictionary
            # with the format: { 'evals' : [], 'wfns' : (X,Y,Z) }


    def getValley(self, index):
		return self.valley_array[index]


    def solveWithCB(self, CB_function, n_eigen=12):
        # Clear solution
        self.eigen_solution = []

		# Helper function for determining the top of the CB
        def maxCB(self, CB_function):
            L = self.fdmesh.L
            W = self.fdmesh.W
            return CB_function(L/2.0, W/2.0)
            	

        # Solve Schroedinger eigenproblem for each valley
        for cur_schroed in self.Sch_array:
            #print 'Starting numerical eigensolve ...'
            cb_max = maxCB(self, CB_function)
            (evals, (X,Y,Z)) = cur_schroed.solve(n_eigen, 
            									CB=CB_function,
            									sigma=cb_max)
            									
            self.eigen_solution.append(
                    { 'evals' : evals, 'wfns' : (X,Y,Z) }
            )
            #print 'Finished with numerical eigensolve ...'
            # Find energy range across all valleys
            all_evals = self.eigen_solution[0]['evals']
            for cur_eig in self.eigen_solution:
                all_evals = np.concatenate((all_evals, cur_eig['evals']))
            self.min_energy = np.min(all_evals)
            self.max_energy = np.max(all_evals)
            
            
            
            

    def getSolution(self, valley_index=-1):
        if valley_index == -1:
            return self.eigen_solution
        else:
            return self.eigen_solution[valley_index]

    def computeDOS(self, n_energy_points=300, e_min='none', e_max='none'):
        # Build energy axis
        if e_min == 'none':
        	e_min = self.min_energy-0.05
        	
        if e_max == 'none':
        	e_max = self.max_energy+0.1
        	
        energy = np.linspace(e_min, e_max, n_energy_points)
        ret_val = []
        for index, cur_valley in enumerate(self.valley_array):
            DOS = cur_valley.computeDOS(energy,
                                self.eigen_solution[index]['evals'])
            ret_val.append(DOS)

        return( (energy, ret_val) )
        
        

    def computeChargePlot(self, fmax, fermi_points=50, n_energy_points=300):
        # energy grid for integration
        energy = np.linspace(self.min_energy-0.05,
                             self.max_energy+0.1, n_energy_points)

        # Fermi energy grid
        fermi_vals  = np.linspace(energy[0]-3.0*phys.kT,
                                  fmax, fermi_points)

        charge_vals = np.zeros((self.nvalleys, fermi_points))
        tot_charge  = np.zeros(fermi_points)

        for i in range(0, fermi_points):
            if i % 5 == 0:
                print 'Working on Fermi level:', fermi_vals[i]
            for valley_index, cur_valley in enumerate(self.valley_array):
                # Compute valley charge
                charge_vals[valley_index,i] = cur_valley.computeQMCharge(
                                fermi_vals[i],
                                self.eigen_solution[valley_index]['evals'])
                # Compute total charge
                tot_charge[i] += charge_vals[valley_index, i]

        return ( fermi_vals, tot_charge, charge_vals )
        
    
    def computeSidewallCharge(self, Fermi):
		# energy grid for integration

		charge_vals = np.zeros(self.nvalleys)
		tot_charge  = 0.0

		for valley_index, cur_valley in enumerate(self.valley_array):
			# Compute valley charge
			charge_vals[valley_index] = cur_valley.computeSidewallCharge(
						    Fermi,
						    self.eigen_solution[valley_index]['evals'])
			# Compute total charge
			tot_charge += charge_vals[valley_index]

		return ( tot_charge, charge_vals )
		
	
    def computeChargeAndVelocity(self, Fermi):
        # Compute total charge
        tot_charge = 0.0
        valley_charge = []

        for valley_index, cur_valley in enumerate(self.valley_array):
            val_dict = self.getSolution(valley_index)
            valley_charge.append(0.0)

            print 'Working on valley ', valley_index, ' = ', valley_charge[valley_index]

            # Energies
            evals   = val_dict['evals']

            # Loop over all states in valley
            # Break from loop when state energy is Ef+5*kT
            for index, cur_Ec in enumerate(evals):

                # Have we included enough states?
                # Don't break until some charge has been added
                # This prevents the problem of getting absolute zero
                # charge under very low bias conditions
                #if (cur_Ec > Fermi + 20*phys.kT 
                #	and
                #	valley_charge[valley_index]  > 0.0):
	            #    break

                # compute the charge in each state
                state_charge = cur_valley.singleStateCharge(Fermi, cur_Ec)
                tot_charge        += state_charge
                valley_charge[valley_index] += state_charge
                print 'state charge=', state_charge, ' valley_charge=', valley_charge[valley_index], ' tot_charge=', tot_charge
			
        # Compute velocity integral, divide by total charge
        velocity = 1e-7*self.computeVelocity(Fermi) / tot_charge

        return (velocity, tot_charge, valley_charge)
		
		
		
    def computeChargeDensity(self, Ef):
        ZCharge = np.zeros((self.fdmesh.ny, self.fdmesh.nx))

        for valley_index, cur_valley in enumerate(self.valley_array):
            val_dict = self.getSolution(valley_index)
            # Wavefunctions - not yet normalized
            (X,Y,Z) = val_dict['wfns']
            
            # Energies
            evals   = val_dict['evals']

            # Loop over all states in valley
            # Break from loop when state energy is Ef+5*kT
            for index, cur_Ec in enumerate(evals):
            	# Have we included enough states?
                if cur_Ec > Ef + 5*phys.kT:
                	break

                    # compute the charge in each state
                    # Normalize wavefunctions to box
                norm_sqr =  self.fdmesh.dx * self.fdmesh.dy * np.sum(np.sum(
		        					Z[:,:,index] * Z[:,:,index] ) )
                state_charge = 1e21*cur_valley.singleStateCharge(Ef, cur_Ec) 
                
                ZCharge += state_charge * Z[:,:,index] * Z[:,:,index] / norm_sqr
                #print 'valley: ', valley_index, ' including state ', cur_Ec, \
                #          ' charge=', state_charge, ' norm_sqr=', norm_sqr
    	return ZCharge

        
    # Compute the total velocity integral across all modes and valleys
    # Must be normalized by total Q to get velocity    
    def computeVelocity(self, Fermi):
        # energy grid for integration

        vel_sum = 0.0
        
        for valley_index, cur_valley in enumerate(self.valley_array):
            # Compute valley velocity for each valley
            vel_sum += cur_valley.computeQMVelocity(
                            Fermi,
                            self.eigen_solution[valley_index]['evals'])

        return vel_sum
        
        
    
                                                                                                                                                                                                                                                                            sparse.py                                                                                           0000755 0002735 0004704 00000001141 13020130274 013465  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #! /opt/local/bin/python2.7

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
                                                                                                                                                                                                                                                                                                                                                                                                                               SplitTable.py                                                                                       0000644 0002735 0004704 00000002573 13032516727 014260  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               # Class for looping over 'split tables'
# The split table is presented as a collecttion of
# iterables
# The iterables are expected to produce a dictionary at every iteration



class SplitTable:
	def __init__(self, iterables, default_dict={}):
		self.iterables    = iterables
		self.default_dict = default_dict
		
	def __iter__(self):
		return self
		
	# Get complete split dictionary
	# No checks are made whether the data for the
	# split is available from the current iterable
	# This is checked in 'next', possibly switching
	# to the next iterable	
	def build_split(self):
	
		# Get new split values
		curVal =  self.iterables[-1].next()
		cur_split = {}
		
		# Accumulate default values and split values
		# into a single return dictionary
		for def_var in self.default_dict:
			cur_split[def_var] = self.default_dict[def_var]
		for cur_var in curVal:
			cur_split[cur_var] = curVal[cur_var]
		return cur_split	
		
		
		
	# Standard iterator 'next' function
	# Switches between iterables as the data
	# from each is exhausted	
	def next(self):
		try:
			# If this iterable is not empty, it will be the next split
			return self.build_split()
		
		except(StopIteration):
			# Done with this iterable, look for next one
			# and return first split of the new iterable
			self.iterables.pop()
			
			if len(self.iterables) > 0:
				return self.build_split()
			else:
				raise StopIteration
			
			
                                                                                                                                     sptest.py                                                                                           0000755 0002735 0004704 00000006333 13020130274 013522  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #! /opt/local/bin/python2.7



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
                                                                                                                                                                                                                                                                                                     SuperLattice.py                                                                                     0000644 0002735 0004704 00000002226 13034465146 014614  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #+---------------------------------------------------------
#|
#| Function for defining a superlattice CB offset
#| This needs to be added to the CB derived from
#| the electrostatic potential.
#|
#+---------------------------------------------------------

import matplotlib.pyplot as plt
import math as math

class SuperLattice:

	def __init__(self, period, thickness, energy, offset):
		self.period    = period
		self.thickness = thickness
		self.energy    = energy
		self.offset    = offset


	def evaluate(self, pos):
		# Determine the position of the evaluation point
		# within the supercell
		frac = ((pos-self.offset)/self.period - 
		               math.floor((pos-self.offset)/self.period))
		if frac < (self.period-self.thickness)/self.period:
			return 0.0
		else:
			return self.energy
		
		
if __name__ == "__main__":
	n_points  = 100
	period    = 5.0
	thickness = 1.0
	energy    = 1.0
	offset    = 3.0
	sl = SuperLattice(period, thickness, energy, offset)
	xvals = [ i*(10.0/(n_points-1)) for i in range(0, n_points) ]
	yvals = [ sl.evaluate(x) for x in xvals ]


	plt.figure(1)
	plt.plot(xvals, yvals, '-b')
	plt.ylim(ymin=-0.1, ymax = 1.5)
	plt.show()
	
	
                                                                                                                                                                                                                                                                                                                                                                          sweep_fin.py                                                                                        0000755 0002735 0004704 00000012674 13024275510 014173  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python



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
                                                                    sweep.py                                                                                            0000755 0002735 0004704 00000012674 13024273147 013343  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python



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
                                                                    template.py                                                                                         0000755 0002735 0004704 00000011620 13051357730 014022  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import ParabolicValley as pvalley
import SP2D as sp2d
import SCF  as SCF
import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys
import SuperLattice as super_l

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors, ticker, cm

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=300)


#+------------------------------------------------------------+
#|
#| Create Structure
#|
#+------------------------------------------------------------+

L        = @L@
W_target = @W@
IL = 0.5

# Superlattice
thick     = @lattice_thick@
energy    = @barrier@
#dy_target = 0.125*thick
dy_target = 0.5*thick
period    = @period_thickness_ratio@*thick
n_periods = np.floor(W_target/period)
if energy > 0.0:
	W         = n_periods * period
else:
	W = W_target
offset    = 0

d2_shift = @d2_shift@


sl = super_l.SuperLattice(period, thick, energy, offset)


# grid

nx = 30
ny = int(W/dy_target)
nh = 20
nv = 15
no = 5

#
# Finite-Difference Mesh: Used by Schroedinger
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = @Eps@
EOT = 0.85
vsp = @vsp@
hsp = IL + (EOT-IL)*Eps_HiK/3.9

dev_type='@dev_type@'

if dev_type=='pgaa':
	pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
								vsp, hsp, nh, nv, no, 
								IL=IL, Eps_HiK=Eps_HiK)
elif dev_type=='fin':
	pgaa_zoner = zoner.finFETZoner(L, W, fdmesh, 
								vsp, hsp, nh, nv, no, 
								IL=IL, Eps_HiK=Eps_HiK)
elif dev_type=='nw':
	vertical = min((vsp, hsp))							
	pgaa_zoner = zoner.NWZoner(L, W, fdmesh, 
								vertical, hsp, nh, nv, no, 
								IL=IL, Eps_HiK=Eps_HiK)

elif dev_type=='mlfet':							
	pgaa_zoner = zoner.MLFETZoner(L, W, fdmesh, 
								vsp, hsp, nh, nv, no, 
								IL=IL, Eps_HiK=Eps_HiK)								

else:
	print 'Unknown device type: ', dev_type
	exit()

#
# Finite-Volume Mesh: used by Poisson
mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			

# Assign the contact, call it 'Gate')			
mesh.assignContact( pgaa_zoner.contact_func, 'Gate' )

#+---------------------------------------------------------+
#|
#| Poisson-related Objects
#|
#+---------------------------------------------------------+

# Create Poisson Object
poisson = pois.FVMPoisson(mesh)

# Extract Laplace matrix; fixed dueing simulation
pmat =  poisson.Lap_mat

# Set Boundary Values on contacts
Vg = @Vg@
bc_val_dict = { 'Gate' : Vg }


#+---------------------------------------------------------+
#|
#| Schroedinger-related Objects
#|
#+---------------------------------------------------------+


# Create Delta4 Valleys
# 'Left' Valleys
mstar = np.array([[0.315, 0.478, 0],[0.478, 0.19, 0],[0, 0, 0.55]])
mult  = 2
Delta4L = pvalley.ParabolicValley('Delta4L', L, W, mstar, mult)

# 'Right' Valleys
mstar = np.array([[0.315, -0.478, 0],[-0.478, 0.19, 0],[0, 0, 0.55]])
mult  = 2
Delta4R = pvalley.ParabolicValley('Delta4R', L, W, mstar, mult)

# Create Delta2 Valley
mstar = np.array([[0.19, 0, 0],[0, 0.91, 0],[0, 0, 0.19]])
mult = 2
Delta2 = pvalley.ParabolicValley('Delta2', L, W, mstar, mult, shift=d2_shift)

# Create Schroedinger Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Delta4L, Delta4R, Delta2], fdmesh)

#+---------------------------------------------------------+
#|
#| SCF Solver: Self-Consistent Solution of QM and Poisson
#|
#+---------------------------------------------------------+

# Include a superlattice if lattice barrier > 0
if @barrier@ > 0.0:
	scf = SCF.SCF( sp, poisson, sl )
else:
	scf = SCF.SCF( sp, poisson )
	
Ef = 0.0
(iterations, converged) = scf.solve(
								bc_val_dict, 
								max_iter = 200, 
								tol      = 2e-3, 
								alpha    = 0.1,
								n_eigen  = 20,
								Ef       = Ef 
							)

print 'Converged=', converged, ' after ', iterations, ' iterations.'

			
# Compute charge and velocity

(vel, tot_charge, valley_charge) = sp.computeChargeAndVelocity(Ef)
sidewall_charge = 1e14 * tot_charge / (2.0*(W+vsp))
nw_charge       = 1e14 * tot_charge / (2.0*(W+vsp+L+hsp))
d2_occupancy = valley_charge[2]/(tot_charge)

# Save everything to 'FOM.dat'
fom = open('FOM.dat', 'w')
fom.write('L ' 			+ str(L) + '\n')
fom.write('W ' 			+ str(W) + '\n')
fom.write('period '     + str(period) + '\n')
fom.write('barrier '    + str(energy) + '\n')
fom.write('vsp ' 		+ str(vsp) + '\n')
fom.write('Eps ' 		+ str(Eps_HiK) + '\n')
fom.write('Vg  ' 		+ str(Vg)  + '\n')
fom.write('dev_type ' 	+ dev_type + '\n')
fom.write('sw_charge ' 	+ str(sidewall_charge) + '\n')
fom.write('nw_charge ' 	+ str(nw_charge) + '\n')
fom.write('tot_charge ' + str(tot_charge) + '\n')
fom.write('d2 ' 		+ str(d2_occupancy) + '\n')
fom.write('vel ' 		+ str(vel) + '\n')
fom.write('converged ' 	+ str(converged) + '\n')
fom.write('d2_shift ' + str(d2_shift) + '\n')

fom.close()






                                                                                                                test.py                                                                                             0000755 0002735 0004704 00000013402 13032155752 013165  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import valley as valley
import SP2D as sp2d
import SCF  as SCF
import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys


import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib import colors, ticker, cm

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=300)


#+------------------------------------------------------------+
#|
#| Create Structure
#|
#+------------------------------------------------------------+

L = 5.0
W = 9.0
IL = 0.5

# grid
nx = 30
ny = 30
nh = 20
nv = 15
no = 5

#
# Finite-Difference Mesh: Used by Schroedinger
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
Eps_HiK  = 22.2
EOT = 0.85
vsp = 2.4923
hsp = IL + (EOT-IL)*Eps_HiK/3.9
print 'hsp: ', hsp

#pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)
							
#pgaa_zoner = zoner.finFETZoner(L, W, fdmesh, 
#							vsp, hsp, nh, nv, no, 
#							IL=IL, Eps_HiK=Eps_HiK)
							
pgaa_zoner = zoner.NWZoner(L, W, fdmesh, 
							vsp, hsp, nh, nv, no, 
							IL=IL, Eps_HiK=Eps_HiK)

#
# Finite-Volume Mesh: used by Poisson
mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			

# Assign the contact, call it 'Gate')			
mesh.assignContact( pgaa_zoner.contact_func, 'Gate' )

#+---------------------------------------------------------+
#|
#| Poisson-related Objects
#|
#+---------------------------------------------------------+

# Create Poisson Object
poisson = pois.FVMPoisson(mesh)

# Extract Laplace matrix; fixed dueing simulation
pmat =  poisson.Lap_mat

# Set Boundary Values on contacts
bc_val_dict = { 'Gate' : 0.7 }


#+---------------------------------------------------------+
#|
#| Schroedinger-related Objects
#|
#+---------------------------------------------------------+


# Create Delta4 Valleys
# 'Left' Valleys
mstar = np.array([[0.315, 0.478, 0],[0.478, 0.19, 0],[0, 0, 0.55]])
mult  = 2
Delta4L = valley.valley('Delta4L', L, W, mstar, mult)

# 'Right' Valleys
mstar = np.array([[0.315, -0.478, 0],[-0.478, 0.19, 0],[0, 0, 0.55]])
mult  = 2
Delta4R = valley.valley('Delta4R', L, W, mstar, mult)

# Create Delta2 Valley
mstar = np.array([[0.19, 0, 0],[0, 0.92, 0],[0, 0, 0.19]])
mult = 2
Delta2 = valley.valley('Delta2', L, W, mstar, mult)

# Create Schroedinger Solver
# Include all valleys
print 'Constructing solver ...'
sp = sp2d.SP2D([Delta4R, Delta4L, Delta2], fdmesh)

#+---------------------------------------------------------+
#|
#| SCF Solver: Self-Consistent Solution of QM and Poisson
#|
#+---------------------------------------------------------+

scf = SCF.SCF( sp, poisson )
Ef = 0.0
(iterations, converged) = scf.solve(
								bc_val_dict, 
								max_iter = 200, 
								tol      = 1e-4, 
								alpha    = 0.1,
								n_eigen  = 20,
								Ef       = Ef 
							)

print 'Converged=', converged, ' after ', iterations, ' iterations.'




# Done with loop, plot results
fig1 = plt.figure(1)

# Potential sub-figure
(Xp, Yp, Z_final_pot) = poisson.pmesh.mapSolutionToGrid(scf.final_potential)

ax = fig1.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z_final_pot), np.amax(Z_final_pot), 30)
cp = plt.contourf(Xp, Yp, Z_final_pot, levels=levels)
plt.contour(Xp, Yp, Z_final_pot, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	

ax.add_patch(
    patches.Rectangle(
        (0.0, 0.0),
        L,
        W,
        fill=False,      # remove background
        linewidth=2
    )
)

ax.add_patch(
    patches.Rectangle(
        (-IL, -IL),
        L+2*IL,
        W+2*IL,
        fill=False,      # remove background
        linewidth=2
    )
)

# QM Charge sub-figure
fig2 = plt.figure(2)
x2 = fig2.add_subplot(111, aspect='equal')
wfn_dict          = sp.getSolution(0)
(Xfd,Yfd,Z_eigen) = wfn_dict['wfns']
ZCharge           = sp.computeChargeDensity(Ef)
levels_c = np.linspace(np.amin(ZCharge), np.amax(ZCharge), 30)
cp_QM = plt.contourf(Xfd, Yfd, ZCharge, levels=levels_c)
plt.colorbar(cp_QM)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')

# Get DOS for each valley
(energy, ret_val) = sp.computeDOS()
D4_DOS = ret_val[0]  # Symmetric D4 valleys, get only one, multiply by 2
D2_DOS = ret_val[2]  # Single D2 valley
plt.figure(3)
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
	nrows = 4
	ncols = 3
	div = (slabs-(slabs % ncols))/ncols + 1
	if div > nrows:
		div = nrows
	print 'div: ', div
	fig, axarr = plt.subplots(div, ncols)
	
	for row in range(0,div):
		for col in range(0,ncols):
			cur_slab = row * ncols + col
			if cur_slab == slabs:
			    break
			plotZ = Z[:,:,cur_slab] * Z[:,:,cur_slab]
			
			levels = np.linspace(np.amin(plotZ), np.amax(plotZ), 30)
			axarr[row, col].contourf(X, Y, plotZ, levels=levels)
			axarr[row, col].set_aspect('equal')
			
			
# Compute charge and velocity

(vel, tot_charge, valley_charge) = sp.computeChargeAndVelocity(Ef)
sidewall_charge = 1e14 * tot_charge / (2.0*(W+vsp))
nw_charge       = 1e14 * tot_charge / (2.0*(W+vsp+L+hsp))
d2_occupancy = valley_charge[2]/(tot_charge)

print 'vel=', vel*1e-7, ' tot_charge=', tot_charge, \
      ' sidewall=', sidewall_charge, ' d2: ', d2_occupancy, \
      ' nw charge: ', nw_charge
#print 'valley charge: ', valley_charge



#sw_total_charge = 1e-14*sw_charge[0]*(2.0*W)
#print 'Standard SW charge: ', sw_charge, ' SW total charge: ', sw_total_charge
plt.show()


                                                                                                                                                                                                                                                              valley.py                                                                                           0000644 0002735 0004704 00000010764 13045206267 013511  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               import numpy as np
import scipy.interpolate as interp
import scipy.integrate   as integral
import matplotlib.pyplot as plt
import phys as phys


#+------------------------------------------------------------------
#|
#| Base valley class
#|
#+------------------------------------------------------------------

class valley(object):
    """ 
        Base valley class
        This is an abstract class which contains generic 
        valley functionality, but nothing specific to the
        valley type. 
        
        The functions which need to be defined by sub-classes:
        1. __init__       : parameters for the valley types are different
        2. singleStateDOS : this depends on the shape of the valley
        3. singleStateVel : this depends on the shape of the valley
        
        Other functions, such as the computation of total charge etc.
        are provided in the generic functionality
    """

    def __init__(self, name, L, W, mstar, mult, shift):
        self.name = name
        self.L = L
        self.W = W
        self.mstar    = mstar
        self.mstarDOS = self.mstar[2,2]
        self.mult     = mult
        self.shift    = shift

    
    def estimateEnergyLevel(self, n, m):
        mxx    = self.mstar[0,0]
        myy    = self.mstar[1,1]
        prefac = phys.hbar * phys.hbar / (2.0 * phys.m0)
        kx     = phys.pi / (self.L*1e-9)
        ky     = phys.pi / (self.W*1e-9)
        Ex     = n * n * prefac * kx * kx / (mxx * phys.q)
        Ey     = m * m * prefac * ky * ky / (myy * phys.q)
        return Ex + Ey

    
    def singleStateDOS(self, E, Ec):
        """ 
        Abstract method 
        Compute DOS for a single 1-D sub-band
        NOTE: singularity at E==Ec, and invalid for
        energies below that.
        Must make sure only energies > Ec are computed
        Result in #/(J-nm)
        """
        raise NotImplementedError

    # Compute DOS, based on a spectrum of energy levels
    # Use the energy argument (array) as the x-axis grid
    def computeDOS(self, energy, energy_levels):
        n_points = energy.shape
        DOS    = np.zeros(n_points)
        for cur_energy in energy_levels:
            DOS[energy > cur_energy] += self.singleStateDOS(
                            energy[energy > cur_energy], cur_energy)
        return DOS * self.mult

    # Fermi-Dirac occupancy of a single energy level
    def occupancy(self, energy, Ef):
        return 1.0/(1+np.exp((energy-Ef)/phys.kT))

    # Compute charge in a single 1-D state
    # This is done by integrating the DOS function
    # Integrand is singula and requires one-sided integration
    # returns charge value
    # Units are #/nm
    def singleStateCharge(self, Ef, Ec, limit=300, epsrel=1e-5):
        integrand = lambda x : self.singleStateDOS(x, Ec) * \
                               self.occupancy(x, Ef)

        (charge, abserr) = integral.quad( integrand,
                                            Ec, np.inf,
                                            limit=limit,
                                            epsrel=epsrel)
        return charge
        
    
    
    
    def singleStateVel(self, Ef, Ec, limit=300, epsrel=1e-5):
       """ Abstract function """
       raise NotImplementedError    



    # Compute QM charge
    # Ef is Fermi level
    # evals is an array of eigen-energies
    # n_energy_points is the size of the energy grid
    # on which to evalue the DOS (basis for charge calc.)
    def computeSidewallCharge(self, Ef, evals,
                        limit=300,
                        epsrel=1e-5):

        # Integrate charge across all bound states,
        # accumulate charge
        charge = 0.0
        for cur_state in evals:
            charge += self.singleStateCharge(Ef, cur_state,
                                                limit=limit,
                                                epsrel=epsrel)

        # Divide total charge by sidewall length
        return  1e14 * charge / (2.0*self.W)
        #return  1e14 * charge / (2.0*self.L)
        
    
    def computeQMVelocity(self, Ef, evals,
                        limit=300,
                        epsrel=1e-5):

        # Integrate velocity*charge product across all bound states,
        # accumulate charge
        vel_sum    = 0.0
        
        for cur_state in evals:
            vel = self.singleStateVel(Ef, cur_state,
                                                limit=limit,
                                                epsrel=epsrel)
            vel_sum    += vel
            
       
        return  vel_sum
            wire.py                                                                                             0000755 0002735 0004704 00000003012 13020130274 013135  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #! /opt/local/bin/python2.7

import valley as valley
import numpy  as np
import scipy.interpolate as interp
import scipy.integrate   as integral
import phys
import matplotlib.pyplot as plt



# Generate DOS for a single valley
n_points = 400
Delta2 = valley.valley(4.0, 20.0, 0.92, 0.19, 2)
#Delta2 = valley(4.0, 40.0, 0.05, 0.05, 1)
energy = np.linspace(0, 0.5, n_points)
DOS    = Delta2.computeDOS(energy, 5, 50)

# Generate charge density for the valley
Ef = 0.1
occDOS = DOS * Delta2.occupancy( energy, Ef )
DOS_interp = interp.interp1d( energy, DOS )


fermi_points = 50
fermi_vals  = np.linspace(0, 0.2, fermi_points)
charge_vals = np.zeros(fermi_points)
for i in range(0, fermi_points):
    if i % 5 == 0:
        print 'Working on Fermi level:', fermi_vals[i]
    (charge, abserr) = integral.quad( lambda x :
                        phys.q * DOS_interp(x) * Delta2.occupancy(x, fermi_vals[i]),
                        0.0, fermi_vals[i] + 5.0 * phys.kT,
                        limit=200,
                        epsrel=1e-3)
    #print 'charge: ', charge, ' with relerr=', abserr/charge
    charge_vals[i] = 1e14 * charge / (2.0*Delta2.L + 2.0*Delta2.W)

#print 'fermi:', fermi_vals
#print 'charge:', charge_vals


plt.figure(1)
plt.plot(energy, DOS, '-b', linewidth=3)
plt.plot(energy, occDOS, '-r', linewidth=3)
plt.xlabel('Energy [eV]')
plt.ylabel('DOS')
plt.grid(True)

plt.figure(2)
plt.plot(fermi_vals, charge_vals, '-k', linewidth=3)
plt.xlabel('Fermi Level [eV]')
plt.ylabel('Electron Density [1/nm]')
plt.grid(True)

plt.show()
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      Zoner.py                                                                                            0000644 0002735 0004704 00000014263 13033516333 013303  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               import FDMesh as fdmesh
import numpy as np

class WireZoner:
	def __init__(self, L, W, fdmesh, vsp, hsp, nh, nv, no, IL=0.5, Eps_HiK=22.2):
		self.L 		= L
		self.W 		= W
		self.fdmesh = fdmesh
		self.vsp 	= vsp
		self.hsp 	= hsp
		self.nh     = nh
		self.nv     = nv
		self.no     = no
		self.IL     = IL
		
		# Property dictionary defines material properties
		self.prop_dict = { 'Silicon' : { 'eps' : 11.9 },
						   'Oxide'   : { 'eps' : 3.9  },
						   'HiK'     : { 'eps' : Eps_HiK },
						   'SiGe'    : { 'eps' : 16.0 } }
		
		self.buildMesh()
		
		
	def buildMesh(self):
		#------ Horizontal grid
		self.dx  = (self.hsp-self.IL) / (self.nh-1)
		self.dxo = self.IL / (self.no-1)

		# Additional vertical grids on left side
		xl = np.linspace(-self.hsp, -self.IL-self.dx, self.nh)
		
		
		# Additional vertical grids in the oxide
		xlo = np.linspace(-self.IL, -self.dxo, self.no)

		# Additional vertical grids on right side IL
		xro = np.linspace(self.fdmesh.L+self.dxo,
				         self.fdmesh.L+self.IL-self.dxo, self.no)
		#print 'xro: ', xro
				         
		xr  = np.linspace(self.L+self.IL, self.L+self.hsp, self.nh)	
		#print 'xr: ', xr

		# FDMesh horizontal grid
		xc = np.linspace(0.0, self.fdmesh.L, self.fdmesh.nx)


		self.xmesh = np.concatenate((xl, xlo, xc, xro, xr))
		#print 'xmesh: ', self.xmesh

		#------ Vertical grid
		self.dy  = (self.vsp-self.IL) / (self.nv-1)
		self.dyo = self.IL / (self.no-1)

		# Additional vertical grids on left side
		yl = np.linspace(-self.vsp, -self.IL-self.dy, self.nv)
		
		# Additional vertical grids in the oxide
		ylo = np.linspace(-self.IL, -self.dyo, self.no)

		# Additional vertical grids on right side IL
		yro = np.linspace(self.fdmesh.W+self.dyo,
				         self.fdmesh.W+self.IL-self.dyo, self.no)
				         
		yr  = np.linspace(self.W+self.IL, self.W+self.vsp, self.nv)	

		# FDMesh vertical grid
		yc = np.linspace(0.0, self.fdmesh.W, self.fdmesh.ny)

		self.ymesh = np.concatenate((yl, ylo, yc, yro, yr))
		#print 'ymesh: ', self.ymesh
	
		
	def getXMesh(self):
		return self.xmesh
		
	def getYMesh(self):
		return self.ymesh
		
	def zoner_func(self, x, y):
		print 'Calling abstract function zoner_func.'
		raise RuntimeException
		
	def contact_func(self, x, y):
		print 'Calling abstract function contact_func.'
		raise RuntimeException
		
	def plotting_rectangle_tuple(self):
		print 'Calling abstract function plotting_rectangle_list'
		raise RuntimeException
		
		
class pGAAZoner(WireZoner):

	# The zoner function has to return a material type string for each x, y
	# coordinate pair
	def zoner_func(self, x, y):
		if (     x>=0.0 
		     and x<=self.L 
		     and y>=0.0
		     and y<=self.W ):
		    return 'Silicon'
		    
		elif (   x>=0.0-self.IL 
		     and x<=self.L+self.IL 
		     and y>=0.0-self.IL
		     and y<=self.W+self.IL ):
			return 'Oxide'
		else:
			return 'HiK'
			
	# The contact function returns True if a given coordinate is contained
	# within the gate electrode, False otherwise		
	def contact_func(self, x, y):
		if (x<-self.hsp+1e-8) or (x>self.L+self.hsp-1e-8):
			return True
		else:
			return False
			
			
	def plotting_rectangle_tuple(self):
		return ( ((0.0,           0.0), self.L,             self.W),
				 ((-self.IL, -self.IL), self.L+2.0*self.IL, self.W+2.0*self.IL)
				)	
			
			
class finFETZoner(WireZoner):

	# The zoner function has to return a material type string for each x, y
	# coordinate pair
	def zoner_func(self, x, y):
		if (     x>=0.0 
		     and x<=self.L 
		     and y>=0.0
		     and y<=self.W ):
		    return 'Silicon'
		    
		elif (   x>=0.0-self.IL 
		     and x<=self.L+self.IL 
		     and y>=0.0-self.IL
		     and y<=self.W+self.IL
		     or  y<=0.0 ):
			return 'Oxide'
		else:
			return 'HiK'
			
	# The contact function returns True if a given coordinate is contained
	# within the gate electrode, False otherwise		
	def contact_func(self, x, y):
		if (   (x<-self.hsp+1e-8         and y > 0.0) 
		    or (x>self.L+self.hsp-1e-8   and y > 0.0)
		    or (y>self.W+self.vsp-1e-8)  and y > 0.0):
			return True
		else:
			return False
	
	def plotting_rectangle_tuple(self):
		return ( ( (0.0,             0.0), self.L,              self.W ),
				 ( (-self.IL,        0.0), self.L+2.0*self.IL,  self.W+self.IL ),
				 ( (-self.hsp, -self.vsp), self.L+2.0*self.hsp, self.vsp )
				)		

class NWZoner(WireZoner):

	# The zoner function has to return a material type string for each x, y
	# coordinate pair
	def zoner_func(self, x, y):
		if (     x>=0.0 
		     and x<=self.L 
		     and y>=0.0
		     and y<=self.W ):
		    return 'Silicon'
		    
		elif (   x>=0.0-self.IL 
		     and x<=self.L+self.IL 
		     and y>=0.0-self.IL
		     and y<=self.W+self.IL ):
			return 'Oxide'
		else:
			return 'HiK'
			
	# The contact function returns True if a given coordinate is contained
	# within the gate electrode, False otherwise		
	def contact_func(self, x, y):
		if (   (x<-self.hsp+1e-8     ) 
		    or (x>self.L+self.hsp-1e-8)
		    or (y>self.W+self.vsp-1e-8) 
		    or (y<-self.vsp+1e-8)     ):
			return True
		else:
			return False
			
	def plotting_rectangle_tuple(self):
		return ( ((0.0,           0.0), self.L,             self.W),
				 ((-self.IL, -self.IL), self.L+2.0*self.IL, self.W+2.0*self.IL)
				)		
			

class MLFETZoner(WireZoner):

	# The zoner function has to return a material type string for each x, y
	# coordinate pair
	def zoner_func(self, x, y):
		if (     x>=0.0 
		     and x<=self.L 
		     and y>=0.0
		     and y<=self.W ):
		    return 'Silicon'
		    
		elif (   x>=0.0
			 and x<=self.L
			 and (y<0.0 or y>self.W)):
			return 'SiGe'
		    
		elif (   x>=0.0-self.IL 
		     and x<=self.L+self.IL):
			return 'Oxide'
		else:
			return 'HiK'
			
	# The contact function returns True if a given coordinate is contained
	# within the gate electrode, False otherwise		
	def contact_func(self, x, y):
		if (x<-self.hsp+1e-8) or (x>self.L+self.hsp-1e-8):
			return True
		else:
			return False
			
	def plotting_rectangle_tuple(self):
		return ( ((0.0,           0.0), self.L,             self.W),
				 ((0.0,        self.W), self.L,           self.vsp),
				 ((0.0,     -self.vsp), self.L,           self.vsp),
				 ((-self.IL, -self.vsp),self.IL, self.W+2.0*self.vsp),
				 ((self.L,   -self.vsp),self.IL, self.W+2.0*self.vsp) 
				)		
		

                                                                                                                                                                                                                                                                                                                                             ztest.py                                                                                            0000755 0002735 0004704 00000003646 13026525704 013371  0                                                                                                    ustar   b.obradovic                     advlogic                                                                                                                                                                                                               #!/home/d.palle/localApps/Enthought/Canopy_64bit/User/bin/python

import Zoner   as zoner
import FDMesh  as fdmesh
import FVMMesh as FVM
import FVMPoisson  as pois
import numpy  as np
import phys

import scipy.sparse.linalg as la_sparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=2)
np.set_printoptions(linewidth=300)


L = 3.0
W = 3.0
IL = 0.5

# grid
nx = 30
ny = 30
nh = 20
nv = 20
fdmesh = fdmesh.FDMesh(L, W, nx, ny)

# Surounding oxide
vsp = 3.0
hsp = 3.0

pgaa_zoner = zoner.pGAAZoner(L, W, fdmesh, vsp, hsp, nh, nv, IL=IL)


mesh = FVM.FVMMesh( pgaa_zoner.getXMesh(), 
					pgaa_zoner.getYMesh(),
					pgaa_zoner.zoner_func, 
					pgaa_zoner.prop_dict ) 
			
mesh.assignContact( lambda x, y: x < -hsp + 1e-8,   'Left' )
mesh.assignContact( lambda x, y: x >  L+hsp-1e-8, 'Right')
					
#print(mesh)

#exit()

poisson = pois.FVMPoisson(mesh)

#print 'Lap mat:'
pmat =  poisson.Lap_mat
#print pmat.toarray()


# Now build rhs mat
bc_val_dict = { 'Left' : 0.5, 'Right' : 1.0 }
rhs = poisson.buildRHS(bc_val_dict)
#print 'RHS: ', rhs
print 'RHS shape: ', rhs.shape
print 'Poisson shape: ', pmat.shape


# Solve the resulting system
z = la_sparse.spsolve(pmat, rhs)
#print 'z:', z

(X, Y, Z) = poisson.pmesh.mapSolutionToGrid(z)


# Make some plots


fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
levels = np.linspace(np.amin(Z), np.amax(Z), 30)
cp = plt.contourf(X, Y, Z, levels=levels)
plt.contour(X, Y, Z, levels=levels, linewidths=0.5, colors='k')
plt.colorbar(cp)
plt.xlabel('X [nm]')
plt.ylabel('Y [nm]')	


ax.add_patch(
    patches.Rectangle(
        (0.0, 0.0),
        L,
        W,
        fill=False,      # remove background
        linewidth=2
    )
)

ax.add_patch(
    patches.Rectangle(
        (-IL, -IL),
        L+2*IL,
        W+2*IL,
        fill=False,      # remove background
        linewidth=2
    )
)


plt.show()	




                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          